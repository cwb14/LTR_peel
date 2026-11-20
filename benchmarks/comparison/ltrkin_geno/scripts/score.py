#!/usr/bin/env python3
"""Score the LTR callers against the pooled real insertion polymorphisms and write the
report tables.

Per dataset and tool: recall (found / partial / missed), exact-end recall, end-error bins, the ltrkin outcome of the
matched call (end-to-end and control runs), false insertion polymorphisms, recall by length and solo counts.
Confidence intervals are block bootstraps (N_BOOT draws, SEED, BLOCK-bp windows of the reference); the paired
difference is metric(tool) - metric(ltr_peel) over the same draws.

Coordinates are 0-based half-open. Blocks come from the ORIGINAL reference .fai: lg_config.DATASETS[ds]["ref"] + ".fai",
or, for tests, a `--fai DS=PATH` override. A tool whose run status is not ok gets status "DNF (<status>)" and empty
metric cells; its other inputs are not read. A run listed with `--running DS:TOOL` is still in progress: it gets
status "DNF (running)" and its resources.json is not read. A tool with status ok but a missing calls/truth/genotypes file is a
pipeline error and aborts.
"""
import argparse
import csv
import math
import os
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
# Appended, never prepended: the main benchmark's scripts dir has a queue.py that would shadow the stdlib module.
sys.path.append(os.path.join(os.environ["BENCH_ROOT"], "scripts"))
import json  # noqa: E402

import numpy as np  # noqa: E402

from calls_to_gff import gff_id  # noqa: E402
from lg_config import (ALL_TOOLS, BIN_BP, BLOCK, DATASETS, GQ_MIN, LEN_BINS, LG, N_BOOT, RO, SECONDARY, SEED,  # noqa: E402
                       TOOLS)
from lg_io import die, final_calls_path, read_calls  # noqa: E402
from pool import read_pool, read_truth  # noqa: E402
from score_elements import Index, match  # noqa: E402

CLASSES = ["missed", "partial", "not_genotypable", "F/E", "F/F", "F/N", "F/S", "no_call", "other"]
ALLELES = "FSEN"
LEN_LABELS = [f"<{LEN_BINS[0]}"] + [f"{a}-{b}" for a, b in zip(LEN_BINS, LEN_BINS[1:])] + [f">={LEN_BINS[-1]}"]


def warn(msg):
    print(f"warning: {msg}", file=sys.stderr)


def locus_status(calls, pool):
    """locus_id -> ("found", call idx) | ("partial", None) | ("missed", None).

    Matching (RO reciprocal overlap, one-to-one greedy) runs on match_start/match_end. Partial: some call overlaps
    the locus by >= 1 bp but none matched."""
    truth = [dict(chrom=p["chrom"], start=p["match_start"], end=p["match_end"]) for p in pool]
    _, t2c = match(calls, truth, RO)
    cidx = Index(calls)
    out = {}
    for ti, p in enumerate(pool):
        if ti in t2c:
            out[p["locus_id"]] = ("found", calls[t2c[ti]]["idx"])
        elif cidx.overlaps(p["chrom"], p["match_start"], p["match_end"]):
            out[p["locus_id"]] = ("partial", None)
        else:
            out[p["locus_id"]] = ("missed", None)
    return out


def end_errors(call, locus):
    """Signed end errors of a matched call against the true ends; None when the locus has no true ends."""
    ts, te = locus["true_start"], locus["true_end"]
    if ts is None or te is None:
        return None
    w = locus["window"] or 0
    lw = ts - w <= call["start"] <= ts
    rw = te <= call["end"] <= te + w
    return dict(left=call["start"] - ts, right=call["end"] - te, left_win=lw, right_win=rw, exact=lw and rw)


def boundary_summary(errs):
    """One row per end (left, right): n frac_exact frac_exact_window frac_le5bp frac_le20bp median_abs_bp.
    frac_exact: error == 0; frac_exact_window: end inside the microhomology window."""
    rows = []
    for end in ("left", "right"):
        n = len(errs)
        ab = [abs(e[end]) for e in errs]
        win = [e[f"{end}_win"] for e in errs]
        f = (lambda x: x / n) if n else (lambda x: None)
        rows.append(dict(end=end, n=n, frac_exact=f(sum(a == 0 for a in ab)), frac_exact_window=f(sum(win)),
                         frac_le5bp=f(sum(a <= BIN_BP[0] for a in ab)), frac_le20bp=f(sum(a <= BIN_BP[1] for a in ab)),
                         median_abs_bp=sorted(ab)[n // 2] if n else None))   # as approach 1 (score_elements.py): sorted |err|[n // 2]
    return rows


def genotype_class(gt):
    """F/E, F/F, F/N, F/S as is (alleles sorted by FSEN); ./. , . and empty -> no_call; anything else -> other."""
    if gt is None or gt.strip() in ("", ".", "./."):
        return "no_call"
    al = gt.strip().split("/")
    if len(al) != 2 or any(a not in ALLELES for a in al):
        return "other"
    g = "/".join(sorted(al, key=ALLELES.index))
    return g if g in ("F/E", "F/F", "F/N", "F/S") else "other"


def ltrkin_outcome(status, genos, gid):
    """(class, gq >= GQ_MIN) of a locus. not_genotypable = found, but the call is absent from the genotypes."""
    st = status[0]
    if st != "found":
        return st, False
    if gid not in genos:
        return "not_genotypable", False
    gt, gq = genos[gid]
    return genotype_class(gt), gq is not None and gq >= GQ_MIN


def false_e(genos, truth_rows):
    """Loop over the genotypes: calls with an E allele (background bg_* elements excluded). Confirmed by truth
    (n_e_confirmed), of which the confirmed truth state is F/Fx/Fn (n_false); not confirmed, including calls with no
    truth row (n_e_unconfirmed)."""
    tr = {t["element_id"]: t for t in truth_rows}
    nc = nf = nu = 0
    for eid, (gt, _) in genos.items():
        if eid.startswith("bg_") or "E" not in (gt or "").split("/"):
            continue
        t = tr.get(eid)
        if t is not None and t["truth_source"] == "confirmed":
            nc += 1
            nf += t["state"] in ("F", "Fx", "Fn")
        else:
            nu += 1
    return dict(n_e_confirmed=nc, n_false=nf, n_e_unconfirmed=nu)


def block_bootstrap(units, blocks, tools, ref_tool="ltr_peel", n=N_BOOT, seed=SEED):
    """Block bootstrap. units: [(block, {tool: (num, den)})]; blocks: all blocks of the reference. Each draw
    resamples len(blocks) blocks with replacement (one rng seeded once, the same draws for every tool); metric =
    sum(num) / sum(den). Per tool: lo hi (2.5/97.5 percentiles) and dlo dhi (same for metric - metric[ref_tool];
    nan when ref_tool is not among tools)."""
    bi = {b: i for i, b in enumerate(blocks)}
    nb = len(blocks)
    num = {t: np.zeros(nb) for t in tools}
    den = {t: np.zeros(nb) for t in tools}
    for b, d in units:
        if b not in bi:
            die(f"bootstrap unit in block {b} outside the reference blocks")
        for t in tools:
            a, c = d.get(t, (0, 0))
            num[t][bi[b]] += a
            den[t][bi[b]] += c
    idx = np.random.default_rng(seed).integers(0, nb, size=(n, nb))
    met = {}
    for t in tools:
        dd = den[t][idx].sum(1)
        with np.errstate(invalid="ignore", divide="ignore"):
            met[t] = np.where(dd > 0, num[t][idx].sum(1) / np.where(dd > 0, dd, 1), np.nan)

    def pct(x):
        if np.isnan(x).all():
            return float("nan"), float("nan")
        lo, hi = np.nanpercentile(x, [2.5, 97.5])
        return float(lo), float(hi)

    out = {}
    for t in tools:
        lo, hi = pct(met[t])
        dlo, dhi = pct(met[t] - met[ref_tool]) if ref_tool in met else (float("nan"), float("nan"))
        out[t] = dict(lo=lo, hi=hi, dlo=dlo, dhi=dhi)
    return out


# ---------------------------------------------------------------------------------------------------------------
def read_genos(path):
    """ltrkin genotypes TSV -> {element_id: (genotype, gq)}. Fails fast on a missing or 0-byte file and on a header
    without element_id, genotype, gq, checked from the header even when the file has no data rows."""
    if not os.path.isfile(path):
        die(f"file not found: {path}")
    if os.path.getsize(path) == 0:
        die(f"{path}: empty genotypes file (0 bytes); expected at least the header element_id, genotype, gq")
    out = {}
    with open(path, newline="") as f:
        rd = csv.DictReader(f, delimiter="\t")
        miss = [c for c in ("element_id", "genotype", "gq") if c not in (rd.fieldnames or [])]
        if miss:
            die(f"{path}: missing column(s) {miss}; expected element_id, genotype, gq")
        rows = list(rd)
    for r in rows:
        try:
            gq = int(float(r["gq"]))
        except (ValueError, KeyError):
            gq = None
        out[r["element_id"]] = (r["genotype"], gq)
    return out


def read_blocks(fai):
    if not os.path.isfile(fai):
        die(f"reference index not found: {fai}")
    blocks = []
    for line in open(fai):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 2:
            blocks += [(f[0], j) for j in range(-(-int(f[1]) // BLOCK))]
    return blocks


def run_status(root, ds, tool, running=()):
    src = "ltr_peel" if tool == "ltr_peel_hc" else tool
    if (ds, src) in running:
        return "running"
    p = f"{root}/runs/{ds}/{src}/resources.json"
    if not os.path.isfile(p):
        die(f"run status not found: {p}")
    try:
        st = json.load(open(p)).get("status")
    except (ValueError, AttributeError) as e:
        die(f"{p}: malformed resources.json ({e})")
    if not isinstance(st, str) or not st:
        die(f"{p}: resources.json lacks a 'status' string")
    return st


def need(path, what):
    if not os.path.isfile(path):
        die(f"{what} file missing for a tool with status ok: {path}")


def bin_of(length):
    return sum(length >= b for b in LEN_BINS)


def fmt(v):
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return ""
    if isinstance(v, float):
        return f"{v:.6g}"
    return str(v)


def write_table(path, cols, rows):
    with open(path, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(fmt(r.get(c)) for c in cols) + "\n")


def score_dataset(root, ds, fai, running=()):
    """Per-tool scoring of one dataset; returns dict of table rows."""
    pool_path = f"{root}/pool/{ds}/pool.tsv"
    pool = read_pool(pool_path)
    blocks = read_blocks(fai)
    N = len(pool)
    ub = [(p["chrom"], p["match_start"] // BLOCK) for p in pool]
    res = {t: dict(status=run_status(root, ds, t, running)) for t in ALL_TOOLS}
    ok = [t for t in ALL_TOOLS if res[t]["status"] == "ok"]
    # per-tool per-locus data
    for t in ok:
        cp = final_calls_path(f"{root}/calls/{ds}", t)
        need(cp, "calls")
        tp = f"{root}/truth/{ds}/{t}.tsv"
        need(tp, "truth")
        genos = {}
        for run, d in (("e2e", "ltrkin"), ("control", "control")):
            gp = f"{root}/{d}/{ds}/{t}.genotypes.tsv"
            need(gp, "genotypes")
            genos[run] = read_genos(gp)
        calls = read_calls(cp, t)
        by_idx = {c["idx"]: c for c in calls}
        st = locus_status(calls, pool)
        truth = read_truth(tp, t)
        r = res[t]
        r["status_rows"] = [st[p["locus_id"]] for p in pool]
        r["errs"] = [end_errors(by_idx[s[1]], p) if s[0] == "found" else None
                     for s, p in zip(r["status_rows"], pool)]
        r["genos"], r["truth"] = genos, truth
        r["outcome"] = {run: [ltrkin_outcome(s, genos[run], gff_id(t, s[1]) if s[0] == "found" else None)
                              for s in r["status_rows"]] for run in genos}
    tabs = dict(recall=[], boundaries=[], ltrkin_yield=[], false_e=[], by_length=[], solo=[])
    base = lambda t: dict(ds=ds, tool=t, role="secondary" if t in SECONDARY else "primary", status=res[t]["status"] if
                          res[t]["status"] == "ok" else f"DNF ({res[t]['status']})")
    if "ltr_peel" not in ok:
        warn(f"{ds}: ltr_peel did not finish; paired differences left empty")
    ci = {}

    def boot(fam, per_tool):
        """per_tool: {tool: [num per locus]}; denominators are all loci."""
        units = [(ub[i], {t: (per_tool[t][i], 1) for t in per_tool}) for i in range(N)]
        ci[fam] = block_bootstrap(units, blocks, list(per_tool)) if N and per_tool else {}

    boot("recall", {t: [int(s[0] == "found") for s in res[t]["status_rows"]] for t in ok})
    boot("exact", {t: [int(e is not None and e["exact"]) for e in res[t]["errs"]] for t in ok})
    for run in ("e2e", "control"):
        boot(f"yield_{run}", {t: [int(o[0] == "F/E") for o in res[t]["outcome"][run]] for t in ok})
    shares_ok = True
    binn = [sum(bin_of(p["match_end"] - p["match_start"]) == j for p in pool) for j in range(len(LEN_LABELS))]
    for t in ALL_TOOLS:
        b = base(t)
        if t not in ok or not N:
            for k in ("recall", "ltrkin_yield", "false_e", "boundaries", "by_length", "solo"):
                if k == "ltrkin_yield":
                    tabs[k] += [dict(b, run=run, n_loci=N or None) for run in ("e2e", "control")]
                elif k == "boundaries":
                    tabs[k] += [dict(b, end=e) for e in ("left", "right")]
                elif k == "by_length":
                    tabs[k] += [dict(b, len_bin=lb, n_loci=binn[j]) for j, lb in enumerate(LEN_LABELS)]
                elif k == "recall":
                    tabs[k].append(dict(b, n_loci=N or None))
                elif k == "false_e":
                    tabs[k].append(dict(b, run="e2e"))
                else:
                    tabs[k].append(dict(b))
            continue
        r = res[t]
        sts = [s[0] for s in r["status_rows"]]
        rec = dict(b, n_loci=N, found=sts.count("found"), partial=sts.count("partial"), missed=sts.count("missed"))
        rec["recall"] = rec["found"] / N
        c = ci["recall"][t]
        rec.update(recall_lo=c["lo"], recall_hi=c["hi"], recall_dlo=c["dlo"], recall_dhi=c["dhi"])
        ne = sum(e is not None and e["exact"] for e in r["errs"])
        rec.update(n_exact=ne, exact_end_recall=ne / N)
        c = ci["exact"][t]
        rec.update(exact_lo=c["lo"], exact_hi=c["hi"], exact_dlo=c["dlo"], exact_dhi=c["dhi"])
        tabs["recall"].append(rec)
        for row in boundary_summary([e for e in r["errs"] if e is not None]):
            tabs["boundaries"].append(dict(b, **row))
        for run in ("e2e", "control"):
            outs = r["outcome"][run]
            y = dict(b, run=run, n_loci=N)
            sh = {c_: sum(o[0] == c_ for o in outs) / N for c_ in CLASSES}
            if abs(sum(sh.values()) - 1) > 1e-9:
                shares_ok = False
                die(f"{ds} {t} {run}: breakdown shares sum to {sum(sh.values())}, not 1")
            y.update({f"share_{c_}": v for c_, v in sh.items()})
            y["yield"] = sh["F/E"]
            y["yield_gq20"] = sum(o[0] == "F/E" and o[1] for o in outs) / N
            c = ci[f"yield_{run}"][t]
            y.update(yield_lo=c["lo"], yield_hi=c["hi"], yield_dlo=c["dlo"], yield_dhi=c["dhi"])
            tabs["ltrkin_yield"].append(y)
        fe = false_e(r["genos"]["e2e"], r["truth"])
        fe["frac_false"] = fe["n_false"] / fe["n_e_confirmed"] if fe["n_e_confirmed"] else None
        tabs["false_e"].append(dict(b, run="e2e", **fe))
        for j, lb in enumerate(LEN_LABELS):
            idx = [i for i, p in enumerate(pool) if bin_of(p["match_end"] - p["match_start"]) == j]
            f = sum(sts[i] == "found" for i in idx)
            tabs["by_length"].append(dict(b, len_bin=lb, n_loci=len(idx), found=f,
                                          recall=f / len(idx) if idx else None))
        tabs["solo"].append(dict(b, n_solo_confirmed=sum(x["state"] == "S" and x["truth_source"] == "confirmed"
                                                         for x in r["truth"])))
    assert shares_ok
    return pool, tabs, res, ok, N


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--ds", nargs="+", default=["camara_hap2", "camara_hap1"], help="datasets to score")
    ap.add_argument("--root", default=str(LG), help="ltrkin_geno root holding calls/ pool/ truth/ ltrkin/ control/ "
                                                    "runs/ (default: production LG)")
    ap.add_argument("--out-dir", default=None, help="table directory (default: <root>/report/tables)")
    ap.add_argument("--fai", action="append", default=[], metavar="DS=PATH",
                    help="override the reference .fai of a dataset (tests); default: DATASETS[ds]['ref'] + '.fai'")
    ap.add_argument("--running", action="append", default=[], metavar="DS:TOOL",
                    help="a tool run still in progress when scoring (repeatable); its rows read \"DNF (running)\" "
                         "and its resources.json is not read")
    ap.add_argument("-v", "--verbose", action="store_true", help="per-tool progress on stderr")
    a = ap.parse_args(argv)
    fais = {}
    for x in a.fai:
        if "=" not in x:
            die(f"--fai expects DS=PATH, got {x!r}")
        k, v = x.split("=", 1)
        fais[k] = v
    running = set()
    for x in a.running:
        ds, sep, t = x.partition(":")
        if not sep or ":" in t or not ds or not t:
            die(f"--running expects DS:TOOL, got {x!r}")
        if ds not in a.ds:
            die(f"--running {x}: dataset {ds!r} is not among --ds {' '.join(a.ds)}")
        if t not in ALL_TOOLS:
            die(f"--running {x}: unknown tool {t!r} (known: {' '.join(ALL_TOOLS)})")
        running.add((ds, "ltr_peel" if t == "ltr_peel_hc" else t))
    out_dir = a.out_dir or f"{a.root}/report/tables"
    os.makedirs(out_dir, exist_ok=True)
    allt = {k: [] for k in ("recall", "boundaries", "ltrkin_yield", "false_e", "by_length", "solo")}
    pool_rows = []
    for ds in a.ds:
        if ds not in DATASETS:
            die(f"unknown dataset {ds!r}")
        fai = fais.get(ds, str(DATASETS[ds]["ref"]) + ".fai")
        pool, tabs, res, ok, N = score_dataset(a.root, ds, fai, running)
        pool_rows += [dict(p, ds=ds) for p in pool]
        for k in allt:
            allt[k] += tabs[k]
        yq = [r for r in tabs["ltrkin_yield"] if r["tool"] == "ltr_peel" and r["run"] == "e2e"]
        print(f"{ds}: pool={N} loci; tools ok={len(ok)}/{len(ALL_TOOLS)}; "
              f"ltr_peel yield={fmt(yq[0].get('yield')) if yq else ''}", file=sys.stderr)
        if a.verbose:
            for t in ALL_TOOLS:
                print(f"  {t}: {res[t]['status']}", file=sys.stderr)
    from pool import POOL_COLS
    write_table(f"{out_dir}/pool.tsv", ["ds"] + POOL_COLS, pool_rows)
    T = ["ds", "tool", "role", "status"]
    ci = lambda p: [f"{p}_lo", f"{p}_hi", f"{p}_dlo", f"{p}_dhi"]
    cols = dict(
        recall=T + ["n_loci", "found", "partial", "missed", "recall"] + ci("recall")
        + ["n_exact", "exact_end_recall"] + ci("exact"),
        boundaries=T + ["end", "n", "frac_exact", "frac_exact_window", "frac_le5bp", "frac_le20bp", "median_abs_bp"],
        ltrkin_yield=T + ["run", "n_loci"] + [f"share_{c}" for c in CLASSES] + ["yield", "yield_gq20"]
        + [f"yield_{s}" for s in ("lo", "hi", "dlo", "dhi")],
        false_e=T + ["run", "n_e_confirmed", "n_false", "n_e_unconfirmed", "frac_false"],
        by_length=T + ["len_bin", "n_loci", "found", "recall"],
        solo=T + ["n_solo_confirmed"])
    for k, c in cols.items():
        write_table(f"{out_dir}/{k}.tsv", c, allt[k])
    print(f"tables written to {out_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
