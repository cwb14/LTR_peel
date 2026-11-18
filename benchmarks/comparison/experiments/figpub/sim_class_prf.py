#!/usr/bin/env python3
"""Per-class precision / recall / F1 of LTR-RT callers in the nested simulations (S1_1..3, S2_1..3).

Truth: sims/truth/<S>.truth.tsv, ext scope (standard intact set + rebuilt hosts), as scripts/score_elements.py.
Matching: score_elements.match, one-to-one, reciprocal overlap >= 0.9 (reused, not reimplemented).
Truth classes (truth.tsv column `class`, set by printe_truth.classify):
  host      = class host or host_tenant (n_pieces > 1: split by at least one tenant, nested or not)
  tenant    = class tenant (nested, n_pieces == 1: inside another feature, not itself split)
  unnested  = class intact (not nested, not split)
Recall per class = matched truth elements of the class / truth elements of the class.
Precision per class: a matched call takes its truth element's class; an unmatched call takes the class of its
false-positive anatomy category (score_elements.anatomy, applied to that one call): partial_host and
partial_host_tenant -> host; partial_tenant -> tenant; partial_intact -> unnested; any other category (LTR_frag,
LTR_solo, gene, nonLTR_TE, nonLTR_frag, no_feature, ...) -> unnested (a standalone call).
Coordinates: calls TSV and truth TSV are both 0-based half-open (score_elements readers); no conversion here.

Checks (fatal): per genome x tool, summed over classes, n_calls, matched calls, n_truth and matched truth equal the
ext row (ro 0.9) of eval/sim/<S>/<tool>.metrics.tsv and the printed precision / recall strings are identical; per-class
matched truth equals the ext:<class> rows; per-call anatomy categories sum to eval/sim/<S>/<tool>.anatomy.tsv.
Outputs (--out): sim_class_prf.tsv, sim_class_prf_summary.tsv, sim_class_prf_merge_rule.tsv.
"""
import argparse
import csv
import math
import os
import statistics
import sys
from collections import Counter, defaultdict
from multiprocessing import Pool

sys.dont_write_bytecode = True
B = os.environ["BENCH_ROOT"]
sys.path.append(f"{B}/scripts")   # appended: scripts/queue.py would shadow the stdlib queue
from score_elements import Index, anatomy, match, read_bed3, read_calls, read_truth  # noqa: E402

TOOLS = ["ltrquest", "ltr_retriever", "ltrharvest", "ltr_finder", "look4ltrs", "hirepeat", "hite", "ltrpred",
         "dante_ltr"]
DIAG = ["ltrquest_hc"]            # developer diagnostic: computed, flagged role=diagnostic, never in figures
SETS = [f"S{s}_{r}" for s in (1, 2) for r in (1, 2, 3)]
CLASSES = ["host", "tenant", "unnested"]
RO = 0.9
TRUTH_CLASS = {"host": "host", "host_tenant": "host", "tenant": "tenant", "intact": "unnested"}
FP_CLASS = {"partial_host": "host", "partial_host_tenant": "host", "partial_tenant": "tenant",
            "partial_intact": "unnested"}   # every other anatomy category -> unnested
MERGE_DIFF, MERGE_MIN_TOOLS = 0.05, 2   # pre-set: "separate" if |F1 unnested - F1 tenant| >= 0.05 for >= 2 of 9 tools


def die(msg):
    sys.exit(f"error: {msg}")


def truth_class(t):
    c = TRUTH_CLASS.get(t["class"])
    if c is None:
        die(f"unknown truth class {t['class']!r} for {t.get('id')}")
    return c


def fp_class(cat):
    return FP_CLASS.get(cat, "unnested")


def call_categories(calls, unmatched, truth, feats, ro):
    """Anatomy category of each unmatched call: score_elements.anatomy on that call alone. anatomy only reads truth
    elements and features overlapping the call, so it gets just those (in their original order, which keeps its
    tie-breaking identical) instead of rebuilding indexes over the whole genome for every call."""
    tidx, fidx = Index(truth), Index(feats)
    out = {}
    for ci in unmatched:
        c = calls[ci]
        ts = [truth[k] for k in sorted(tidx.overlaps(c["chrom"], c["start"], c["end"]))]
        fs = [feats[k] for k in sorted(fidx.overlaps(c["chrom"], c["start"], c["end"]))]
        cnt = anatomy([c], [0], ts, fs, ro)
        if sum(cnt.values()) != 1:
            die("anatomy returned != 1 category for one call")
        out[ci] = next(iter(cnt))
    return out


def class_counts(calls, truth, feats, ro=RO):
    """Per class: n_truth, tp_truth, n_calls, tp_calls; plus the per-call anatomy categories and match sizes."""
    c2t, t2c = match(calls, truth, ro)
    k = {c: dict(n_truth=0, tp_truth=0, n_calls=0, tp_calls=0) for c in CLASSES}
    for ti, t in enumerate(truth):
        c = truth_class(t)
        k[c]["n_truth"] += 1
        k[c]["tp_truth"] += ti in t2c
    unmatched = [ci for ci in range(len(calls)) if ci not in c2t]
    cats = call_categories(calls, unmatched, truth, feats, ro)
    for ci in range(len(calls)):
        if ci in c2t:
            c = truth_class(truth[c2t[ci]])
            k[c]["tp_calls"] += 1
        else:
            c = fp_class(cats[ci])
        k[c]["n_calls"] += 1
    return k, cats, c2t, t2c, unmatched


def prf(k):
    r = k["tp_truth"] / k["n_truth"] if k["n_truth"] else float("nan")
    p = k["tp_calls"] / k["n_calls"] if k["n_calls"] else float("nan")
    f = 2 * p * r / (p + r) if math.isfinite(p) and math.isfinite(r) and p + r > 0 else float("nan")
    return p, r, f


def fmt(x):
    return "NA" if x is None or (isinstance(x, float) and not math.isfinite(x)) else (
        f"{x:.4f}" if isinstance(x, float) else str(x))


def read_metrics(path):
    rows = {}
    for r in csv.DictReader(open(path), delimiter="\t"):
        if r["ro"] == "0.9":
            rows[r["scope"]] = r
    return rows


def run_one(args):
    ds, tool = args
    cp = f"{B}/calls/{ds}/{tool}.tsv"
    if not os.path.isfile(cp):
        die(f"calls not found: {cp}")
    calls = read_calls(cp)
    truth = read_truth(f"{B}/sims/truth/{ds}.truth.tsv")
    feats = read_bed3(f"{B}/sims/truth/{ds}.features.bed")
    k, cats, c2t, t2c, unmatched = class_counts(calls, truth, feats)
    bad = []
    # check 1: summed over classes == ext metrics (ro 0.9), integers and printed strings
    m = read_metrics(f"{B}/eval/sim/{ds}/{tool}.metrics.tsv")
    ext = m["ext"]
    tot = {x: sum(k[c][x] for c in CLASSES) for x in ("n_truth", "tp_truth", "n_calls", "tp_calls")}
    if (tot["n_calls"], tot["n_truth"], tot["tp_calls"], tot["tp_truth"]) != \
            (int(ext["n_calls"]), int(ext["n_truth"]), int(ext["tp"]), int(ext["tp"])):
        bad.append(f"{ds} {tool}: totals {tot} vs ext n_calls {ext['n_calls']} n_truth {ext['n_truth']} tp {ext['tp']}")
    p = tot["tp_calls"] / tot["n_calls"] if tot["n_calls"] else float("nan")
    r = tot["tp_truth"] / tot["n_truth"]
    if f"{p:.4f}" != ext["precision"] or f"{r:.4f}" != ext["recall"]:
        bad.append(f"{ds} {tool}: P/R {p:.4f}/{r:.4f} vs ext {ext['precision']}/{ext['recall']}")
    # check 2: per-class matched truth == ext:<class> rows
    exp = {"host": int(m["ext:host"]["tp"]) + int(m["ext:host_tenant"]["tp"]), "tenant": int(m["ext:tenant"]["tp"]),
           "unnested": int(m["ext:intact"]["tp"])}
    for c in CLASSES:
        if k[c]["tp_truth"] != exp[c]:
            bad.append(f"{ds} {tool} {c}: tp_truth {k[c]['tp_truth']} vs metrics {exp[c]}")
    # check 3: per-call categories == anatomy on all unmatched calls == eval anatomy.tsv
    per_call = Counter(cats.values())
    whole = anatomy(calls, unmatched, truth, feats, RO)
    ev = Counter({r["category"]: int(r["n"]) for r in csv.DictReader(open(f"{B}/eval/sim/{ds}/{tool}.anatomy.tsv"),
                                                                       delimiter="\t")})
    if per_call != whole or per_call != ev:
        bad.append(f"{ds} {tool}: anatomy per call {dict(per_call)} vs whole {dict(whole)} vs eval {dict(ev)}")
    rows = []
    for c in CLASSES:
        pp, rr, ff = prf(k[c])
        rows.append(dict(set=ds[:2], seed=ds[3], tool=tool, cls=c, n_truth=k[c]["n_truth"], tp_truth=k[c]["tp_truth"],
                         recall=rr, n_calls=k[c]["n_calls"], tp_calls=k[c]["tp_calls"], precision=pp, f1=ff,
                         role="diagnostic" if tool in DIAG else "primary"))
    return rows, bad, (ds, tool, dict(per_call))


def summarize(rows):
    g = defaultdict(list)
    for r in rows:
        g[(r["set"], r["tool"], r["cls"])].append(r)
    out = []
    for (st, tool, c), rr in g.items():
        d = dict(set=st, tool=tool, cls=c)
        for x in ("recall", "precision", "f1"):
            v = [r[x] for r in rr if math.isfinite(r[x])]
            d[f"{x}_mean"] = statistics.mean(v) if v else float("nan")
            d[f"{x}_sd"] = statistics.stdev(v) if len(v) > 1 else float("nan")
            d[f"{x}_n_genomes"] = len(v)
        d["n_truth"] = sum(r["n_truth"] for r in rr)
        d["n_calls"] = sum(r["n_calls"] for r in rr)
        d["role"] = rr[0]["role"]
        out.append(d)
    return out


def merge_rule(summ):
    s = {(r["tool"], r["cls"]): r for r in summ if r["set"] == "S2"}
    rows = []
    for t in TOOLS:
        diff = s[(t, "unnested")]["f1_mean"] - s[(t, "tenant")]["f1_mean"]
        rows.append(dict(tool=t, f1_unnested=s[(t, "unnested")]["f1_mean"], f1_tenant=s[(t, "tenant")]["f1_mean"],
                         diff=diff, large=int(abs(diff) >= MERGE_DIFF)))
    n = sum(r["large"] for r in rows)
    return rows, n, "separate" if n >= MERGE_MIN_TOOLS else "merge"


def write(path, rows, cols):
    with open(path, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(fmt(r[c]) for c in cols) + "\n")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", default=os.path.dirname(os.path.abspath(__file__)), help="output dir (default: script dir)")
    ap.add_argument("-t", "--threads", type=int, default=4, help="worker processes (default 4)")
    ap.add_argument("-v", "--verbose", action="store_true", help="per-run anatomy categories on stderr")
    a = ap.parse_args(argv)
    jobs = [(ds, t) for ds in SETS for t in TOOLS + DIAG]
    print(f"sim_class_prf: {len(jobs)} genome x tool runs", file=sys.stderr)
    with Pool(a.threads) as pool:
        res = pool.map(run_one, jobs)
    bad = [b for _, bb, _ in res for b in bb]
    if bad:
        die("checks failed:\n  " + "\n  ".join(bad))
    print(f"checks passed: {len(jobs)} runs reproduce ext n_calls/n_truth/tp, precision and recall strings, "
          "ext:<class> tp and anatomy.tsv counts", file=sys.stderr)
    if a.verbose:
        for _, _, (ds, t, c) in res:
            print(ds, t, c, file=sys.stderr)
    rows = [r for rr, _, _ in res for r in rr]
    order = {t: i for i, t in enumerate(TOOLS + DIAG)}
    rows.sort(key=lambda r: (r["set"], r["seed"], order[r["tool"]], CLASSES.index(r["cls"])))
    write(f"{a.out}/sim_class_prf.tsv", rows, ["set", "seed", "tool", "cls", "n_truth", "tp_truth", "recall", "n_calls",
                                               "tp_calls", "precision", "f1", "role"])
    summ = summarize(rows)
    summ.sort(key=lambda r: (r["set"], order[r["tool"]], CLASSES.index(r["cls"])))
    write(f"{a.out}/sim_class_prf_summary.tsv", summ,
          ["set", "tool", "cls"] + [f"{x}_{y}" for x in ("recall", "precision", "f1") for y in ("mean", "sd", "n_genomes")]
          + ["n_truth", "n_calls", "role"])
    mr, n, verdict = merge_rule(summ)
    write(f"{a.out}/sim_class_prf_merge_rule.tsv", mr, ["tool", "f1_unnested", "f1_tenant", "diff", "large"])
    with open(f"{a.out}/sim_class_prf_merge_rule.tsv", "a") as f:
        f.write(f"# S2 rule (pre-set): separate if |diff| >= {MERGE_DIFF} for >= {MERGE_MIN_TOOLS} of {len(TOOLS)} "
                f"tools; {n} tools -> {verdict}\n")
    print(f"S2 merge rule: {n} of {len(TOOLS)} tools with |F1 unnested - F1 tenant| >= {MERGE_DIFF} -> {verdict}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
