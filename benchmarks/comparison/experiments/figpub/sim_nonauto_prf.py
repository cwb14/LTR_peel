#!/usr/bin/env python3
"""Coding / partial / protein-less recall, precision and F1 in the simulations (S1, S2), nonauto rule.

Why recomputed: experiments/nonauto/sim_callclass.tsv classes every call "proteinless" (all 60 genome x tool rows).
nonauto_sim.call_class looks the call's chromosome up in runs/<S>/dante_ltr/dante.gff3, whose sequence IDs are the
prep_genome.py names (c0001..), while the calls carry the original PrinTE names (chr1..); doms.get('chr1') is None,
so every call falls to "proteinless" and the row equals all calls (precision = std-scope precision). Here the DANTE
chromosome IDs are mapped back with data/<S>/names.tsv (new_id -> original_id) and the rest is the nonauto code:
  truth, element class, has_ltr filter: nonauto_sim.truth (std == 1 rows; class of the PrinTE library entry);
  recall per element class: elements whose library entry has an LTR pair (has_ltr), matched one-to-one at RO >= 0.9
  (score_elements.match) against the std == 1 rows (as nonauto_sim.do_main); the same in both output files;
  call class: nonauto_sim.call_class on the renamed DANTE domains (coding = coding unit in the call's middle 80%;
  proteinless = no Class_I|LTR domain fully inside the call span; partial = the rest);
  precision per call class = matched calls of the class / calls of the class, on two bases:
    std  (sim_nonauto_prf.tsv): matched one-to-one to a std == 1 truth row (as sim_callclass.tsv; calls on rebuilt
         hosts count as unmatched);
    ext  (sim_nonauto_prf_ext.tsv): matched one-to-one to any ext truth element (every truth.tsv row, rebuilt hosts
         included; score_elements ext scope, the basis of sim_class_prf.py).
F1 = 2PR/(P+R) from recall (element class) and precision (call class). Seeds pooled = summed counts.
Coordinates: DANTE GFF3 1-based closed -> 0-based half-open in evidence.dante_domains (start - 1, end); calls and
truth 0-based half-open.
Checks (fatal): recall n/tp equal experiments/nonauto/sim_class.tsv for every genome x tool x class; summed over call
classes, std basis: n_calls and matched equal the old sim_callclass.tsv "proteinless" row (= all calls); ext basis:
n_calls and matched equal the ext row (ro 0.9) of eval/sim/<S>/<tool>.metrics.tsv and matched / calls printed to 4
digits equals its precision; every DANTE sequence ID is in names.tsv.
Outputs (--out): sim_nonauto_prf.tsv (std basis), sim_nonauto_prf_ext.tsv (ext basis); diagnostics on stderr
(DANTE_LTR calls classed coding; matched calls whose call class equals the element class).
"""
import argparse
import csv
import math
import os
import sys
from collections import defaultdict
from multiprocessing import Pool

sys.dont_write_bytecode = True
B = os.environ["BENCH_ROOT"]
NA = f"{B}/experiments/nonauto"
sys.path.insert(0, NA)
import nonauto_sim as N  # noqa: E402  (imports queue before scripts/ goes on the path; see its gotcha)
from evidence import dante_domains  # noqa: E402
from score_elements import match, read_calls, read_truth  # noqa: E402

TOOLS = ["ltrquest", "ltr_retriever", "ltrharvest", "ltr_finder", "look4ltrs", "hirepeat", "hite", "ltrpred",
         "dante_ltr"]
DIAG = ["ltrquest_hc"]
SETS = [f"S{s}_{r}" for s in (1, 2) for r in (1, 2, 3)]
CLS = ["coding", "partial", "proteinless"]
BASES = ["std", "ext"]


def die(msg):
    sys.exit(f"error: {msg}")


def rename_doms(doms, names_path):
    """evidence.dante_domains dict keyed by prep_genome IDs -> keyed by original IDs (names.tsv new_id -> original_id)."""
    nm = {r["new_id"]: r["original_id"] for r in csv.DictReader(open(names_path), delimiter="\t")}
    miss = sorted(set(doms) - set(nm))
    if miss:
        die(f"{names_path}: DANTE sequence IDs missing from the names map: {miss[:5]}")
    return {nm[k]: v for k, v in doms.items()}


def run_ds(args):
    ds, lib = args
    T = N.truth(ds, lib)
    X = read_truth(f"{B}/sims/truth/{ds}.truth.tsv")          # ext truth (all rows)
    doms = rename_doms(dante_domains(f"{B}/runs/{ds}/dante_ltr/dante.gff3"), f"{B}/data/{ds}/names.tsv")
    out = {}
    for tool in TOOLS + DIAG:
        calls = read_calls(f"{B}/calls/{ds}/{tool}.tsv")
        c2t, t2c = match(calls, T, 0.9)
        c2x, _ = match(calls, X, 0.9)
        el = defaultdict(lambda: [0, 0])
        for ti, t in enumerate(T):
            if t["has_ltr"]:
                el[t["cls"]][0] += 1
                el[t["cls"]][1] += ti in t2c
        cc = defaultdict(lambda: [0, 0, 0])      # calls, matched (std), matched (ext)
        agree = [0, 0]
        for ci, c in enumerate(calls):
            k = N.call_class(doms, c)
            cc[k][0] += 1
            cc[k][2] += ci in c2x
            if ci in c2t:
                cc[k][1] += 1
                agree[0] += 1
                agree[1] += T[c2t[ci]]["cls"] == k
        out[tool] = dict(el={c: tuple(el[c]) for c in CLS}, cc={c: tuple(cc[c]) for c in CLS}, agree=tuple(agree))
    return ds, out


def prf(n_t, tp_t, n_c, tp_c):
    r = tp_t / n_t if n_t else float("nan")
    p = tp_c / n_c if n_c else float("nan")
    f = 2 * p * r / (p + r) if math.isfinite(p) and math.isfinite(r) and p + r > 0 else float("nan")
    return r, p, f


def fmt(x):
    return "NA" if isinstance(x, float) and not math.isfinite(x) else (f"{x:.4f}" if isinstance(x, float) else str(x))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", default=os.path.dirname(os.path.abspath(__file__)), help="output dir (default: script dir)")
    ap.add_argument("-t", "--threads", type=int, default=4, help="worker processes (default 4)")
    a = ap.parse_args(argv)
    lib = N.load_lib(NA)
    with Pool(a.threads) as pool:
        res = dict(pool.map(run_ds, [(d, lib) for d in SETS]))
    # checks against the existing nonauto tables and the ext metrics
    bad = []
    old_el = {(r["ds"], r["tool"], r["cls"]): (int(r["n"]), int(r["tp"]))
              for r in csv.DictReader(open(f"{NA}/sim_class.tsv"), delimiter="\t")}
    old_cc = {(r["ds"], r["tool"], r["call_cls"]): (int(r["n_calls"]), int(r["matched"]))
              for r in csv.DictReader(open(f"{NA}/sim_callclass.tsv"), delimiter="\t")}
    for ds in SETS:
        for tool in TOOLS + DIAG:
            o = res[ds][tool]
            for c in CLS:
                if o["el"][c] != old_el[(ds, tool, c)]:
                    bad.append(f"{ds} {tool} {c}: recall n/tp {o['el'][c]} vs sim_class.tsv {old_el[(ds, tool, c)]}")
            tot = tuple(sum(o["cc"][c][j] for c in CLS) for j in (0, 1, 2))
            if tot[:2] != old_cc[(ds, tool, "proteinless")] or any(old_cc[(ds, tool, c)] != (0, 0) for c in CLS[:2]):
                bad.append(f"{ds} {tool}: calls/matched {tot[:2]} vs sim_callclass.tsv proteinless row "
                           f"{old_cc[(ds, tool, 'proteinless')]}")
            ext = [r for r in csv.DictReader(open(f"{B}/eval/sim/{ds}/{tool}.metrics.tsv"), delimiter="\t")
                   if r["scope"] == "ext" and r["ro"] == "0.9"]
            if len(ext) != 1:
                die(f"{ds} {tool}: {len(ext)} ext ro 0.9 rows in metrics.tsv")
            e = ext[0]
            if (tot[0], tot[2]) != (int(e["n_calls"]), int(e["tp"])) or f"{tot[2] / tot[0]:.4f}" != e["precision"]:
                bad.append(f"{ds} {tool}: ext calls/matched {tot[0]}/{tot[2]} vs metrics {e['n_calls']}/{e['tp']} "
                           f"(precision {e['precision']})")
    if bad:
        die("checks failed:\n  " + "\n  ".join(bad))
    print("checks passed: recall n/tp = sim_class.tsv (60 genome x tool x 3 classes); std call totals = old "
          "sim_callclass.tsv all-calls row; ext calls/matched and precision string = eval/sim ext row (60 runs)",
          file=sys.stderr)
    rows = {b: [] for b in BASES}
    for st in ("S1", "S2"):
        for tool in TOOLS + DIAG:
            for seed in ("1", "2", "3", "pooled"):
                dss = [f"{st}_{seed}"] if seed != "pooled" else [f"{st}_{r}" for r in (1, 2, 3)]
                for c in CLS:
                    n_t = sum(res[d][tool]["el"][c][0] for d in dss)
                    tp_t = sum(res[d][tool]["el"][c][1] for d in dss)
                    n_c = sum(res[d][tool]["cc"][c][0] for d in dss)
                    for bi, b in enumerate(BASES):
                        tp_c = sum(res[d][tool]["cc"][c][1 + bi] for d in dss)
                        r, p, f = prf(n_t, tp_t, n_c, tp_c)
                        rows[b].append(dict(set=st, seed=seed, tool=tool, cls=c, n_truth=n_t, tp_truth=tp_t, recall=r,
                                            n_calls=n_c, tp_calls=tp_c, precision=p, f1=f, precision_basis=b,
                                            role="diagnostic" if tool in DIAG else "primary"))
            ag = [sum(res[f"{st}_{r}"][tool]["agree"][j] for r in (1, 2, 3)) for j in (0, 1)]
            dl = sum(res[f"{st}_{r}"][tool]["cc"]["coding"][0] for r in (1, 2, 3))
            al = sum(res[f"{st}_{r}"][tool]["cc"][c][0] for r in (1, 2, 3) for c in CLS)
            print(f"{st} {tool}: calls coding {dl}/{al} = {dl / al:.3f}; matched (std) calls whose call class = element "
                  f"class {ag[1]}/{ag[0]} = {ag[1] / ag[0]:.3f}", file=sys.stderr)
    cols = ["set", "seed", "tool", "cls", "n_truth", "tp_truth", "recall", "n_calls", "tp_calls", "precision", "f1",
            "precision_basis", "role"]
    for b, name in (("std", "sim_nonauto_prf.tsv"), ("ext", "sim_nonauto_prf_ext.tsv")):
        with open(f"{a.out}/{name}", "w") as f:
            f.write("\t".join(cols) + "\n")
            for r in rows[b]:
                f.write("\t".join(fmt(r[c]) for c in cols) + "\n")


if __name__ == "__main__":
    main()
