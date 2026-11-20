#!/usr/bin/env python3
"""Non-autonomous benchmark, simulations (approach 1 S1/S2 and approach 1b divergence series), from existing calls.

Element class = the class of its PrinTE library entry (lib/lib_class.tsv, classify_lib.py): coding (>= 2 REXdb
LTR-RT domains of one superfamily on one strand, internal region), proteinless (no REXdb LTR-RT domain record),
partial (rest). Only elements whose entry has an LTR pair (has_ltr = 1) are scored, as in div_score.py.
Matching: score_elements.match, one-to-one, reciprocal overlap >= 0.9, against the PrinTE standard set (std == 1).
Calls are also classed by their own REXdb domains in the simulated genome (runs/<ds>/dante_ltr/dante.gff3, the same
rule on the call's middle 80%, i.e. approach 4's fallback internal region), for precision per call class.
Writes (to --out): sim_class.tsv, sim_len.tsv, sim_callclass.tsv, div_class.tsv, div_class_d50.tsv, div_class_bins.tsv
"""
import argparse
import csv
import glob
import math
import os
import statistics
import sys
from collections import defaultdict
import queue  # noqa: F401  (stdlib, cached before div_score puts scripts/ first: scripts/queue.py shadows it)
import multiprocessing.pool  # noqa: F401
from multiprocessing import Pool

B = os.environ["BENCH_ROOT"]
sys.path.append(f"{B}/scripts")   # appended: scripts/queue.py would shadow the stdlib queue
from score_elements import match, read_calls  # noqa: E402
from div_score import element_stats, status, d50, snap, wilson, tsd_adjusted  # noqa: E402
from evidence import dante_domains, coding_unit  # noqa: E402

TOOLS = ["ltr_peel", "ltrharvest", "ltr_finder", "ltr_retriever", "look4ltrs", "hirepeat", "hite", "ltrpred",
         "dante_ltr", "ltr_peel_hc"]
CLS = ["coding", "partial", "proteinless"]
LBINS = [0, 2000, 3000, 4000, 5000, 6000, 8000, 10 ** 9]
KBINW = 0.02


def lbin(L):
    for lo, hi in zip(LBINS, LBINS[1:]):
        if lo <= L < hi:
            return f"{lo}-{hi}"


def load_lib(od):
    return {r["family"]: r for r in csv.DictReader(open(f"{od}/lib/lib_class.tsv"), delimiter="\t")}


def truth(ds, lib):
    out = []
    for r in csv.DictReader(open(f"{B}/sims/truth/{ds}.truth.tsv"), delimiter="\t"):
        if r["std"] != "1":
            continue
        e = lib[r["family"]]
        out.append(dict(chrom=r["chrom"], start=int(r["start"]), end=int(r["end"]), cls=e["cls"],
                        has_ltr=e["has_ltr"] == "1", ltr_len=int(r["ltr_len"]), tsd=r["tsd"], fam=r["family"]))
    return out


def call_class(doms, c):
    s, e = c["start"], c["end"]
    m = (e - s) // 10
    d = doms.get(c["chrom"])
    if not d:
        return "proteinless"
    starts, recs = d
    import bisect
    i = bisect.bisect_left(starts, s - 20000)
    n = 0
    while i < len(recs) and recs[i][0] < e:
        if recs[i][0] >= s and recs[i][1] <= e:
            n += 1
        i += 1
    if coding_unit(doms, c["chrom"], s + m, e - m):
        return "coding"
    return "proteinless" if n == 0 else "partial"


def do_main(args):
    ds, lib = args
    T = truth(ds, lib)
    doms = dante_domains(f"{B}/runs/{ds}/dante_ltr/dante.gff3")
    rows, lrows, crows, erows = [], [], [], []
    for tool in TOOLS:
        p = f"{B}/calls/{ds}/{tool}.tsv"
        if not os.path.exists(p):
            continue
        calls = read_calls(p)
        c2t, t2c = match(calls, T, 0.9)
        n = defaultdict(lambda: [0, 0])
        nl = defaultdict(lambda: [0, 0])
        for ti, t in enumerate(T):
            if not t["has_ltr"]:
                continue
            erows.append((ds, tool, ti, t["fam"], t["cls"], t["end"] - t["start"], int(ti in t2c)))
            n[t["cls"]][0] += 1
            n[t["cls"]][1] += ti in t2c
            k = (t["cls"], lbin(t["end"] - t["start"]))
            nl[k][0] += 1
            nl[k][1] += ti in t2c
        for c in CLS:
            rows.append(dict(ds=ds, set=ds[:2], seed=ds[3], tool=tool, cls=c, n=n[c][0], tp=n[c][1],
                             recall=f"{n[c][1] / n[c][0]:.4f}" if n[c][0] else "NA"))
        for (c, b), (a, h) in nl.items():
            lrows.append(dict(ds=ds, set=ds[:2], seed=ds[3], tool=tool, cls=c, len_bin=b, n=a, tp=h))
        cc = defaultdict(lambda: [0, 0, 0])   # calls, matched, matched to a proteinless element
        for ci, c in enumerate(calls):
            k = call_class(doms, c)
            cc[k][0] += 1
            if ci in c2t:
                cc[k][1] += 1
                cc[k][2] += T[c2t[ci]]["cls"] == "proteinless"
        for k in CLS:
            crows.append(dict(ds=ds, set=ds[:2], seed=ds[3], tool=tool, call_cls=k, n_calls=cc[k][0],
                              matched=cc[k][1], matched_proteinless=cc[k][2],
                              precision=f"{cc[k][1] / cc[k][0]:.4f}" if cc[k][0] else "NA"))
    return rows, lrows, crows, erows


def do_div(args):
    ds, lib = args
    rep, k = int(ds[2]), int(ds.split("_")[1])
    T = truth(ds, lib)
    tl = [len(t["tsd"]) if t["tsd"] not in ("", "NA") else 0 for t in T]
    es = element_stats(f"{B}/sims/div/D{rep}/{snap(k)}.fasta",
                       [(t["chrom"], t["start"], t["end"], t["ltr_len"], l) for t, l in zip(T, tl)])
    nan3 = (float("nan"), None, float("nan"))
    ev = [es.get((t["chrom"], t["start"], t["end"]), nan3) for t in T]
    tk = [x[0] for x in ev]
    qc = {}
    for c in CLS:   # share of TSD pairs expected identical at the element's own LTR divergence (div_score.tsd_adjusted)
        kn = [(1 - x[2]) ** l for x, l, t in zip(ev, tl, T) if t["cls"] == c and math.isfinite(x[0]) and x[1] is not None]
        qc[c] = sum(kn) / len(kn) if kn else float("nan")
    fin = [x for x in tk if math.isfinite(x)]
    med = statistics.median(fin)
    medc = {c: statistics.median([x for x, t in zip(tk, T) if math.isfinite(x) and t["cls"] == c] or [float("nan")])
            for c in CLS}
    rows, bins, erows = [], [], []
    for tool in TOOLS:
        st = status(ds, tool) if os.getcwd() == B else None
        if st not in ("ok", "no_candidates"):
            continue
        t2c = {}
        if st == "ok":
            _, t2c = match(read_calls(f"{B}/calls/{ds}/{tool}.tsv"), T, 0.9)
        for c in CLS:
            idx = [i for i, t in enumerate(T) if t["cls"] == c and math.isfinite(tk[i])]
            tp = sum(1 for i in idx if i in t2c)
            strata = [0, 0, 0, 0]
            for i in idx:
                if ev[i][1] is not None:
                    j = 0 if ev[i][1] else 2
                    strata[j] += 1
                    strata[j + 1] += i in t2c
            adj = tsd_adjusted(strata, qc[c])
            rows.append(dict(ds=ds, seed=rep, timepoint=k, median_k2p=f"{med:.4f}", cls_median_k2p=f"{medc[c]:.4f}",
                             tool=tool, status=st, cls=c, n=len(idx), tp=tp,
                             recall=f"{tp / len(idx):.4f}" if idx else "NA",
                             n_tsd_same=strata[0], tp_tsd_same=strata[1], n_tsd_diff=strata[2], tp_tsd_diff=strata[3],
                             tsd_same_expected=f"{qc[c]:.4f}" if math.isfinite(qc[c]) else "NA",
                             recall_tsd_adj=f"{adj:.4f}" if math.isfinite(adj) else "NA"))
            for i in idx:
                bins.append((tool, c, int(tk[i] / KBINW), i in t2c))
                erows.append((ds, rep, k, tool, i, T[i]["fam"], c, T[i]["end"] - T[i]["start"], f"{tk[i]:.4f}",
                               "NA" if ev[i][1] is None else int(ev[i][1]), int(i in t2c)))
    return rows, bins, erows


def write(path, rows):
    with open(path, "w") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", default=os.path.dirname(os.path.abspath(__file__)), help="output dir (default: script dir)")
    ap.add_argument("-t", "--threads", type=int, default=8, help="worker processes (default 8)")
    a = ap.parse_args()
    os.chdir(B)
    lib = load_lib(a.out)
    main_ds = [f"S{s}_{r}" for s in (1, 2) for r in (1, 2, 3)]
    div_ds = sorted(os.path.basename(p)[:-len(".truth.tsv")] for p in glob.glob(f"{B}/sims/truth/SD*_*.truth.tsv"))
    with Pool(a.threads) as pool:
        res = pool.map(do_main, [(d, lib) for d in main_ds])
        dres = pool.map(do_div, [(d, lib) for d in div_ds])
    write(f"{a.out}/sim_class.tsv", [r for x in res for r in x[0]])
    write(f"{a.out}/sim_len.tsv", [r for x in res for r in x[1]])
    write(f"{a.out}/sim_callclass.tsv", [r for x in res for r in x[2]])
    with open(f"{a.out}/sim_elements.tsv", "w") as f:
        f.write("ds\ttool\telement\tfamily\tcls\tlen\thit\n")
        for x in res:
            for e in x[3]:
                f.write("\t".join(map(str, e)) + "\n")
    with open(f"{a.out}/div_elements.tsv", "w") as f:
        f.write("ds\tseed\ttimepoint\ttool\telement\tfamily\tcls\tlen\tk2p\ttsd_same\thit\n")
        for x in dres:
            for e in x[2]:
                f.write("\t".join(map(str, e)) + "\n")
    drows = [r for x in dres for r in x[0]]
    write(f"{a.out}/div_class.tsv", drows)
    agg = defaultdict(lambda: [0, 0])
    for x in dres:
        for tool, c, b, hit in x[1]:
            agg[(tool, c, b)][0] += 1
            agg[(tool, c, b)][1] += hit
    brows = []
    for (tool, c, b), (n, h) in sorted(agg.items()):
        lo, hi = wilson(h, n)
        brows.append(dict(tool=tool, cls=c, k2p_lo=f"{b * KBINW:.2f}", k2p_hi=f"{(b + 1) * KBINW:.2f}", n=n, tp=h,
                          recall=f"{h / n:.4f}", ci_lo=f"{lo:.4f}", ci_hi=f"{hi:.4f}"))
    write(f"{a.out}/div_class_bins.tsv", brows)
    series = defaultdict(list)
    for r in drows:
        if r["recall"] != "NA":
            series[(r["tool"], r["seed"], r["cls"])].append((r["timepoint"], float(r["cls_median_k2p"]),
                                                             float(r["recall"])))
    sadj = defaultdict(list)
    for r in drows:
        if r["recall_tsd_adj"] != "NA":
            sadj[(r["tool"], r["seed"], r["cls"])].append((r["timepoint"], float(r["cls_median_k2p"]),
                                                           float(r["recall_tsd_adj"])))
    d50rows = []
    for (tool, seed, c), pts in sorted(series.items()):
        pts.sort()
        pa = sorted(sadj.get((tool, seed, c), []))
        d50rows.append(dict(tool=tool, seed=seed, cls=c, recall_t0=f"{pts[0][2]:.4f}" if pts[0][0] == 0 else "NA",
                            d50_k2p=d50(pts), d50_k2p_tsd_adj=d50(pa) if pa else "NA"))
    write(f"{a.out}/div_class_d50.tsv", d50rows)


if __name__ == "__main__":
    main()
