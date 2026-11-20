#!/usr/bin/env python3
"""Non-autonomous benchmark, real genomes, from existing approach-4 tables (eval/adj/<ds>/<tool>.calls.tsv).

Call class from REXdb Class_I|LTR domain records (runs/<ds>/dante_ltr/dante.gff3, genome-wide DANTE scan):
  coding      = approach-4 coding unit (column 'coding': >= 2 distinct domains of one superfamily, one strand, inside
                the internal region)
  proteinless = no Class_I|LTR domain record overlapping the call span (any score)
  partial     = the rest
Per tool x class: calls, library-supported share (approach-4 rule: >= 80% of span curated LTR-RT) and the same share
for the span moved 20-200 kb (null_supported, chance level), both hallmarks (TG..CA + exact 4-6 bp TSD), median K2P.
Reference loci: library-supported calls of the primary tools, union-find at reciprocal overlap >= 0.9 (approach 4a);
locus class = coding if any member is coding, proteinless if every member is, else partial; precise = some member
has both hallmarks; locus K2P = median member K2P. Per tool: share of loci matched (RO >= 0.9), per class, age bin.
95% CIs: 1,000 bootstrap resamples of 5 Mb genome blocks (numpy default_rng(1); approach 4's block rule).
Locus length = median member span; locus start/end = median member start/end (0-based, half-open, renamed
chromosome IDs of data/<ds>/genome.fa); bins < 2 kb, 2-4 kb, >= 4 kb (real_loci_len.tsv).
Writes real_calls.tsv, real_loci.tsv, real_loci_age.tsv, real_loci_len.tsv, real_locus_table.tsv to --out.
"""
import argparse
import bisect
import csv
import math
import os
import statistics
import sys
from collections import defaultdict
from multiprocessing import Pool

import numpy as np

B = os.environ["BENCH_ROOT"]
sys.path.append(f"{B}/scripts")
from adj_summary import loci, found_loci, load  # noqa: E402
from evidence import dante_domains  # noqa: E402

TOOLS = ["ltr_peel", "ltrharvest", "ltr_finder", "ltr_retriever", "look4ltrs", "hirepeat", "hite", "ltrpred",
         "dante_ltr"]
DS = ["arab", "rice", "finch", "maize", "gorilla"]
CLS = ["coding", "partial", "proteinless"]
BLOCK, NB = 5_000_000, 1000
LENB = [(0, 2000, "<2kb"), (2000, 4000, "2-4kb"), (4000, 10 ** 9, ">=4kb")]
AGEB = [(0, 0.01, "<0.01"), (0.01, 0.03, "0.01-0.03"), (0.03, 0.1, "0.03-0.1"), (0.1, 99, ">=0.1")]


def n_dom(doms, c, s, e):
    d = doms.get(c)
    if not d:
        return 0
    starts, recs = d
    i = bisect.bisect_left(starts, s - 20000)
    n = 0
    while i < len(recs) and recs[i][0] < e:
        n += recs[i][1] > s
        i += 1
    return n


def k2p(r):
    try:
        x = float(r["k2p"])
        return x if math.isfinite(x) else None
    except ValueError:
        return None


def agebin(x):
    if x is None:
        return "NA"
    for lo, hi, lab in AGEB:
        if lo <= x < hi:
            return lab
    return "NA"


def run(ds):
    doms = dante_domains(f"{B}/runs/{ds}/dante_ltr/dante.gff3")
    rows = {}
    for t in TOOLS:
        p = f"{B}/eval/adj/{ds}/{t}.calls.tsv"
        if os.path.exists(p):
            rows[t] = load(p)
            for r in rows[t]:
                nd = n_dom(doms, r["chrom"], r["start"], r["end"])
                r["cls"] = "coding" if r["coding"] == "1" else ("proteinless" if nd == 0 else "partial")
                r["n_dom"] = nd
                r["hall"] = r["tgca"] == "1" and r["tsd"] == "1"
                r["k"] = k2p(r)
    blocks = sorted({(r["chrom"], r["start"] // BLOCK) for rr in rows.values() for r in rr})
    bpos = {b: i for i, b in enumerate(blocks)}
    rng = np.random.default_rng(1)
    W = np.stack([np.bincount(d, minlength=len(blocks)) for d in rng.integers(0, len(blocks), (NB, len(blocks)))])

    def boot_share(sel, key):
        if not sel:
            return "NA", "NA"
        bi = np.array([bpos[(r["chrom"], r["start"] // BLOCK)] for r in sel])
        x = np.array([r[key] == "1" for r in sel], float)
        num, den = W @ np.bincount(bi, weights=x, minlength=len(blocks)), W @ np.bincount(bi, minlength=len(blocks))
        v = num[den > 0] / den[den > 0]
        return f"{np.percentile(v, 2.5):.4f}", f"{np.percentile(v, 97.5):.4f}"
    out_calls = []
    for t, rr in rows.items():
        for c in CLS + ["all"]:
            sel = [r for r in rr if c == "all" or r["cls"] == c]
            n = len(sel)
            ks = [r["k"] for r in sel if r["k"] is not None]
            slo, shi = boot_share(sel, "supported")
            nlo, nhi = boot_share(sel, "null_supported")
            out_calls.append(dict(sup_lo=slo, sup_hi=shi, null_lo=nlo, null_hi=nhi,
                ds=ds, tool=t, cls=c, n_calls=n, share=f"{n / len(rr):.4f}" if rr else "NA",
                supported=sum(r["supported"] == "1" for r in sel),
                supported_share=f"{sum(r['supported'] == '1' for r in sel) / n:.4f}" if n else "NA",
                null_supported_share=f"{sum(r['null_supported'] == '1' for r in sel) / n:.4f}" if n else "NA",
                junk_share=f"{sum(r['junk'] == '1' for r in sel) / n:.4f}" if n else "NA",
                null_junk_share=f"{sum(r['null_junk'] == '1' for r in sel) / n:.4f}" if n else "NA",
                unsup_junk_share=(f"{sum(r['junk'] == '1' for r in sel if r['supported'] != '1') / max(1, sum(r['supported'] != '1' for r in sel)):.4f}"),
                hallmarks_share=f"{sum(r['hall'] for r in sel) / n:.4f}" if n else "NA",
                sup_hallmarks_share=(f"{sum(r['hall'] for r in sel if r['supported'] == '1') / max(1, sum(r['supported'] == '1' for r in sel)):.4f}"),
                median_k2p=f"{statistics.median(ks):.4f}" if ks else "NA",
                median_len=f"{statistics.median([r['end'] - r['start'] for r in sel]):.0f}" if n else "NA"))
    # reference loci
    mem = []
    for t, rr in rows.items():
        for k, r in enumerate(rr):
            if r["supported"] == "1":
                mem.append((r["chrom"], r["start"], r["end"], (t, k)))
    lab, nl = loci(mem)
    info = defaultdict(lambda: dict(cls=set(), hall=False, ks=[], lens=[], tools=set(), chrom=None, ss=[], es=[]))
    for c, s, e, (t, k) in mem:
        r = rows[t][k]
        d = info[lab[(t, k)]]
        d["chrom"] = c
        d["ss"].append(s)
        d["es"].append(e)
        d["cls"].add(r["cls"])
        d["hall"] |= r["hall"]
        d["tools"].add(t)
        d["lens"].append(e - s)
        if r["k"] is not None:
            d["ks"].append(r["k"])
    lcls, lprec, lage, llen = {}, {}, {}, {}
    for L, d in info.items():
        ml = statistics.median(d["lens"])
        llen[L] = next(lab_ for lo, hi, lab_ in LENB if lo <= ml < hi)
        lcls[L] = "coding" if "coding" in d["cls"] else ("proteinless" if d["cls"] == {"proteinless"} else "partial")
        lprec[L] = d["hall"]
        lage[L] = agebin(statistics.median(d["ks"]) if d["ks"] else None)
    member_spans = [(c, s, e, lab[key]) for c, s, e, key in mem]
    out_loci, out_age, out_len = [], [], []
    ntot = defaultdict(int)
    for L in info:
        ntot[(lcls[L], "len", llen[L])] += 1
        ntot[(lcls[L], "all")] += 1
        if lprec[L]:
            ntot[(lcls[L], "precise")] += 1
        ntot[(lcls[L], "all", lage[L])] += 1
        if lprec[L]:
            ntot[(lcls[L], "precise", lage[L])] += 1
    lblock = {}
    for c, s_, e, key in mem:
        lblock.setdefault(lab[key], bpos[(c, s_ // BLOCK)])
    for t, rr in rows.items():
        hit = found_loci(member_spans, rr, 0.9)
        rci = {}
        for c in CLS:
            for sc in ("all", "precise"):
                Ls = [L for L in info if lcls[L] == c and (sc == "all" or lprec[L])]
                if not Ls:
                    rci[(c, sc)] = ("NA", "NA")
                    continue
                bi = np.array([lblock[L] for L in Ls])
                x = np.array([L in hit for L in Ls], float)
                num, den = W @ np.bincount(bi, weights=x, minlength=len(blocks)), W @ np.bincount(bi, minlength=len(blocks))
                v = num[den > 0] / den[den > 0]
                rci[(c, sc)] = (f"{np.percentile(v, 2.5):.4f}", f"{np.percentile(v, 97.5):.4f}")
        h = defaultdict(int)
        for L in hit:
            h[(lcls[L], "len", llen[L])] += 1
            h[(lcls[L], "all")] += 1
            h[(lcls[L], "all", lage[L])] += 1
            if lprec[L]:
                h[(lcls[L], "precise")] += 1
                h[(lcls[L], "precise", lage[L])] += 1
        for c in CLS:
            for _, _, lb in LENB:
                n = ntot[(c, "len", lb)]
                out_len.append(dict(ds=ds, tool=t, cls=c, len_bin=lb, n_loci=n, found=h[(c, "len", lb)],
                                    recall=f"{h[(c, 'len', lb)] / n:.4f}" if n else "NA"))
            for sc in ("all", "precise"):
                n = ntot[(c, sc)]
                out_loci.append(dict(ds=ds, tool=t, pool_tools=len(rows), cls=c, scope=sc, n_loci=n,
                                     found=h[(c, sc)], recall=f"{h[(c, sc)] / n:.4f}" if n else "NA",
                                     recall_lo=rci[(c, sc)][0], recall_hi=rci[(c, sc)][1]))
                for _, _, ab in AGEB + [(0, 0, "NA")]:
                    n = ntot[(c, sc, ab)]
                    out_age.append(dict(ds=ds, tool=t, cls=c, scope=sc, age=ab, n_loci=n, found=h[(c, sc, ab)],
                                        recall=f"{h[(c, sc, ab)] / n:.4f}" if n else "NA"))
    # locus-level table for drawing
    lrows = [dict(ds=ds, locus=L, chrom=d["chrom"], start=int(statistics.median(d["ss"])),
                  end=int(statistics.median(d["es"])), cls=lcls[L], precise=int(lprec[L]), age=lage[L],
                  k2p=f"{statistics.median(d['ks']):.4f}" if d["ks"] else "NA",
                  len=int(statistics.median(d["lens"])), n_tools=len(d["tools"]), tools=",".join(sorted(d["tools"])))
             for L, d in info.items()]
    print(f"nonauto_real: {ds} done ({nl} loci)", file=sys.stderr, flush=True)
    return out_calls, out_loci, out_age, lrows, out_len


def write(path, rows):
    with open(path, "w") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", default=os.path.dirname(os.path.abspath(__file__)), help="output dir (default: script dir)")
    ap.add_argument("--ds", nargs="+", default=DS, help="datasets (default: %(default)s)")
    a = ap.parse_args()
    with Pool(len(a.ds)) as pool:
        res = pool.map(run, a.ds)
    write(f"{a.out}/real_calls.tsv", [r for x in res for r in x[0]])
    write(f"{a.out}/real_loci.tsv", [r for x in res for r in x[1]])
    write(f"{a.out}/real_loci_age.tsv", [r for x in res for r in x[2]])
    write(f"{a.out}/real_locus_table.tsv", [r for x in res for r in x[3]])
    write(f"{a.out}/real_loci_len.tsv", [r for x in res for r in x[4]])


if __name__ == "__main__":
    main()
