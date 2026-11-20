#!/usr/bin/env python3
"""How much could curated-library gaps inside real elements explain masking false positives? (post hoc diagnostic)

Real elements = calls with an LTR-RT protein coding unit (approach 4 'coding'), pooled over the 9 primary callers;
'precise' = those calls that also have an exact TSD and TG..CA ends. Gap = bases inside those spans that the curated
reference annotation does not mark as LTR (proxies.py LTR pattern on mask/<ds>/curated/<ds>.fa.out). If a caller
masked every gap base, the masking test would count them as false positives although they lie in real elements, so
the gap (Mb) bounds the masking-FP penalty that library incompleteness inside protein-confirmed elements can cause.
Compared with each caller's masking FP (report/tables/masking.tsv, intact-call libraries).
usage: annotation_gap.py [-o experiments/a4_diagnostics]
"""
import argparse
import csv
import os
import sys

import numpy as np

B = os.environ["BENCH_ROOT"]
sys.path.insert(0, f"{B}/scripts")
import adj_summary as su  # noqa: E402
import proxies  # noqa: E402

TOOLS = ["ltr_peel", "ltr_retriever", "ltrharvest", "ltr_finder", "look4ltrs", "hirepeat", "hite", "ltrpred",
         "dante_ltr"]


def main(argv=None):
    ap = argparse.ArgumentParser(description="Curated-annotation gap inside protein-confirmed LTR-RT calls.")
    ap.add_argument("-o", "--outdir", default=f"{B}/experiments/a4_diagnostics", help="output directory")
    a = ap.parse_args(argv)
    fp = {(r["ds"], r["lib"]): r for r in csv.DictReader(open(f"{B}/report/tables/masking.tsv"), delimiter="\t")}
    rows = [("ds", "set", "n_calls", "span_mb", "curated_ltr_mb_in_span", "gap_mb", "genome_curated_ltr_mb")]
    cmp_rows = [("ds", "tool", "masking_FP_mb", "gap_all_mb", "gap_precise_mb", "gap_all_over_FP")]
    for ds in ["arab", "rice", "finch", "maize", "gorilla"]:
        rb = su.load_ds(f"{B}/eval/adj/{ds}", TOOLS)
        if len(rb) < 5:
            continue
        lens = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(f"{B}/data/{ds}/genome.fa.fai")}
        ltr = {c: np.zeros(n, bool) for c, n in lens.items()}
        for c, iv in proxies.ref_ltr(f"{B}/mask/{ds}/curated/{ds}.fa.out").items():
            for s, e in iv:
                ltr[c][s:e] = True
        genome_ltr = sum(int(v.sum()) for v in ltr.values())
        gaps = {}
        for tag, keep in (("all", lambda r: r["coding"] == "1"),
                          ("precise", lambda r: r["coding"] == "1" and r["tsd"] == "1" and r["tgca"] == "1")):
            span = {c: np.zeros(n, bool) for c, n in lens.items()}
            n = 0
            for rs in rb.values():
                for r in rs:
                    if keep(r):
                        span[r["chrom"]][r["start"]:r["end"]] = True
                        n += 1
            sp = sum(int(v.sum()) for v in span.values())
            inl = sum(int((span[c] & ltr[c]).sum()) for c in lens)
            gaps[tag] = (sp - inl) / 1e6
            rows.append((ds, tag, n, f"{sp / 1e6:.2f}", f"{inl / 1e6:.2f}", f"{gaps[tag]:.2f}", f"{genome_ltr / 1e6:.2f}"))
            del span
        for t in TOOLS:
            r = fp.get((ds, f"{t}.intact"))
            if r and r["FP_mb"] not in (".", ""):
                f = float(r["FP_mb"])
                cmp_rows.append((ds, t, f"{f:.2f}", f"{gaps['all']:.2f}", f"{gaps['precise']:.2f}",
                                 f"{gaps['all'] / f:.3f}" if f else "."))
        print(f"annotation_gap: {ds} gap all {gaps['all']:.2f} Mb, precise {gaps['precise']:.2f} Mb "
              f"(curated LTR {genome_ltr / 1e6:.1f} Mb)", file=sys.stderr)
    for fn, rr in (("annotation_gap.tsv", rows), ("annotation_gap_vs_fp.tsv", cmp_rows)):
        with open(os.path.join(a.outdir, fn), "w") as f:
            f.write("".join("\t".join(map(str, r)) + "\n" for r in rr))
    print(f"annotation_gap: {a.outdir}/annotation_gap.tsv, annotation_gap_vs_fp.tsv")


if __name__ == "__main__":
    main()
