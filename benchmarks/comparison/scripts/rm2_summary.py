#!/usr/bin/env python3
"""BM_RM2 summary: curated families found perfect / good / present by a test library, LTR families and all families.

Input: the working directory of scripts/vendor/get_family_summary_paper.sh (RepeatModeler2 benchmark script by J.
Flynn, as shipped with HiTE f5f2c7f bin/), run on the RepeatMasker .out of the TEST library masked with the CURATED
library (-nolow). Its definitions: perfect = one test sequence and one family cover each other > 95% at < 5%
divergence; good = test sequences together cover > 95% of the family at < 5% divergence (not perfect); present =
together cover > 80% at < 20% divergence. LTR families = curated entries whose '#class' matches the lib-test LTR regex.
usage: rm2_summary.py --workdir DIR --curated lib.fa --out OUT.tsv
"""
import argparse
import os
import re

LTR = re.compile("RLG|RLC|RLB|RLR|RLE|LTR|RLX|Gypsy|Copia", re.I)


def families(path):
    return set(open(path).read().split()) if os.path.exists(path) else set()


def main(argv=None):
    ap = argparse.ArgumentParser(description="Count BM_RM2 perfect/good/present curated families (LTR and all).")
    ap.add_argument("--workdir", required=True, help="directory where get_family_summary_paper.sh ran")
    ap.add_argument("--curated", required=True, help="curated library FASTA (headers name#class)")
    ap.add_argument("--out", required=True, help="output TSV (one row)")
    a = ap.parse_args(argv)
    cls = {}
    for line in open(a.curated):
        if line.startswith(">"):
            name, _, c = line[1:].split()[0].partition("#")
            cls[name] = c
    ltr = {n for n, c in cls.items() if LTR.search(c)}
    perfect = families(os.path.join(a.workdir, "perfect.families"))
    good = families(os.path.join(a.workdir, "good.families"))
    present = families(os.path.join(a.workdir, "present.all.families"))
    with open(a.out, "w") as f:
        f.write("set\tn_curated\tperfect\tgood\tperfect_or_good\tpresent\n")
        for tag, fam in (("LTR", ltr), ("all", set(cls))):
            p, g, pr = len(perfect & fam), len(good & fam), len(present & fam)
            f.write(f"{tag}\t{len(fam)}\t{p}\t{g}\t{p + g}\t{pr}\n")
    print(f"rm2_summary: {a.out} LTR perfect {len(perfect & ltr)} good {len(good & ltr)} present {len(present & ltr)}"
          f" of {len(ltr)}")


if __name__ == "__main__":
    main()
