#!/usr/bin/env python3
"""Fig. S4 figS4_ltrkin (publication set): genotyping real insertion polymorphisms in one Cardamine
amara plant from long reads, details. Panels A and B follow the genotyping benchmark's own yield figure (same example and rule), redrawn as fig_benchmark
A draws calls; C and D break fig_benchmark H down by element class, in its style.

Panels:
  A  one insertion polymorphism on haplotype 2, chosen at run time by the source figure's rule: haplotype 1 with the
     empty site on top, the aligned flanks as ribbons, the element with both LTRs and TSDs on haplotype 2 (gray, truth),
     the plant's HiFi reads that span the site (carrying or skipping the element), and every caller's calls in the
     window in its tool color with the features that caller reports (LTR end blocks, TSD boxes, strand chevron,
     protein-domain boxes) from its call table or native output; a feature not reported is not drawn.
  B  both junctions at base-pair resolution: the TSDs on haplotype 2, the single TSD at the haplotype 1 empty site, and
     each finding caller's called end and the TSD it reports at these positions.
  C  share of loci genotyped correctly by nesting class (non-nested / tenant / host), one 100% bar per caller and
     haplotype: genotyped correctly in the caller's color, split into loci other callers also nominated (solid) and
     loci nominated by this caller alone (nominators = this caller only; white hatch over its color, true width), the
     rest light gray, correct share with its 95% CI at the bar's end (experiments/figpub/ltrkin_breakdown.tsv).
  D  the same by protein content (coding / partly coding / protein-less).
Numbers are read at run time from ltrkin_geno/report/tables/{ltrkin_yield,recall}.tsv, ltrkin_geno/report/checks/
{subsets,ends}.tsv and experiments/figpub/{ltrkin_breakdown,ltrkin_locus_class}.tsv; the per-locus outcomes come from
ltrkin_geno/scripts/score.py (score.score_dataset) and are self-checked against those tables (and the class table's
per-locus outcomes) before anything is drawn. Missing results: tool-side failure = dark gray "tool failed".
Coordinates: pool, calls and pysam are 0-based half-open; truth TSVs and ltrkin genotypes are 1-based closed (converted
by the project readers); PAF 0-based half-open; native GFF3 / pass list / scn 1-based closed (start - 1 on reading).
Outputs (this directory only): figS4_ltrkin.pdf (figure + legend page), figS4_ltrkin.png (300 dpi, figure only),
figS4_ltrkin_numbers.tsv. Nothing is written in ltrkin_geno/ (sys.dont_write_bytecode; run with PYTHONDONTWRITEBYTECODE=1).
Source data: the default (full) mode writes every value drawn or printed to source_data/figS4_ltrkin/*.tsv and draws
from those tables read back; --from-source-data redraws from them alone (no BENCH_ROOT, no pysam; the checks against
the raw data run in full mode only).
Usage: figS4_ltrkin.py [-v] [--draft PNG] [--from-source-data]
"""
import argparse
import collections
import csv
import glob
import gzip
import math
import os
import re
import statistics
import sys
from decimal import ROUND_HALF_UP, Decimal
from pathlib import Path

sys.dont_write_bytecode = True
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np  # noqa: E402

BENCH = Path(os.environ["BENCH_ROOT"]) if os.environ.get("BENCH_ROOT") else None
LG = BENCH / "ltrkin_geno" if BENCH else None           # data root of the genotyping benchmark
LGM = Path(__file__).resolve().parent.parent / "ltrkin_geno"   # its modules
PRIMARY = None                                           # lg_config.TOOLS (full mode)


def import_raw_modules():
    """Modules needed only to read the full dataset: the ltrkin_geno scripts, pysam and parasail."""
    global score, lg_config, gff_id, final_calls_path, read_calls, read_truth, revcomp, parasail, pysam, PRIMARY
    sys.path.insert(0, str(LGM / "scripts"))
    import score
    import lg_config
    from calls_to_gff import gff_id
    from lg_io import final_calls_path, read_calls
    from pool import read_truth, revcomp
    import parasail
    import pysam
    PRIMARY = list(lg_config.TOOLS)                      # the secondary row (ltr_peel_hc) is never drawn


sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # poafig.py sits beside this script
import poafig  # noqa: E402
import matplotlib as mpl  # noqa: E402
import matplotlib.font_manager as fm  # noqa: E402

# the shared matplotlib font cache on this host lists no system fonts: register Liberation in-process (cache untouched)
for _f in sorted(glob.glob("/usr/share/fonts/truetype/liberation2/LiberationSans-*.ttf")
                 + glob.glob("/usr/share/fonts/truetype/liberation2/LiberationMono-*.ttf")):
    fm.fontManager.addfont(_f)
poafig.setup_mpl()
mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def _advance_pt_mpl311(text, style):
    """poafig._advance_pt for matplotlib >= 3.11 (layout() lost kern_mode; LayoutItem.glyph_idx -> glyph_index).
    Same measurement as backend_pdf: x of the last glyph + its unhinted advance, font set to the size in points."""
    from matplotlib import _text_helpers
    from matplotlib.ft2font import LoadFlags
    font, _ = poafig._legend_font(style)
    last = None
    for last in _text_helpers.layout(text, font):
        pass
    if last is None:
        return 0.0
    g = last.ft_object.load_glyph(last.glyph_index, flags=LoadFlags.NO_HINTING)
    return last.x + g.linearHoriAdvance / 65536


def _needs_patch():
    import inspect
    from matplotlib import _text_helpers
    return "kern_mode" not in inspect.signature(_text_helpers.layout).parameters


if _needs_patch():
    poafig._advance_pt = _advance_pt_mpl311          # in-process only; poafig.py is not modified
from matplotlib.collections import PolyCollection  # noqa: E402
from matplotlib.colors import to_rgb  # noqa: E402
from matplotlib.patches import Polygon, Rectangle  # noqa: E402

STEM = "figS4_ltrkin"
HERE = Path(__file__).resolve().parent
SRC_DIR = HERE / "source_data" / STEM                # the tables the figure is drawn from
FIGPUB = BENCH / "experiments" / "figpub" if BENCH else None
TAB = LG / "report" / "tables" if LG else None
CHK = LG / "report" / "checks" if LG else None
MM = poafig.MM
DS = ["camara_hap2", "camara_hap1"]                  # primary first
HAPN = {"camara_hap2": "haplotype 2", "camara_hap1": "haplotype 1"}
NAME = {"ltr_peel": "LTR_peel", "hirepeat": "HiRepEat", "ltr_retriever": "LTR_retriever", "hite": "HiTE",
        "ltrharvest": "LTRharvest", "ltr_finder": "LTR_Finder", "dante_ltr": "DANTE_LTR", "ltrpred": "LTRpred",
        "look4ltrs": "Look4LTRs"}
# shared tool order and palette of the publication set
ORDER = ["ltr_peel", "ltr_retriever", "ltrharvest", "ltr_finder", "look4ltrs", "hirepeat", "hite", "ltrpred",
         "dante_ltr"]
COL = {"ltr_peel": "#D55E00", "ltr_retriever": "#0072B2", "ltrharvest": "#56B4E9", "ltr_finder": "#009E73",
       "look4ltrs": "#CC79A7", "hirepeat": "#E69F00", "hite": "#000000", "ltrpred": "#A6761D", "dante_ltr": "#6A3D9A"}

# colors: the true element in grays (truth context; vermillion names LTR_peel in this set); calls in each caller's tool
# color as fig_benchmark A (body = tint, LTR end blocks = full color, TSD = white box, strand = chevron, protein domain =
# dark box); breakdown bars as fig_benchmark H (correct in the caller's color, rest light gray, tool failed
# dark gray)
EL_LTR, EL_INT, EL_READ = "#3C3C3C", "#8C8C8C", "#3C3C3C"
GRAY, LGRAY, INK, MUTED = "#8C8C8C", "#BDBDBD", poafig.INK, poafig.MUTED
READC = "#A6A6A6"
RIBBON = "#E3E3E3"
RESTC = "#DDDDDD"                                    # C, D: not genotyped correctly
DOMC = "#2E2E2E"                                     # A: protein domains reported by a caller
ROWLINE = "#E4E4E4"
OUTN = {"correct": "genotyped correctly", "wrong": "found, wrong genotype", "partial": "wrong extent",
        "missed": "missed"}

# example selection rule (thresholds below)
FLANK_READ = 1000          # bp of aligned flank a read must cover on both sides of the element
MIN_READS = 4              # reads carrying and reads skipping the element, each
LTR_LEN_TOL = 0.05         # LTR lengths alike: |l5 - l3| <= 5% of the longer
LTR_MIN_ID = 0.95          # LTR identity (global alignment, matches / longer LTR)
JX_FLANK, JX_MISM = 30, 3  # empty-site check: 30 bp each side of the TSD, <= 3 mismatches each, ungapped
WIN_CHECK = 5000           # "no other pooled locus" within this many bp of the element
VIEW_FLANK = 4000          # panel A: flank drawn on each side (clipped to the colinear alignment blocks)
BIG_INDEL = 30             # alignment indels >= this split the ribbons
JX_HALF = (16, 6)          # panel B: bp drawn outside | inside the element at each junction



# ---------------------------------------------------------------------------------------------------- helpers
def log(msg, verbose=True):
    if verbose:
        print(msg, file=sys.stderr, flush=True)


def die(msg):
    sys.exit(f"error: {msg}")


def rd(path):
    if not os.path.isfile(path):
        die(f"table not found: {path}")
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def pct(x, nd=0):
    """Round half up (table values carry 6 significant digits)."""
    q = Decimal("1") if nd == 0 else Decimal("1." + "0" * nd)
    return str(Decimal(str(100 * float(x))).quantize(q, rounding=ROUND_HALF_UP))


NUMS = []
LGREL = ("report/", "calls/", "truth/", "pool/")


def num(panel, item, value, source, ds="", tool=""):
    """Record a number with its source; paths relative to ltrkin_geno/ get that prefix (bench-relative otherwise)."""
    if source.startswith(LGREL):
        source = "ltrkin_geno/" + source
    source = source.replace("score.score_dataset", "ltrkin_geno/scripts/score.py score_dataset")
    NUMS.append(dict(panel=panel, item=item, dataset=ds, caller=tool, value=value, source=source))
    return value


def wilson(k, n, z=1.959964):
    if n == 0:
        return float("nan"), float("nan")
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - h), min(1.0, c + h)


# ---------------------------------------------------------------------------------------------------- data
def load_tables():
    T = dict(yield_=rd(TAB / "ltrkin_yield.tsv"), recall=rd(TAB / "recall.tsv"), subsets=rd(CHK / "subsets.tsv"),
             ends=rd(CHK / "ends.tsv"), gap=rd(CHK / "gap.tsv"))
    T["Y"] = {(r["ds"], r["tool"], r["run"]): r for r in T["yield_"]}
    T["R"] = {(r["ds"], r["tool"]): r for r in T["recall"]}
    T["S"] = {(r["ds"], r["subset"], r["metric"], r["tool"]): r for r in T["subsets"]}
    T["E"] = {(r["ds"], r["tool"], r["end_bin"]): r for r in T["ends"]}
    return T


def running_set(T):
    """Runs still in progress, taken from the table, not from run markers."""
    return {(r["ds"], r["tool"]) for r in T["yield_"] if r["status"] == "DNF (running)"}


def cat_of(st, o):
    if st[0] in ("missed", "partial"):
        return st[0]
    return "correct" if o[0] == "F/E" else "wrong"


def end_bin4(err):
    """lead_checks.end_bin: exact (microhomology window), le5, le20, gt20; no_true_ends."""
    if err is None:
        return "no_true_ends"
    if err["exact"]:
        return "exact"
    w = max(abs(err["left"]), abs(err["right"]))
    return "le5" if w <= 5 else "le20" if w <= 20 else "gt20"


def score_all(T, running, verbose):
    S = {}
    bad = []
    for ds in DS:
        fai = str(lg_config.DATASETS[ds]["ref"]) + ".fai"
        pool, _, res, ok, n = score.score_dataset(str(LG), ds, fai, running)
        prim = [t for t in ok if t in PRIMARY]
        for t in prim:
            outs = {run: res[t]["outcome"][run] for run in ("e2e", "control")}
            for run in ("e2e", "control"):
                fe = sum(o[0] == "F/E" for o in outs[run])
                y = T["Y"].get((ds, t, run))
                if y is None or y["status"] != "ok" or round(float(y["yield"]) * int(y["n_loci"])) != fe \
                        or int(y["n_loci"]) != n:
                    bad.append(f"{ds} {t} {run}: F/E {fe} vs ltrkin_yield.tsv {y and y['yield']} x {n}")
            r = T["R"][(ds, t)]
            found = sum(s[0] == "found" for s in res[t]["status_rows"])
            nex = sum(e is not None and e["exact"] for e in res[t]["errs"])
            if (found, nex) != (int(r["found"]), int(r["n_exact"])):
                bad.append(f"{ds} {t}: found/exact {found}/{nex} vs recall.tsv {r['found']}/{r['n_exact']}")
            cnt = collections.Counter()
            for (st, _), err, o in zip(res[t]["status_rows"], res[t]["errs"], outs["e2e"]):
                if st == "found":
                    b = end_bin4(err)
                    cnt[(b, "n")] += 1
                    cnt[(b, "fe")] += o[0] == "F/E"
            for b in ("exact", "le5", "le20", "gt20", "no_true_ends"):
                e = T["E"][(ds, t, b)]
                if (cnt[(b, "n")], cnt[(b, "fe")]) != (int(e["n_found"]), int(e["n_fe"])):
                    bad.append(f"{ds} {t} {b}: {cnt[(b, 'n')]}/{cnt[(b, 'fe')]} vs ends.tsv {e['n_found']}/{e['n_fe']}")
        S[ds] = dict(pool=pool, res=res, ok=prim, n=n,
                     cats={t: [cat_of(s, o) for s, o in zip(res[t]["status_rows"], res[t]["outcome"]["e2e"])]
                           for t in prim})
        log(f"{ds}: {n} pooled loci; callers scored: {', '.join(prim)}", verbose)
    if bad:
        die("self-check failed, the per-locus outcomes do not reproduce the report tables:\n  " + "\n  ".join(bad))
    log("self-check passed: F/E (e2e, control), found, exact and end-error bins equal ltrkin_yield.tsv, recall.tsv "
        "and ends.tsv for every scored caller", verbose)
    return S


def caller_order(T, S):
    """The shared tool order of the publication set; fails if the scored callers differ."""
    if set(ORDER) != set(PRIMARY):
        die(f"caller set differs from the shared order: {sorted(PRIMARY)} vs {sorted(ORDER)}")
    return list(ORDER)


# ---------------------------------------------------------------------------------------------------- example
def read_classes(bam, chrom, ts, te, flank=FLANK_READ, tol=30):
    """MAPQ>=20 primary/supplementary alignments around [ts, te): reads carrying the element (one alignment over
    [ts-flank, te+flank] without a deletion of >= half the element) and reads skipping it (a deletion of >= half the
    element inside one such alignment, or a split read whose parts end/start within tol + TSD of ts / te and cover
    the flanks). Returns {read: (kind, [(ref_start, ref_end), ...] sorted)} for informative spanning reads."""
    segs = collections.defaultdict(list)
    for a in bam.fetch(chrom, max(0, ts - flank), te + flank):
        if a.is_secondary or a.is_unmapped or a.mapping_quality < 20:
            continue
        segs[a.query_name].append(a)
    L = te - ts
    out = {}
    for q, al in segs.items():
        kind, parts = None, None
        for a in al:
            if a.reference_start <= ts - flank and a.reference_end >= te + flank:
                pos, dels = a.reference_start, []
                for op, ln in a.cigartuples:
                    if op in (2, 3) and ln >= 0.5 * L and pos < te and pos + ln > ts:
                        dels.append((pos, pos + ln))
                    if op in (0, 2, 3, 7, 8):
                        pos += ln
                if dels:
                    d0, d1 = dels[0]
                    kind, parts = "skip", [(a.reference_start, d0), (d1, a.reference_end)]
                else:
                    kind, parts = "carry", [(a.reference_start, a.reference_end)]
        if kind is None and len(al) >= 2:
            left = [a for a in al if a.reference_start <= ts - flank and abs(a.reference_end - ts) <= tol]
            right = [a for a in al if a.reference_end >= te + flank and abs(a.reference_start - te) <= tol]
            if left and right and left[0].is_reverse == right[0].is_reverse:
                kind, parts = "skip", [(left[0].reference_start, left[0].reference_end),
                                       (right[0].reference_start, right[0].reference_end)]
        if kind:
            out[q] = (kind, sorted(parts))
    return out


def ltr_identity(a, b):
    r = parasail.nw_stats_scan_16(a, b, 5, 1, parasail.dnafull)
    return r.matches / max(len(a), len(b))


def empty_site(oth, o, left30, tsd, right30):
    """Offsets in the other haplotype's window (oriented like the reference) where flank30 + TSD + flank30 align
    ungapped with the TSD exact and <= JX_MISM mismatches per flank; returns (hits, segment, segment start, strand)."""
    oc, os_, ol, orr = o
    s0 = max(0, ol - 300)
    seg = oth.fetch(oc, s0, orr + 300).upper()
    if os_ == "-":
        seg = revcomp(seg)
    d, hits = len(tsd), []
    for k in range(0, len(seg) - 2 * JX_FLANK - d + 1):
        if seg[k + JX_FLANK:k + JX_FLANK + d] != tsd:
            continue
        ml = sum(x != y for x, y in zip(seg[k:k + JX_FLANK], left30))
        mr = sum(x != y for x, y in zip(seg[k + JX_FLANK + d:k + 2 * JX_FLANK + d], right30))
        if ml <= JX_MISM and mr <= JX_MISM:
            hits.append(k)
    return hits, seg, s0, os_


def select_example(S, verbose):
    """The selection rule on haplotype 2; returns the chosen locus and the funnel."""
    ds = DS[0]
    D = lg_config.DATASETS[ds]
    pool, res, cats, ok = S[ds]["pool"], S[ds]["res"], S[ds]["cats"], S[ds]["ok"]
    mode = {t: collections.Counter(cats[t]).most_common(1)[0][0] for t in ok}
    lens = [p["match_end"] - p["match_start"] for p in pool]
    med = statistics.median(lens)
    ref, oth = pysam.FastaFile(str(D["ref"])), pysam.FastaFile(str(D["other"]))
    bam = pysam.AlignmentFile(str(D["bam"]))
    lq_calls = read_calls(final_calls_path(str(LG / "calls" / ds), "ltr_peel"), "ltr_peel")
    lq_truth = {r["element_id"]: r for r in read_truth(str(LG / "truth" / ds / "ltr_peel.tsv"), "ltr_peel")}
    rows = []
    for i, p in enumerate(pool):
        if p["basis"] != "true_confirmed" or p["te_flag"] != "ok":
            continue
        ch, ts, te, d, w = p["chrom"], p["true_start"], p["true_end"], p["dup_len"], p["window"]
        tsd = (p["dup_seq"] or "").upper()
        r = dict(i=i, locus=p["locus_id"], chrom=ch, ts=ts, te=te, len=te - ts, d=d, w=w, tsd=tsd)
        r["tgca"] = ref.fetch(ch, ts, ts + 2).upper() == "TG" and ref.fetch(ch, te - 2, te).upper() == "CA"
        r["tgca_win"] = any(ref.fetch(ch, ts - k, ts - k + 2).upper() == "TG"
                            and ref.fetch(ch, te + w - k - 2, te + w - k).upper() == "CA" for k in range(w + 1))
        r["tsd_both"] = ref.fetch(ch, ts - d, ts).upper() == tsd == ref.fetch(ch, te, te + d).upper()
        st = res["ltr_peel"]["status_rows"][i] if "ltr_peel" in res and "status_rows" in res["ltr_peel"] else None
        r["ltr_ok"], r["empty"], r["o"] = False, False, None
        if st and st[0] == "found":
            c = lq_calls[st[1]]
            r["lq_call"] = c
            if c["lltr_end"] and c["rltr_start"] and c["rltr_start"] > c["lltr_end"]:
                l5, l3 = c["lltr_end"] - c["start"], c["end"] - c["rltr_start"]
                ida = ltr_identity(ref.fetch(ch, c["start"], c["lltr_end"]).upper(),
                                   ref.fetch(ch, c["rltr_start"], c["end"]).upper())
                r.update(l5=l5, l3=l3, ltr_id=ida)
                r["ltr_ok"] = l5 > 0 and l3 > 0 and abs(l5 - l3) <= LTR_LEN_TOL * max(l5, l3) and ida >= LTR_MIN_ID
            tr = lq_truth.get(gff_id("ltr_peel", st[1]))
            if tr and tr["o_left"] is not None and tr["o_right"] is not None:
                r["o"] = (tr["o_chrom"], tr["o_strand"], tr["o_left"], tr["o_right"])
                hits, _, _, _ = empty_site(oth, r["o"], ref.fetch(ch, ts - d - JX_FLANK, ts - d).upper(), tsd,
                                           ref.fetch(ch, te + d, te + d + JX_FLANK).upper())
                r["empty"] = len(hits) == 1
        r["typical"] = all(cats[t][i] == mode[t] for t in ok)
        r["others"] = sum(1 for q in pool if q is not p and q["chrom"] == ch and q["match_start"] < te + WIN_CHECK
                          and q["match_end"] > ts - WIN_CHECK)
        r["carry"] = r["skip"] = None
        if r["tgca_win"] and r["tsd_both"] and r["ltr_ok"] and r["empty"]:
            rc = read_classes(bam, ch, ts, te)
            r["carry"] = sum(k == "carry" for k, _ in rc.values())
            r["skip"] = sum(k == "skip" for k, _ in rc.values())
        rows.append(r)
    reads_ok = lambda r: r["carry"] is not None and r["carry"] >= MIN_READS and r["skip"] >= MIN_READS  # noqa: E731
    funnel = [("haplotype 2 pooled loci", len(pool)), ("true ends known, exact TSD (true_confirmed, te_flag ok)", len(rows))]
    f = [r for r in rows if r["tgca"] and r["tsd_both"]]
    funnel.append(("TG...CA at the true ends", len(f)))
    f = [r for r in f if r["ltr_ok"]]
    funnel.append(("both LTRs, lengths within 5%, identity >= 0.95 (LTR_peel's LTRs)", len(f)))
    f = [r for r in f if r["empty"]]
    funnel.append(("TSD once at the haplotype 1 empty site", len(f)))
    f = [r for r in f if reads_ok(r)]
    funnel.append(("reads: >= 4 carrying and >= 4 skipping, >= 1 kb flank each side", len(f)))
    n_pre_typ = len(f)
    f = [r for r in f if r["typical"]]
    funnel.append(("typical: every finished caller at its most frequent outcome", len(f)))
    n_typ = len(f)
    f_clean = [r for r in f if r["others"] == 0]
    funnel.append((f"no other pooled locus within {WIN_CHECK} bp", len(f_clean)))
    elig = sorted(f_clean or f, key=lambda r: (abs(r["len"] - med), r["locus"]))
    if not elig:
        die("no haplotype 2 locus passes the example rule")
    relaxed = sorted([r for r in rows if r["tgca_win"] and r["tsd_both"] and r["ltr_ok"] and r["empty"] and reads_ok(r)
                      and r["typical"] and r["others"] == 0], key=lambda r: (abs(r["len"] - med), r["locus"]))
    for k, v in funnel:
        log(f"  example funnel: {k}: {v}", verbose)
    return dict(chosen=elig[0], eligible=elig, relaxed=relaxed, funnel=funnel, median=med, mode=mode,
                n_pre_typical=n_pre_typ, n_typical=n_typ, n_excluded_by_others=n_typ - len(f_clean))


def paf_blocks(paf, qchrom, tchrom, q0, q1):
    """Aligned blocks [(q_start, q_end, t_start, t_end)] of the + strand PAF records qchrom->tchrom overlapping
    [q0, q1), split at indels >= BIG_INDEL (smaller indels and mismatches stay inside a block). 0-based half-open."""
    blocks = []
    with gzip.open(paf, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] != qchrom or f[5] != tchrom or int(f[3]) <= q0 or int(f[2]) >= q1:
                continue
            if f[4] != "+":
                die(f"{paf}: minus-strand record at the example; not drawn by this script")
            cg = next(x[5:] for x in f[12:] if x.startswith("cg:Z:"))
            q, t = int(f[2]), int(f[7])
            bq, bt = q, t
            for n, op in re.findall(r"(\d+)([=XIDM])", cg):
                n = int(n)
                if op in "=XM":
                    q += n
                    t += n
                elif op == "I":
                    if n >= BIG_INDEL:
                        blocks.append((bq, q, bt, t))
                        bq, bt = q + n, t
                    q += n
                else:
                    if n >= BIG_INDEL:
                        blocks.append((bq, q, bt, t))
                        bq, bt = q, t + n
                    t += n
            blocks.append((bq, q, bt, t))
    return sorted(b for b in blocks if b[1] > q0 and b[0] < q1)


def example_data(S, ex, verbose):
    ds = DS[0]
    D = lg_config.DATASETS[ds]
    r = ex["chosen"]
    i, ch, ts, te, d = r["i"], r["chrom"], r["ts"], r["te"], r["d"]
    ref, oth = pysam.FastaFile(str(D["ref"])), pysam.FastaFile(str(D["other"]))
    oc, ostr, _, _ = r["o"]
    blocks = paf_blocks(str(D["paf"]), ch, oc, ts - 50_000, te + 50_000)
    # the colinear stretch around the element: the block ending at ts and the chain of blocks from te
    lb = [b for b in blocks if b[1] == ts or (b[0] < ts <= b[1])]
    rb = [b for b in blocks if b[0] == te or (b[0] <= te < b[1])]
    if not lb or not rb:
        die("example: no alignment block ends at the element's left end or starts at its right end")
    right_end = rb[0][1]
    for b in blocks:                                   # extend through blocks split only by indels (contiguous in q)
        if b[0] == right_end:
            right_end = b[1]
    left_start = lb[0][0]
    for b in sorted(blocks, reverse=True):
        if b[1] == left_start:
            left_start = b[0]
    fl = min(VIEW_FLANK, ts - left_start, right_end - te)
    W0, W1 = ts - fl, te + fl
    vblocks = [(max(b[0], W0), min(b[1], W1), b[2] + max(0, W0 - b[0]), b[3] - max(0, b[1] - W1))
               for b in blocks if b[1] > W0 and b[0] < W1]
    gaps1, gaps2 = [], []                              # indels >= BIG_INDEL in the drawn flanks (not the element)
    vb = sorted(vblocks)
    for u, w in zip(vb, vb[1:]):
        qg, tg = w[0] - u[1], w[2] - u[3]
        if u[1] <= ts and w[0] >= te:
            continue
        if tg >= BIG_INDEL and qg < BIG_INDEL:
            gaps1.append((u[1], tg))
        if qg >= BIG_INDEL and tg < BIG_INDEL:
            gaps2.append((u[1], qg))
    # haplotype 1 TSD (empty site) from the junction check
    hits, seg, s0, os_ = empty_site(oth, r["o"], ref.fetch(ch, ts - d - JX_FLANK, ts - d).upper(), r["tsd"],
                                    ref.fetch(ch, te + d, te + d + JX_FLANK).upper())
    if os_ != "+":
        die("example: the empty site lies on the minus strand of haplotype 1; not drawn by this script")
    h_tsd = s0 + hits[0] + JX_FLANK                     # 0-based start of the single TSD on haplotype 1
    # calls of every caller in the window, the matched call flagged
    calls = {}
    for t in PRIMARY:
        if t not in S[ds]["ok"]:
            calls[t] = None
            continue
        cl = read_calls(final_calls_path(str(LG / "calls" / ds), t), t)
        st = S[ds]["res"][t]["status_rows"][i]
        cc = []
        for c in cl:
            if c["chrom"] == ch and c["start"] < W1 and c["end"] > W0:
                cc.append(dict(c, matched=(st[0] == "found" and c["idx"] == st[1])))
        calls[t] = dict(calls=cc, status=st[0], cat=S[ds]["cats"][t][i],
                        geno=S[ds]["res"][t]["outcome"]["e2e"][i], err=S[ds]["res"][t]["errs"][i])
    reads = read_classes(pysam.AlignmentFile(str(D["bam"])), ch, ts, te)
    lq = r["lq_call"]
    E = dict(r=r, W0=W0, W1=W1, fl=fl, blocks=vblocks, h_tsd=h_tsd, ochrom=oc, calls=calls, reads=reads,
             gaps1=gaps1, gaps2=gaps2,
             l5=lq["lltr_end"] - lq["start"], l3=lq["end"] - lq["rltr_start"], strand=lq["strand"], cls=lq["cls"],
             seq2=lambda a, b: ref.fetch(ch, a, b).upper(), seq1=lambda a, b: oth.fetch(oc, a, b).upper())
    log(f"example {r['locus']} {ch}:{ts}-{te} (0-based) len {te - ts} TSD {r['tsd']}; window {W0}-{W1}; "
        f"haplotype 1 TSD at {oc}:{h_tsd}; reads carry {sum(k == 'carry' for k, _ in reads.values())} "
        f"skip {sum(k == 'skip' for k, _ in reads.values())}", verbose)
    return E


# ---------------------------------------------------------------------------------------------------- breakdown data
DIMS = {"nesting": ("nesting", ["unnested", "tenant", "host"]),
        "protein": ("protein_class", ["coding", "partial", "proteinless"])}
CLSN = {"unnested": "non-nested", "tenant": "tenant", "host": "host", "coding": "coding", "partial": "partly coding",
        "proteinless": "protein-less"}      # data key "partial" (ltrkin_breakdown.tsv cls) is named "partly coding"
CATS = ["correct", "wrong", "partial", "missed"]
MIN_CI = 20                # classes with fewer loci: counts only


def load_breakdown(S, verbose):
    """experiments/figpub/ltrkin_breakdown.tsv (per ds x tool x class: counts, shares, block-bootstrap CIs) and
    ltrkin_locus_class.tsv (per locus: class, nominators, LTR_peel-only flag, every caller's outcome). Self-checks
    (fatal): the class table's per-locus outcomes equal score.score_dataset's for every finished caller; its per-class
    counts equal the breakdown table's; every finished caller's class counts sum to all loci; every nominator is a
    finished caller; the LTR_peel-only flag equals nominators == ltr_peel."""
    rows = rd(FIGPUB / "ltrkin_breakdown.tsv")
    lc = rd(FIGPUB / "ltrkin_locus_class.tsv")
    B = {(r["ds"], r["tool"], r["dimension"], r["cls"]): r for r in rows if r["role"] == "primary"}
    bad = []
    LC = collections.defaultdict(list)
    for r in lc:
        LC[r["ds"]].append(r)
    for ds in DS:
        pool = S[ds]["pool"]
        loci = LC[ds]
        if [p["locus_id"] for p in pool] != [r["locus_id"] for r in loci]:
            bad.append(f"{ds}: ltrkin_locus_class.tsv loci differ from the pool (order or ids)")
            continue
        for t in S[ds]["ok"]:
            if [r[f"out_{t}"] for r in loci] != S[ds]["cats"][t]:
                bad.append(f"{ds} {t}: per-locus outcomes in ltrkin_locus_class.tsv differ from score_dataset")
        for r in loci:                                 # hatching rule: own loci = nominators exactly one caller
            noms = r["nominators"].split(",") if r["nominators"] else []
            if not noms or any(x not in S[ds]["ok"] for x in noms):
                bad.append(f"{ds} {r['locus_id']}: nominators '{r['nominators']}' empty or not a finished caller")
            if (r["ltr_peel_only"] == "1") != (r["nominators"] == "ltr_peel"):
                bad.append(f"{ds} {r['locus_id']}: ltr_peel_only {r['ltr_peel_only']} vs nominators {r['nominators']}")
        for dim, (col, classes) in DIMS.items():
            for t in ORDER:
                tot = 0
                for c in classes:
                    b = B.get((ds, t, dim, c))
                    sub = [r for r in loci if r[col] == c]
                    if b is None or int(b["n_loci"]) != len(sub):
                        bad.append(f"{ds} {t} {dim} {c}: n_loci {b and b['n_loci']} vs {len(sub)}")
                        continue
                    if t not in S[ds]["ok"]:
                        if b["status"] == "ok":
                            bad.append(f"{ds} {t}: unfinished here but status ok in the breakdown table")
                        continue
                    for k in CATS:
                        n_k = sum(r[f"out_{t}"] == k for r in sub)
                        if int(b[f"n_{k}"]) != n_k:
                            bad.append(f"{ds} {t} {dim} {c} {k}: {b[f'n_{k}']} vs {n_k}")
                    tot += len(sub)
                if t in S[ds]["ok"] and tot != S[ds]["n"]:
                    bad.append(f"{ds} {t} {dim}: classes sum to {tot}, not {S[ds]['n']}")
    if bad:
        die("breakdown self-check failed:\n  " + "\n  ".join(bad))
    log("breakdown self-check passed: class table outcomes = score_dataset; class counts = ltrkin_breakdown.tsv; "
        "nominators are finished callers; ltr_peel_only = (nominators == ltr_peel)", verbose)
    return B, LC


def own_loci_raw(LC, ds, t, col=None, c=None):
    """Hatching rule: (loci whose only nominator is caller t, i.e. nominated by this caller alone,
    of them genotyped correctly from t's call); one class (col == c) or all loci (col None). Nominators = callers with
    a member call at the pooled locus (ltrkin_geno lead_checks.nominators)."""
    sub = [r for r in LC[ds] if r["nominators"] == t and (col is None or r[col] == c)]
    return len(sub), sum(r[f"out_{t}"] == "correct" for r in sub)


def own_loci(OWN, ds, t, col=None, c=None):
    """own_loci_raw's counts as stored in the source data (CD_own_loci.tsv)."""
    return OWN[(ds, t, col or "", c or "")]




# ---------------------------------------------------------------------------------------------------- A: reported features
def gff_rows(path, seqid, w0, w1):
    """GFF3 rows on seqid overlapping [w0, w1) (0-based); GFF3 is 1-based closed -> start - 1, end."""
    out = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9 or c[0] != seqid:
                continue
            a, b = int(c[3]) - 1, int(c[4])
            if b > w0 and a < w1:
                att = dict(x.split("=", 1) for x in c[8].strip().rstrip(";").split(";") if "=" in x)
                out.append(dict(type=c[2], s=a, e=b, strand=c[6], att={k.strip(): v for k, v in att.items()}))
    return out


def reported_features(E, ds=DS[0]):
    """Per caller, every call in panel A's window with what the caller itself reports (strand, LTR bounds, TSDs,
    protein domains), from its call table (ltrkin_geno/calls, converted from the native output) or its native output
    (ltrkin_geno/runs/<ds>/<tool>/); a feature the caller does not report stays None / empty. 0-based half-open.
    TSD: 'tsd' = intervals to draw in A (sequence-only reports drawn as the adjacent flanks, fig_benchmark A's rule);
    'tsd_bp' = intervals placeable at base pairs in B (reported coordinates, or a reported sequence that equals the
    adjacent flank). Returns ({tool: [feat] or None}, {tool: {feature: source}})."""
    ch, W0, W1 = E["r"]["chrom"], E["W0"], E["W1"]
    nm = rd(LG / "data" / ds / "names.tsv")
    seq = {r["original_id"]: r["new_id"] for r in nm}[ch]               # renamed ID in the native outputs (c0006)
    seqn = f"seq{[r['new_id'] for r in nm].index(seq)}"                  # LTRharvest's seqN = Nth FASTA record
    run = LG / "runs" / ds
    ref = pysam.FastaFile(str(lg_config.DATASETS[ds]["ref"]))
    raw_hr = {c["id"]: c for c in read_calls(str(LG / "calls" / ds / "hirepeat.tsv"), "hirepeat")}
    F, SRC = {}, {}
    for t in ORDER:
        info = E["calls"][t]
        if info is None:
            F[t] = None
            continue
        F[t] = []
        for c in info["calls"]:
            f = dict(s=c["start"], e=c["end"], matched=c["matched"], id=c["id"],
                     strand=c["strand"] if c["strand"] in ("+", "-") else None, lltr=c["lltr_end"],
                     rltr=c["rltr_start"], tsd_seq=c["tsd"], tsd=[], tsd_bp=[], dom=[])
            if t == "hirepeat":                       # LTR edges in hirepeat.ltr.tsv come from Kmer2LTR, not the tool
                r0 = raw_hr[c["id"]]
                if (r0["start"], r0["end"]) != (c["start"], c["end"]) or r0["lltr_end"] is not None:
                    die(f"HiRepEat call {c['id']}: raw call table disagrees with hirepeat.ltr.tsv")
                f["lltr"] = f["rltr"] = None
            F[t].append(f)

    def adjacent(f):
        n = len(f["tsd_seq"])
        fl = [(f["s"] - n, f["s"]), (f["e"], f["e"] + n)]
        f["tsd"] = fl
        seqs = [ref.fetch(ch, a, b).upper() for a, b in fl]
        if all(x in (f["tsd_seq"], revcomp(f["tsd_seq"])) for x in seqs):
            f["tsd_bp"] = fl
    # LTR_peel: strand, LTRs, TSD sequence (calls table <- clean_ltr.tsv, offsets 0,0); domains = protein_match children
    g = gff_rows(run / "ltr_peel" / f"{ds}_LTRs_all_depth_LTR_cleaned.gff3", seq, W0, W1)
    ids = {(r["s"], r["e"]): r["att"]["ID"] for r in g if r["type"] == "LTR_retrotransposon"}
    for f in F["ltr_peel"] or []:
        if f["tsd_seq"]:
            adjacent(f)
        pid = ids.get((f["s"], f["e"]))
        f["dom"] = [(r["s"], r["e"], r["att"].get("Name", "")) for r in g
                    if r["type"] == "protein_match" and r["att"].get("Parent") == pid]
    SRC["ltr_peel"] = dict(strand="calls table (clean_ltr.tsv strand)", LTR="calls table lltr_end / rltr_start",
                           TSD="calls table tsd sequence (clean_ltr.tsv, offset 0,0) -> adjacent flanks",
                           domains=f"runs/{ds}/ltr_peel/{ds}_LTRs_all_depth_LTR_cleaned.gff3 protein_match (1-based)")
    # LTR_retriever: strand, LTRs (calls table); TSD coordinates from the pass list (5_TSD, 3_TSD; 1-based closed)
    pl = {}
    with open(run / "ltr_retriever" / f"{ds}.fa.mod.pass.list") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.split("\t")
            m = re.match(r"(\S+):(\d+)\.\.(\d+)", c[0])
            if not m or m.group(1) != seq:
                continue
            a, b = sorted((int(m.group(2)), int(m.group(3))))
            t5, t3 = c[4], c[5]
            if t5.startswith("NA") or t3.startswith("NA"):
                continue
            iv = [tuple(int(x) for x in v.split("..")) for v in (t5, t3)]
            pl[(a - 1, b)] = sorted((min(u) - 1, max(u)) for u in iv)
    for f in F["ltr_retriever"] or []:
        f["tsd"] = f["tsd_bp"] = pl.get((f["s"], f["e"]), [])
    SRC["ltr_retriever"] = dict(strand="calls table (pass.list)", LTR="calls table", domains="none",
                                TSD=f"runs/{ds}/ltr_retriever/{ds}.fa.mod.pass.list 5_TSD / 3_TSD (1-based)")
    # LTRharvest: no strand ('?'); LTRs (calls table); TSD features in its GFF3 (seqN = Nth FASTA record)
    g = gff_rows(run / "ltrharvest" / "ltrharvest.gff3", seqn, W0 - 100, W1 + 100)
    for f in F["ltrharvest"] or []:
        el = [r for r in g if r["type"] == "LTR_retrotransposon" and (r["s"], r["e"]) == (f["s"], f["e"])]
        if el:
            par = el[0]["att"]["Parent"]
            f["tsd"] = f["tsd_bp"] = sorted((r["s"], r["e"]) for r in g if r["type"] == "target_site_duplication"
                                            and r["att"].get("Parent") == par)
    SRC["ltrharvest"] = dict(strand="not reported ('?')", LTR="calls table", domains="none",
                             TSD=f"runs/{ds}/ltrharvest/ltrharvest.gff3 target_site_duplication ({seqn}; 1-based)")
    # LTR_Finder: strand, LTRs, TSR sequence (calls table <- scn); domains from scn RT / IN / RH columns
    scn = {}
    with open(run / "ltr_finder" / f"{ds}.fa.finder.combine.scn") as fh:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 13 or c[1] != seq or "-" not in c[2]:
                continue
            a, b = (int(x) for x in c[2].split("-"))
            dom = []
            for name, v in zip(("RT", "IN (core)", "IN (c-term)", "RH"), c[8:12]):
                if v != "N-N":
                    x, y = (int(z) for z in v.split("-"))
                    dom.append((min(x, y) - 1, max(x, y), name))
            scn[(a - 1, b)] = dom
    for f in F["ltr_finder"] or []:
        if f["tsd_seq"]:
            adjacent(f)
        f["dom"] = scn.get((f["s"], f["e"]), [])
    SRC["ltr_finder"] = dict(strand="calls table (scn Strand)", LTR="calls table", TSD="scn TSR -> adjacent flanks",
                             domains=f"runs/{ds}/ltr_finder/{ds}.fa.finder.combine.scn RT / IN / RH columns")
    # HiRepEat: span only; HiTE: LTRs and TSD sequence (coordinates NA), no strand
    for f in F["hite"] or []:
        if f["tsd_seq"]:
            adjacent(f)
    SRC["hirepeat"] = dict(strand="not reported", LTR="not reported (Kmer2LTR edges are a benchmark instrument)",
                           TSD="not reported", domains="none")
    SRC["hite"] = dict(strand="not reported ('.')", LTR="calls table (intact_LTR.list)", domains="none",
                       TSD="calls table tsd sequence (intact_LTR.list, coordinates NA) -> adjacent flanks")
    for t in ("ltrpred", "dante_ltr"):            # their domain readers were not needed: no call in the window
        if F[t]:
            die(f"{t} has calls in panel A's window; add its native TSD / domain reader first")
        SRC[t] = dict(strand="-", LTR="-", TSD="-", domains="- (no call in the window)")
    return F, SRC


# ---------------------------------------------------------------------------------------------------- drawing
class Page:
    """Figure in mm from the top-left corner."""

    def __init__(self, w, h):
        self.W, self.H = w, h
        self.fig = plt.figure(figsize=(w * MM, h * MM), dpi=300)   # text measured at the output resolution
        self.can = poafig.mm_axes(self.fig, 0, 0, w, h)
        self.can.set_xlim(0, w)
        self.can.set_ylim(h, 0)
        self.can.axis("off")

    def ax(self, x, y, w, h, xlim=None, ylim=None):
        a = poafig.mm_axes(self.fig, x, y, w, h)
        a.set_xlim(*(xlim or (0, w)))
        a.set_ylim(*(ylim or (y + h, y)))          # y data = page mm, downward
        a.axis("off")
        return a

    def text(self, x, y, s, **kw):
        kw.setdefault("fontsize", 6.0)
        return self.fig.text(x / self.W, 1 - y / self.H, s, **kw)

    def letter(self, x, y, s):
        poafig.panel_letter(self.fig, x, y, s)

    def line(self, pts, **kw):
        poafig.fig_line(self.fig, pts, **kw)

    def width(self, s, fontsize):
        """width (mm) of a text at fontsize"""
        t = self.can.text(0, 0, s, fontsize=fontsize)
        bb = t.get_window_extent(renderer=self.fig.canvas.get_renderer())
        t.remove()
        return bb.width / self.fig.dpi * 25.4


FAILC = "#5A5A5A"          # set-wide missing-result rule: tool-side failure (timeout, memory cap, crash)
HATCH = "//////"           # loci nominated by this caller alone (as fig_benchmark): white lines 4 pt apart along the bar
THIN_MM = 0.3              # hatches narrower than ~2 hatch-line widths (2 x 0.4 pt = 0.28 mm) cannot show a white line
                           # with color beside it (single-pixel specks at 300 dpi): the legend gives their counts
HATCH_MM = {}              # (panel, ds, caller, class) -> (hatched loci, class n, drawn width mm)


def tint(hexc, a):
    """hexc blended with white; a = share of the color (fig_benchmark)"""
    r, g, b = to_rgb(hexc)
    return (1 - a + a * r, 1 - a + a * g, 1 - a + a * b)


def tool_failed(status):
    """Tool-side failure (did not finish within 36 h, exceeded 200 GB, crashed) vs evaluation-side gap."""
    return status.startswith("DNF") or "fail" in status or "timeout" in status


def missing_rect(ax, x0, y0, w, h, tx, ty, status, fs=6):
    """Missing result, set-wide style: tool failed = dark gray fill, white "tool failed"; evaluation gap = light
    hatch, "not evaluated"."""
    if tool_failed(status):
        ax.add_patch(Rectangle((x0, y0), w, h, fc=FAILC, lw=0))
        ax.text(tx, ty, "tool failed", ha="center", va="center", fontsize=fs, color="white")
    else:
        ax.add_patch(Rectangle((x0, y0), w, h, fc="#F2F2F2", ec="none", lw=0.0, hatch="//////",
                               hatchcolor="#B8B8B8"))
        ax.text(tx, ty, "not evaluated", ha="center", va="center", fontsize=fs, color=MUTED)
    return "tool failed" if tool_failed(status) else "not evaluated"


def lab_left(P, x, y, s, **kw):
    kw.setdefault("ha", "right")
    kw.setdefault("va", "center")
    return P.text(x, y, s, **kw)


def tool_lab(P, x, y, t, fontsize=6.0, **kw):
    """A caller's row label, set-wide style (fig_benchmark): dark name (LTR_peel bold) and a 1.3-mm square of the
    tool's color; the square's right edge at x (mm)."""
    kw.setdefault("ha", "right")
    kw.setdefault("va", "center")
    P.can.add_patch(Rectangle((x - 1.3, y - 0.65), 1.3, 1.3, fc=COL[t], ec="none", lw=0))
    return P.text(x - 2.0, y, NAME[t], fontsize=fontsize, color=INK,
                  fontweight="bold" if t == "ltr_peel" else "normal", **kw)


def chevron(ax, x, y, h, strand, col, bp_per_mm):
    d = 0.55 * bp_per_mm * (1 if strand == "+" else -1)
    ax.plot([x - d / 2, x + d / 2, x - d / 2], [y - h * 0.28, y, y + h * 0.28], color=col, lw=0.5,
            solid_capstyle="butt", solid_joinstyle="miter", zorder=5)


def free_center(s0, e0, blocks):
    """center of the widest stretch of [s0, e0) not covered by blocks (for the strand chevron; fig_benchmark)"""
    cuts = sorted((max(a, s0), min(b, e0)) for a, b in blocks if b > s0 and a < e0)
    best, x = (0, (s0 + e0) / 2), s0
    for a, b in cuts + [(e0, e0)]:
        if a - x > best[0]:
            best = (a - x, (x + a) / 2)
        x = max(x, b)
    return best[1]


def hkey(P, x, y, items, fs=6.0, gap=2.6):
    """one-row key on the page canvas: items = (draw(x, y) -> glyph width mm, label)"""
    for draw, lab in items:
        w = draw(x, y)
        P.can.text(x + w + 0.8, y, lab, ha="left", va="center", fontsize=fs)
        x += w + 0.8 + P.width(lab, fs) + gap
    return x


def key_A(P, y, x0):
    """A and B glyph key: what a call can carry (fig_benchmark A) and the read glyphs."""
    c = P.can

    def ltrs(x, yy):
        c.add_patch(Rectangle((x, yy - 0.7), 3.2, 1.4, fc=tint("#7A7A7A", 0.38), ec="none"))
        c.add_patch(Rectangle((x, yy - 0.7), 0.8, 1.4, fc="#7A7A7A", ec="none"))
        c.add_patch(Rectangle((x + 2.4, yy - 0.7), 0.8, 1.4, fc="#7A7A7A", ec="none"))
        return 3.2

    def tsd(x, yy):
        c.add_patch(Rectangle((x, yy - 0.5), 0.5, 1.0, fc="white", ec=INK, lw=0.4))
        return 0.5

    def strand(x, yy):
        c.plot([x, x + 0.6, x], [yy - 0.45, yy, yy + 0.45], color=INK, lw=0.5, solid_joinstyle="miter")
        return 0.6

    def dom(x, yy):
        c.add_patch(Rectangle((x, yy - 0.4), 1.6, 0.8, fc=DOMC, ec="none"))
        return 1.6

    def other(x, yy):
        c.add_patch(Rectangle((x, yy - 0.7), 3.2, 1.4, fc="white", ec="#7A7A7A", lw=0.5))
        return 3.2

    def read(x, yy):
        c.plot([x, x + 3.0], [yy, yy], color=READC, lw=0.9, solid_capstyle="butt")
        return 3.0

    def skip(x, yy):
        c.plot([x, x + 1.0], [yy, yy], color=READC, lw=0.9, solid_capstyle="butt")
        c.plot([x + 1.0, x + 2.6], [yy, yy], color=READC, lw=0.3)
        c.plot([x + 2.6, x + 3.6], [yy, yy], color=READC, lw=0.9, solid_capstyle="butt")
        return 3.6
    return hkey(P, x0, y, [(ltrs, "LTRs"), (tsd, "TSD"), (strand, "strand"), (dom, "protein domain"),
                           (other, "other call"), (read, "HiFi read"), (skip, "read skipping the element")])


def key_C(P, y, x0):
    """C and D key, as fig_benchmark H; the hatch entry ("nominated by this caller alone") is a white
    hatch over the tool-color strip of the first entry (chosen over a hatch on neutral gray, which read as a
    third gray category beside "not genotyped correctly" and "tool failed")."""
    c = P.can

    def strip(hatch=None):
        def f(x, yy):
            for i, t in enumerate(("ltr_peel", "ltr_retriever", "dante_ltr")):
                c.add_patch(Rectangle((x + i * 1.0, yy - 0.75), 1.0, 1.5, fc=COL[t], ec="none", lw=0))
            if hatch:                                  # one hatched patch over the strip: continuous stripes
                c.add_patch(Rectangle((x, yy - 0.75), 3.0, 1.5, fc="none", ec="none", lw=0, hatch=hatch,
                                      hatchcolor="white"))
            return 3.0
        return f

    def box(**st):
        def f(x, yy):
            c.add_patch(Rectangle((x, yy - 0.75), 3.0, 1.5, ec="none", lw=0, **st))
            return 3.0
        return f
    return hkey(P, x0, y, [(strip(), "genotyped correctly (caller's color)"), (strip(HATCH), "nominated by this caller alone"),
                           (box(fc=RESTC), "not genotyped correctly"), (box(fc=FAILC), "tool failed")], gap=3.0)


def draw_A(P, E, FE, S, order, y0, xL, xR, verbose):
    """Panel A: the example locus (kb scale). Returns the bottom y, the x (mm) of the element ends and read counts."""
    r = E["r"]
    ts, te, W0, W1 = r["ts"], r["te"], E["W0"], E["W1"]
    k = (xR - xL) / (W1 - W0)                          # mm per bp
    bppm = 1 / k
    X = lambda b: xL + (b - W0) * k                    # noqa: E731
    c2 = (ts + te) / 2
    off = c2 - (E["h_tsd"] + r["d"] / 2)              # haplotype 1 drawn so the empty site sits over the element centre
    ax = P.ax(xL, y0, xR - xL, 92, xlim=(W0, W1), ylim=(y0 + 92, y0))
    tsd_w = 0.5 * bppm                                 # TSDs (5 bp) drawn 0.5 mm wide so they can be seen
    y_h1, y_h2 = y0 + 6.6, y0 + 17.2
    b_left = [b for b in E["blocks"] if b[1] <= ts + 1]
    b_right = [b for b in E["blocks"] if b[0] >= te - 1]
    h1_lo = min(b[2] for b in b_left)
    h1_hi = max(b[3] for b in b_right)
    ax.plot([h1_lo + off, h1_hi + off], [y_h1, y_h1], color=INK, lw=0.6, solid_capstyle="butt")
    hx = E["h_tsd"] + r["d"] / 2 + off
    ax.add_patch(Rectangle((hx - tsd_w / 2, y_h1 - 0.6), tsd_w, 1.2, fc="white", ec=INK, lw=0.4, zorder=4))
    tick0 = math.ceil(h1_lo / 2000) * 2000
    for t in range(tick0, h1_hi + 1, 2000):
        ax.plot([t + off] * 2, [y_h1 - 0.6, y_h1 - 1.4], color=INK, lw=0.5)
        ax.text(t + off, y_h1 - 1.8, f"{t / 1e6:.3f}", ha="center", va="bottom", fontsize=6)
    lab_left(P, xL - 1.5, y_h1 - 0.1, "haplotype 1", fontsize=6.5)
    lab_left(P, xL - 1.5, y_h1 - 2.9, f"{E['ochrom']} (Mb)", fontsize=6, color=MUTED)
    P.text(xL, y0 + 0.6, "$\\it{Cardamine\\ amara}$", fontsize=6.5, ha="left", va="center")
    # ribbons (aligned flanks) and the insertion wedge (sequence present only on haplotype 2)
    polys = []
    for q0, q1, t0, t1 in E["blocks"]:
        polys.append([(t0 + off, y_h1 + 0.6), (t1 + off, y_h1 + 0.6), (q1, y_h2 - 2.0), (q0, y_h2 - 2.0)])
    ax.add_collection(PolyCollection(polys, fc=RIBBON, ec="none", lw=0))
    ax.add_patch(Polygon([(hx, y_h1 + 0.6), (te, y_h2 - 2.0), (ts, y_h2 - 2.0)], closed=True, fc="#F7F7F7",
                         ec="#A6A6A6", lw=0.4))
    # haplotype 2 line, the true element (LTRs from LTR_peel's call; same LTR ends as LTR_retriever's and LTR_Finder's)
    ax.plot([W0, W1], [y_h2, y_h2], color=INK, lw=0.6, solid_capstyle="butt")
    l5, l3 = E["l5"], E["l3"]
    h = 1.5
    ax.add_patch(Rectangle((ts, y_h2 - h), te - ts, 2 * h, fc=EL_INT, lw=0, zorder=3))
    ax.add_patch(Rectangle((ts, y_h2 - h), l5, 2 * h, fc=EL_LTR, lw=0, zorder=3))
    ax.add_patch(Rectangle((te - l3, y_h2 - h), l3, 2 * h, fc=EL_LTR, lw=0, zorder=3))
    for xb in (ts - tsd_w, te):                        # TSDs, same glyph as the calls' TSDs
        ax.add_patch(Rectangle((xb, y_h2 - h * 0.72), tsd_w, 2 * h * 0.72, fc="white", ec=INK, lw=0.4, zorder=4))
    chevron(ax, c2, y_h2, 2 * h, E["strand"], "white", bppm)
    fam = (E["cls"] or "").split("/")
    ax.text(c2, y_h2 - h - 0.9, f"{(te - ts) / 1000:.1f}-kb {fam[1] if len(fam) >= 2 else 'LTR'} element",
            ha="center", va="bottom", fontsize=6, color=INK)
    lab_left(P, xL - 1.5, y_h2, "haplotype 2", fontsize=6.5)
    num("A", "element length (bp)", te - ts, "report/tables/pool.tsv true_start/true_end (0-based)", DS[0])
    num("A", "element classification as labelled", E["cls"], "calls/camara_hap2/ltr_peel.tsv class (LTR_peel's call)",
        DS[0], "ltr_peel")
    num("A", "left/right LTR lengths as drawn (bp)", f"{l5}/{l3}", "calls/camara_hap2/ltr_peel.tsv lltr_end/rltr_start",
        DS[0], "ltr_peel")
    ya = y_h2 + 2.8
    ax.plot([W0, W1], [ya, ya], color=INK, lw=0.5)
    for t in range(math.ceil(W0 / 2000) * 2000, W1 + 1, 2000):
        ax.plot([t, t], [ya, ya + 0.8], color=INK, lw=0.5)
        if (t - W0) * k > 4.0 and (W1 - t) * k > 4.0:
            ax.text(t, ya + 1.1, f"{t / 1e6:.3f}", ha="center", va="top", fontsize=6)
    lab_left(P, xL - 1.5, ya + 2.3, f"{r['chrom']} (Mb)", fontsize=6, color=MUTED)
    # reads
    yr = ya + 5.6
    pitch = 0.55
    groups = [("carry", "HiFi reads\nwith element"), ("skip", "HiFi reads\nwithout")]
    nread = {}
    for kind, glab in groups:
        rr = sorted([v[1] for v in E["reads"].values() if v[0] == kind], key=lambda s: (s[0][0], s[-1][1]))
        nread[kind] = len(rr)
        ystart = yr
        for parts in rr:
            for s, e in parts:
                s, e = max(s, W0), min(e, W1)
                if e <= s:
                    continue
                if kind == "carry":
                    segs = [(s, min(e, ts), READC), (max(s, ts), min(e, te), EL_READ), (max(s, te), e, READC)]
                    for a, b, col in segs:
                        if b > a:
                            ax.plot([a, b], [yr, yr], color=col, lw=0.8, solid_capstyle="butt")
                else:
                    ax.plot([s, e], [yr, yr], color=READC, lw=0.8, solid_capstyle="butt")
            if kind == "skip":
                for (a0, a1), (b0, b1) in zip(parts, parts[1:]):
                    ax.plot([a1, b0], [yr, yr], color=READC, lw=0.25)
            yr += pitch
        lab_left(P, xL - 1.5, (ystart + yr - pitch) / 2, f"{glab} ({len(rr)})", fontsize=6, linespacing=1.05)
        yr += 1.1
    num("A", "HiFi reads spanning the site with the element", nread["carry"],
        "hap2_d20_s1.bam (pysam; MAPQ>=20, >=1 kb flank both sides)", DS[0])
    num("A", "HiFi reads spanning the site without the element", nread["skip"],
        "hap2_d20_s1.bam (pysam; split or deleted over the element, >=1 kb flank both sides)", DS[0])
    # caller tracks: every call in the window, in the caller's color, with the features it reports
    yc = yr + 1.0
    row_h, bh = 2.3, 1.3
    sub_h = bh + 0.4
    for t in order:
        info, feats = E["calls"][t], FE[t]
        if feats is None:
            missing_rect(ax, W0, yc - 0.25, W1 - W0, bh + 0.5, c2, yc + bh / 2, S[DS[0]]["res"][t]["status"])
            tool_lab(P, xL - 1.5, yc + bh / 2, t)
            yc += row_h
            continue
        ax.plot([W0, W1], [yc + bh / 2, yc + bh / 2], color=ROWLINE, lw=0.4, zorder=0)
        cc = sorted(feats, key=lambda f: (not f["matched"], f["s"]))
        rows = poafig.stack_rows(np.array([f["s"] for f in cc]), np.array([f["e"] for f in cc]), 1.5 * bppm) \
            if cc else np.array([], dtype=int)
        if any(f["matched"] for f in cc):              # the call at the locus on the first sub-row
            mrow = rows[[f["matched"] for f in cc].index(True)]
            rows = np.where(rows == mrow, -1, rows)
            rows = np.where(rows < mrow, rows + 1, rows)
            rows = np.where(rows == -1, 0, rows)
        nsub = int(rows.max()) + 1 if len(rows) else 1
        for f, rw in zip(cc, rows):
            yy = yc + rw * sub_h + bh / 2
            if not f["matched"]:                       # a call not at the locus: outline only (key "other call")
                ax.add_patch(Rectangle((f["s"], yy - bh / 2), f["e"] - f["s"], bh, fc="white", ec=COL[t], lw=0.5,
                                       zorder=2))
                num("A", f"other call {f['s']}-{f['e']} (0-based), drawn as outline", "not at the locus",
                    f"calls/{DS[0]}/{'hirepeat.tsv' if t == 'hirepeat' else t + '.tsv'}", DS[0], t)
                continue
            body = tint(COL[t], 0.38) if f["lltr"] is not None else tint(COL[t], 0.75)
            ax.add_patch(Rectangle((f["s"], yy - bh / 2), f["e"] - f["s"], bh, fc=body, lw=0, zorder=2))
            for a, b in ((f["s"], f["lltr"]), (f["rltr"], f["e"])):
                if a is not None and b is not None and b > a:
                    ax.add_patch(Rectangle((a, yy - bh / 2), b - a, bh, fc=COL[t], lw=0, zorder=2.1))
            for a, b, _ in f["dom"]:
                ax.add_patch(Rectangle((a, yy - bh * 0.28), b - a, bh * 0.56, fc=DOMC, ec="none", zorder=3))
            for a, b in f["tsd"]:
                xa = f["s"] - tsd_w if b <= f["s"] + 1 else f["e"]
                ax.add_patch(Rectangle((xa, yy - bh * 0.36), tsd_w, bh * 0.72, fc="white", ec=INK, lw=0.4,
                                       zorder=4))
            if f["strand"]:
                lo = f["lltr"] if f["lltr"] is not None else f["s"]
                hi = f["rltr"] if f["rltr"] is not None else f["e"]
                chevron(ax, free_center(lo, hi, [(a, b) for a, b, _ in f["dom"]]), yy, bh, f["strand"], INK, bppm)
            num("A", f"call {f['s']}-{f['e']} (0-based){' at the locus' if f['matched'] else ''}: strand / LTR ends / "
                     f"TSDs / protein domains",
                f"{f['strand'] or 'not reported'} / "
                f"{'not reported' if f['lltr'] is None else str(f['lltr']) + ',' + str(f['rltr'])} / "
                f"{f['tsd'] or 'not reported'} / {[(a, b, n) for a, b, n in f['dom']] or 'none'}",
                "; ".join(f"{kk}: {vv}" for kk, vv in FE['_src'][t].items()), DS[0], t)
        tool_lab(P, xL - 1.5, yc + bh / 2, t)
        num("A", "outcome at the example locus",
            OUTN[info["cat"]] + (f" ({info['geno'][0]})" if info["status"] == "found" else ""),
            "score.score_dataset per-locus (status_rows, outcome e2e)", DS[0], t)
        if info["err"] is not None:
            num("B", "called end minus true end, left/right (bp)", f"{info['err']['left']}/{info['err']['right']}",
                "score.score_dataset per-locus errs (calls TSV vs pool true ends)", DS[0], t)
        yc += row_h + (nsub - 1) * sub_h
    y_bot = yc - row_h + bh
    for xb in (ts, te):
        ax.plot([xb, xb], [y_h2 + h, y_bot + 0.6], color=GRAY, lw=0.4, dashes=(2.2, 1.6), zorder=0.5)
    return y_bot, X(ts), X(te), nread


def draw_B(P, E, FE, order, y0, xA_ts, xA_te, yA_bot, cell=1.36, half=JX_HALF):
    """Panel B: both junctions at base-pair resolution; each finding caller's call (its color, LTR end block full
    color, body tint) and the TSDs it reports at placeable positions (white boxes)."""
    r = E["r"]
    ts, te, d = r["ts"], r["te"], r["d"]
    nb = half[0] + half[1]
    wwin = nb * cell
    xl = min(max(xA_ts - wwin / 2, 22.0), P.W / 2 - wwin - 8)
    xr = max(min(xA_te - wwin / 2, P.W - 2 - wwin), P.W / 2 + 8)
    wins = [(xl, ts - half[0], "left"), (xr, te - half[1], "right")]
    found = [t for t in order if E["calls"][t] and E["calls"][t]["status"] == "found"]
    row_h = 2.1
    y_calls = y0 + 3.0
    y_h2 = y_calls + (len(found) - 1) * row_h + 4.2
    y_h1 = y_h2 + 2.8
    mono = dict(family="Liberation Mono", fontsize=6, ha="center", va="center")
    for x0, g0, side in wins:
        ax = P.ax(x0, y0, wwin, y_h1 + 2 - y0, xlim=(0, nb), ylim=(y_h1 + 2, y0))
        s2 = E["seq2"](g0, g0 + nb)
        if side == "left":
            h1 = E["seq1"](E["h_tsd"] - (half[0] - d), E["h_tsd"] + d) + "-" * (half[1])
            tsd_cols = range(half[0] - d, half[0])
            el_cols = range(half[0], nb)
            h1_tsd_cols = tsd_cols
        else:
            h1 = "-" * (half[1] + d) + E["seq1"](E["h_tsd"] + d, E["h_tsd"] + d + (half[0] - d))
            tsd_cols = range(half[1], half[1] + d)
            el_cols = range(0, half[1])
            h1_tsd_cols = range(0)
        for j in range(nb):
            inel = j in el_cols
            if inel:
                ax.add_patch(Rectangle((j, y_h2 - 1.0), 1, 2.0, fc=EL_LTR, lw=0))
            ax.text(j + 0.5, y_h2, s2[j], color="white" if inel else INK, **mono)
            b1 = h1[j]
            col = MUTED if (b1 == "-" or b1 == s2[j]) else "black"
            ax.text(j + 0.5, y_h1, b1, color=col, fontweight="bold" if col == "black" else "normal", **mono)
        for cols, yy in ((tsd_cols, y_h2), (h1_tsd_cols, y_h1)):
            cols = list(cols)
            if cols:
                ax.add_patch(Rectangle((cols[0] + 0.02, yy - 1.0), len(cols) - 0.04, 2.0, fc="none", ec=INK,
                                       lw=0.6))
        if list(tsd_cols):
            ax.text(list(tsd_cols)[0] + d / 2, y_h2 - 1.25, "TSD", ha="center", va="bottom", fontsize=6)
        xe = half[0] if side == "left" else half[1]
        ax.plot([xe, xe], [y_calls - 1.0, y_h1 + 1.2], color=GRAY, lw=0.4, dashes=(2.2, 1.6), zorder=0.5)
        for kk, t in enumerate(found):
            f = next(f for f in FE[t] if f["matched"])
            yy = y_calls + kk * row_h
            if side == "left":
                a, b = f["s"] - g0, nb
            else:
                a, b = 0, f["e"] - g0
            a, b = max(a, 0), min(b, nb)
            lo_ltr = f["lltr"] if side == "left" else f["rltr"]
            if lo_ltr is None:
                ax.add_patch(Rectangle((a, yy - 0.5), b - a, 1.0, fc=tint(COL[t], 0.75), lw=0))
            else:                                       # LTR part full color (the whole window: LTRs are > 200 bp)
                la, lb = (a, min(b, lo_ltr - g0)) if side == "left" else (max(a, lo_ltr - g0), b)
                ax.add_patch(Rectangle((a, yy - 0.5), b - a, 1.0, fc=tint(COL[t], 0.38), lw=0))
                if lb > la:
                    ax.add_patch(Rectangle((la, yy - 0.5), lb - la, 1.0, fc=COL[t], lw=0))
            for ta, tb in f["tsd_bp"]:
                ca, cb = max(ta - g0, 0), min(tb - g0, nb)
                if cb > ca:
                    ax.add_patch(Rectangle((ca + 0.05, yy - 0.42), cb - ca - 0.1, 0.84, fc="white", ec=INK, lw=0.4,
                                           zorder=3))
        xa = xA_ts if side == "left" else xA_te
        P.line([(xa, yA_bot + 0.6), (x0, y0 + 0.6)], color=LGRAY, lw=0.4)
        P.line([(xa, yA_bot + 0.6), (x0 + wwin, y0 + 0.6)], color=LGRAY, lw=0.4)
    for kk, t in enumerate(found):
        tool_lab(P, xl - 1.5, y_calls + kk * row_h, t)
        f = next(f for f in FE[t] if f["matched"])
        num("B", "TSDs drawn at base pairs (0-based)", f["tsd_bp"] or "none placeable",
            "reported coordinates, or a reported sequence equal to the adjacent flank", DS[0], t)
    lab_left(P, xl - 1.5, y_h2, "haplotype 2", fontsize=6.5)
    lab_left(P, xl - 1.5, y_h1, "haplotype 1", fontsize=6.5)
    inner = (te - half[1]) - (ts + half[1])
    P.text((xl + wwin + xr) / 2, y_h2, f"… {inner:,} bp …", ha="center", va="center", fontsize=6)
    P.text((xl + wwin + xr) / 2, y_h1, "empty site", ha="center", va="center", fontsize=6, color=MUTED)
    num("B", "TSD sequence and length", f"{r['tsd']} ({d} bp)", "report/tables/pool.tsv dup_seq/dup_len", DS[0])
    num("B", "element bases between the two windows (bp)", inner, "pool.tsv true ends", DS[0])
    return y_h1 + 1.5


def draw_breakdown(P, S, B, OWN, dim, pnl, y0, xL, xR, geo, axis_label):
    """Panels C (nesting) and D (protein content), fig_benchmark H's design per class: per caller (shared
    order) one 100% bar per haplotype (2 above 1): genotyped correctly in the caller's color, split into loci other
    callers also nominated (solid) and loci nominated by this caller alone (white hatch over its color, true width;
    every caller), the rest light gray; correct share and 95% CI (block bootstrap,
    ltrkin_breakdown.tsv) at the bar's end, zero counts as "0%" (a zero count's bootstrap interval has no width);
    tool failed in dark gray."""
    col, classes = DIMS[dim]
    ncls = len(classes)
    cw = (xR - xL - (ncls - 1) * geo["cgap"]) / ncls
    bw = cw - geo["tw"]
    bp, bh, pg = geo["bar_pitch"], geo["bar_h"], geo["pair_gap"]
    y_hdr = y0 + 1.0
    y_rows = y0 + geo["hdr"]
    pair = 2 * bp + pg
    y_end = y_rows + len(ORDER) * pair - pg - (bp - bh)
    for ci, c in enumerate(classes):
        x0 = xL + ci * (cw + geo["cgap"])
        ns = [int(B[(ds, "ltr_peel", dim, c)]["n_loci"]) for ds in DS]
        P.text(x0, y_hdr, CLSN[c], ha="left", va="center", fontsize=6.5, fontweight="bold")
        P.text(x0, y_hdr + 2.6, f"n = {ns[0]:,} | {ns[1]:,}", ha="left", va="center", fontsize=6, color=MUTED)
        P.text(x0 + cw, y_hdr + 2.6, "correct (95% CI)", ha="right", va="center", fontsize=6, color=MUTED)
        for ds, n_ in zip(DS, ns):
            num(pnl, f"loci in class {CLSN[c]}", n_, "experiments/figpub/ltrkin_breakdown.tsv n_loci", ds)
        ax = P.ax(x0, y_rows - 0.5, bw, len(ORDER) * pair + 1, xlim=(0, 1),
                  ylim=(y_rows - 0.5 + len(ORDER) * pair + 1, y_rows - 0.5))
        for k, t in enumerate(ORDER):
            for j, ds in enumerate(DS):
                yb = y_rows + k * pair + j * bp
                b = B[(ds, t, dim, c)]
                yc = yb + bh / 2
                if b["status"] != "ok":
                    reason = missing_rect(ax, 0, yb, 1, bh, 0.5, yc, b["status"])
                    num(pnl, f"{CLSN[c]}: status", f"{reason} ({b['status']})",
                        "experiments/figpub/ltrkin_breakdown.tsv status", ds, t)
                    continue
                n_ = int(b["n_loci"])
                kc = int(b["n_correct"])
                if n_ < MIN_CI:                        # counts only
                    P.text(x0, yc, f"{kc} of {n_} correct", ha="left", va="center", fontsize=6)
                    num(pnl, f"{CLSN[c]}: genotyped correctly (loci of n; counts only, n < {MIN_CI})",
                        f"{kc}/{n_}", "experiments/figpub/ltrkin_breakdown.tsv n_correct", ds, t)
                    continue
                w = kc / n_
                # hatching rule, every caller: its correct loci that no other caller nominated, at true width
                nown, kown = own_loci(OWN, ds, t, col, c)
                num(pnl, f"{CLSN[c]}: nominated by this caller alone (loci) | of them genotyped correctly (hatched) | "
                    f"genotyped correctly without them (solid, of class n)",
                    f"{nown} | {kown} | {kc - kown}/{n_} = {pct((kc - kown) / n_)}%",
                    "experiments/figpub/ltrkin_locus_class.tsv nominators (= this caller only), out_<caller>; "
                    "ltrkin_breakdown.tsv n_correct, n_loci", ds, t)
                wr = (kc - kown) / n_
                if wr > 0:
                    ax.add_patch(Rectangle((0, yb), wr, bh, fc=COL[t], lw=0))
                if kown:
                    ax.add_patch(Rectangle((wr, yb), w - wr, bh, fc=COL[t], hatch=HATCH, hatchcolor="white",
                                           lw=0))
                    HATCH_MM[(pnl, ds, t, c)] = (kown, n_, kown / n_ * bw)
                    num(pnl, f"{CLSN[c]}: hatch width at true scale (mm; bar {bw:.1f} mm; too thin to see if "
                        f"< {THIN_MM} mm)", f"{kown / n_ * bw:.2f}", "drawn: hatched loci / class n x bar width",
                        ds, t)
                ax.add_patch(Rectangle((w, yb), 1 - w, bh, fc=RESTC, lw=0))
                if abs(w - float(b["share_correct"])) > 5e-6:
                    die(f"{ds} {t} {c}: share {w} vs table {b['share_correct']}")
                # a zero count's block-bootstrap interval has no width: "0%" alone
                lab = "0%" if kc == 0 else f"{pct(w)}% ({pct(b['correct_lo'])}–{pct(b['correct_hi'])})"
                P.text(x0 + bw + 0.8, yc + 0.05, lab, ha="left", va="center", fontsize=6,
                       fontweight="bold" if t == "ltr_peel" else "normal")
                num(pnl, f"{CLSN[c]}: genotyped correctly, % (95% CI)", lab,
                    "experiments/figpub/ltrkin_breakdown.tsv share_correct, correct_lo, correct_hi", ds, t)
                num(pnl, f"{CLSN[c]}: loci correct | wrong genotype | wrong extent | missed (not drawn apart)",
                    " | ".join(b[f"n_{kk}"] for kk in CATS),
                    "experiments/figpub/ltrkin_breakdown.tsv n_correct, n_wrong, n_partial, n_missed", ds, t)
        # share axis under each column
        yax = y_end + 0.5
        P.can.plot([x0, x0 + bw], [yax, yax], color=INK, lw=0.5)
        for v in (0, 50, 100):
            xx = x0 + bw * v / 100
            P.can.plot([xx, xx], [yax, yax + 0.7], color=INK, lw=0.5)
            P.can.text(xx, yax + 1.0, f"{v}", ha="center", va="top", fontsize=6)
        if ci == 0:
            for k, t in enumerate(ORDER):
                yp = y_rows + k * pair
                tool_lab(P, xL - 2.6, yp + (bp + bh) / 2, t)
                for j, hn in enumerate(("2", "1")):
                    P.text(xL - 0.6, yp + j * bp + bh / 2, hn, ha="right", va="center", fontsize=6, color=MUTED)
            P.text(xL - 0.6, y_hdr + 2.6, "hap.", ha="right", va="center", fontsize=6, color=MUTED)
    yb_ = y_end + 0.5 + 3.4
    if axis_label:
        P.text(xL + (xR - xL) / 2, yb_ + 0.6, "loci genotyped correctly (%)",
               ha="center", va="top", fontsize=6.5)
        yb_ += 3.6
    return yb_


# ---------------------------------------------------------------------------------------------------- legend
def bshare(B, ds, t, dim, c, k="correct"):
    b = B[(ds, t, dim, c)]
    return int(b[f"n_{k}"]) / int(b["n_loci"])


def legend_text(T, S, B, OWN, E, FE, order, nread):
    """Legend page (markup for poafig.write_figure); every number is read from the tables or the per-locus data at
    run time, recorded in the numbers table, and every worded claim is checked first (check_claims)."""
    ds2, ds1 = DS
    Y = T["Y"]
    yv = lambda ds, t: Y[(ds, t, "e2e")]                                     # noqa: E731
    run2 = [t for t in order if t not in S[ds2]["ok"]]
    r = E["r"]
    found_ex = [t for t in order if E["calls"][t] and E["calls"][t]["status"] == "found"]
    miss_ex = [t for t in order if E["calls"][t] and E["calls"][t]["status"] != "found"]
    for t in miss_ex:
        num("legend", "share of haplotype 2 loci missed (no overlapping call)", pct(yv(ds2, t)["share_missed"]),
            f"report/tables/ltrkin_yield.tsv ({ds2}, {t}, e2e) share_missed", ds2, t)
        if FE[t]:
            die(f"{t} misses the example but has calls in the window; reword 'made no call here'")
    miss_lo = min(float(yv(ds2, t)["share_missed"]) for t in miss_ex) if miss_ex else None
    worst = max(found_ex, key=lambda t: max(abs(E["calls"][t]["err"]["left"]), abs(E["calls"][t]["err"]["right"])))
    we = E["calls"][worst]["err"]
    all_ok = all(E["calls"][t]["cat"] == "correct" for t in found_ex)
    fin = {ds: [t for t in order if t in S[ds]["ok"]] for ds in DS}
    others = {ds: [t for t in fin[ds] if t != "ltr_peel"] for ds in DS}
    host = {ds: pct(bshare(B, ds, "ltr_peel", "nesting", "host")) for ds in DS}
    host_next = {ds: max(others[ds], key=lambda t: bshare(B, ds, t, "nesting", "host")) for ds in DS}
    lq_rest = {}
    for dim, c in (("nesting", "host"), ("protein", "proteinless")):
        col = DIMS[dim][0]
        for ds in DS:
            nlq, klq = own_loci(OWN, ds, "ltr_peel", col, c)
            n_ = int(B[(ds, "ltr_peel", dim, c)]["n_loci"])
            lq_rest[(ds, c)] = (int(B[(ds, "ltr_peel", dim, c)]["n_correct"]) - klq, n_, nlq, klq)
    pl = {ds: pct(bshare(B, ds, "ltr_peel", "protein", "proteinless")) for ds in DS}
    pl_next = {ds: max(others[ds], key=lambda t: bshare(B, ds, t, "protein", "proteinless")) for ds in DS}
    for ds in DS:
        num("legend", "host: LTR_peel genotyped correctly (%) | best other caller (%)",
            f"{host[ds]} | {NAME[host_next[ds]]} {pct(bshare(B, ds, host_next[ds], 'nesting', 'host'))}",
            "experiments/figpub/ltrkin_breakdown.tsv", ds)
        num("legend", "protein-less: LTR_peel genotyped correctly (%) | best other caller (%)",
            f"{pl[ds]} | {NAME[pl_next[ds]]} {pct(bshare(B, ds, pl_next[ds], 'protein', 'proteinless'))}",
            "experiments/figpub/ltrkin_breakdown.tsv", ds)
    own_txt, thin_txt = own_sentences(OWN, S, order)
    check_claims(S, B, order, lq_rest, host_next, pl_next)
    gap_txt = ""
    if E["gaps1"]:
        q, sz = E["gaps1"][0]
        num("A", "indel in the drawn flanks: present only on haplotype 1 (bp)", sz, "hap2_on_hap1.paf.gz cg CIGAR", ds2)
        gap_txt = f" (the gap in the ribbons is {sz} bp present only on haplotype 1)"
    dnf_st = {Y[(ds2, t, "e2e")]["status"] for t in run2 if (ds2, t, "e2e") in Y}
    lim = " within 36 hours" if any("timeout" in x for x in dnf_st) else ""
    for t in run2:
        num("legend", "haplotype 2 run status", f"{Y[(ds2, t, 'e2e')]['status']}; cap 36 h",
            "report/tables/ltrkin_yield.tsv status (e2e); cap: ltrkin_geno/memo.md roll-in line (36.0 h)", ds2, t)
    run_txt = (f" {_and_list(NAME[t] for t in run2)} did not finish on haplotype 2{lim} (tool failed)."
               if run2 else "")
    head = ("**Fig. S4. Most callers lose the real insertion polymorphisms that carry another element inside them "
            "or no proteins of their own, so long-read genotyping misses these insertions.** One *Cardamine amara* "
            "plant, each haplotype in turn the reference. Only insertions that some caller's calls pointed to are "
            f"scored, so recall is relative.{run_txt} Fig. 1H shows the totals per caller and these class bars without "
            "intervals or the partly coding class.")
    pA = (f"**A**, One insertion polymorphism, a {(r['te'] - r['ts']) / 1000:.1f}-kb LTR retrotransposon on "
          f"haplotype 2 (gray) with its empty site on haplotype 1{gap_txt}. {_num_words(nread['carry']).capitalize()} "
          f"spanning reads carry the element and {nread['skip']} lack it, so the plant is heterozygous. Each caller's "
          f"calls are drawn with the features it reports. The {_num_words(len(found_ex))} callers that found the "
          f"element {'all ' if all_ok else ''}led to the correct genotype. {_and_list(NAME[t] for t in miss_ex)}, which "
          f"miss {'about half' if miss_lo and miss_lo >= 0.45 else 'many'} of all insertions, made no call here.")
    pB = (f"**B**, Both ends at base-pair resolution. A {r['d']}-bp target-site duplication (TSD) flanks the element "
          f"and occurs once at the empty site. White boxes on the call rows are the TSDs callers report at these "
          f"positions. {NAME[worst]}'s ends lie {abs(we['left'])} and {abs(we['right'])} bp "
          f"{'outside' if we['left'] < 0 and we['right'] > 0 else 'off'} the true ends, yet its call still led to "
          f"the correct genotype.")
    h2, h1 = lq_rest[(ds2, "host")], lq_rest[(ds1, "host")]
    hatch_txt = " ".join(x for x in (own_txt, thin_txt) if x)
    pC = (f"**C**, Insertions genotyped correctly by nesting, haplotype 2 above haplotype 1, with 95% confidence "
          f"intervals. {hatch_txt} Hosts contain an element that at least two callers agree on. Every other caller "
          f"genotypes hosts far less often than non-nested insertions. LTR_peel genotypes {host[ds2]}% and {host[ds1]}% of hosts (next best "
          f"{NAME[host_next[ds2]]} {pct(bshare(B, ds2, host_next[ds2], 'nesting', 'host'))}% and "
          f"{NAME[host_next[ds1]]} {pct(bshare(B, ds1, host_next[ds1], 'nesting', 'host'))}%), still {h2[0]} and "
          f"{h1[0]} of {h2[1]} without its hatched loci.")
    p2, p1 = lq_rest[(ds2, "proteinless")], lq_rest[(ds1, "proteinless")]
    pD = (f"**D**, The same by the LTR retrotransposon protein domains DANTE finds in each insertion. Coding "
          f"insertions carry at least two kinds from one superfamily, protein-less ones none, partly coding ones the "
          f"rest. Every caller genotypes protein-less insertions less often than coding ones, and DANTE_LTR, which "
          f"builds elements from these domains, genotypes none. LTR_peel genotypes {pl[ds2]}% and {pl[ds1]}%, but {p2[3]} of "
          f"{p2[1]} and {p1[3]} of {p1[1]} are its hatched loci. Without them it reaches {pct(p2[0] / p2[1])}% and "
          f"{pct(p1[0] / p1[1])}%, against {pct(bshare(B, ds2, pl_next[ds2], 'protein', 'proteinless'))}% and "
          f"{pct(bshare(B, ds1, pl_next[ds1], 'protein', 'proteinless'))}% for the best other caller "
          f"({NAME[pl_next[ds2]]}, {NAME[pl_next[ds1]]}). The genotyper was developed with LTR_peel annotations.")
    return "\n\n".join([head, pA, pB, pC, pD])


def own_sentences(OWN, S, order):
    """Legend sentences for the hatching rule, worded from the per-locus table and the drawn hatch widths: each
    caller's correctly genotyped loci found by it alone (totals per haplotype), and the hatches too thin to see
    (< THIN_MM). Exits if the wording no longer fits (another caller with as many own correct loci as LTR_peel; the
    thin hatches not exactly the hatches of one or two loci in the classes named)."""
    ds2, ds1 = DS
    k, bad = {}, []
    for ds in DS:
        for t in order:
            if t not in S[ds]["ok"]:
                num("legend", "nominated by this caller alone, all loci | of them genotyped correctly (hatched)",
                    "tool failed", "report/tables/ltrkin_yield.tsv status", ds, t)
                continue
            nown, kown = own_loci(OWN, ds, t)
            k[(ds, t)] = kown
            num("legend", "nominated by this caller alone, all loci | of them genotyped correctly (hatched)",
                f"{nown} | {kown}", "experiments/figpub/ltrkin_locus_class.tsv nominators (= this caller only), "
                "out_<caller>", ds, t)
    oth = [t for t in order if t != "ltr_peel"]
    for ds in DS:
        if any(k.get((ds, t), 0) >= k[(ds, "ltr_peel")] for t in oth):
            bad.append(f"{ds}: another caller has as many correct own loci as LTR_peel")
    both, single, none = [], collections.OrderedDict(), []
    for t in oth:
        a, b = k.get((ds2, t), 0), k.get((ds1, t), 0)
        if a and b:
            both.append(f"{NAME[t]} {a} and {b}")
        elif a or b:
            single.setdefault((("2", a) if a else ("1", b)), []).append(NAME[t])
        else:
            none.append(t)
    parts = both + [f"{_and_list(ts)} {_num_words(v)}{' each' if len(ts) > 1 else ''} on haplotype {h}"
                    for (h, v), ts in single.items()]
    own_txt = (f"Hatched, loci genotyped correctly that only this caller's calls pointed to (LTR_peel "
               f"{k[(ds2, 'ltr_peel')]} and {k[(ds1, 'ltr_peel')]}, " + ", ".join(parts)
               + (", the others none)." if none else ")."))
    # "Hatches of one or two non-nested or coding loci are too thin to see": the thin hatches must be exactly the
    # hatches of those classes with that many loci (every one of them thin, no other hatch thin)
    thin = {key for key, v in HATCH_MM.items() if v[2] < THIN_MM}
    if not thin:
        return own_txt, ""
    cls_thin = [c for dim in ("nesting", "protein") for c in DIMS[dim][1] if any(key[3] == c for key in thin)]
    cnt = sorted({HATCH_MM[key][0] for key in thin})
    named = {key for key, v in HATCH_MM.items() if key[3] in cls_thin and v[0] <= cnt[-1]}
    if named != thin or len(cnt) > 2 or cnt[0] != 1:
        bad.append(f"thin hatches {sorted(thin)} differ from those the wording names {sorted(named)} "
                   f"(classes {cls_thin}, loci {cnt})")
    if bad:
        die("a legend claim about the hatching no longer holds:\n  " + "\n  ".join(bad))
    for key in sorted(thin):
        num("legend", f"hatch too thin to see (< {THIN_MM} mm): {CLSN[key[3]]}, loci (drawn width mm)",
            f"{HATCH_MM[key][0]} ({HATCH_MM[key][2]:.2f})",
            "experiments/figpub/ltrkin_locus_class.tsv; width = loci / class n x bar width", key[1], key[2])
    nw = " or ".join(_num_words(x) for x in cnt)
    thin_txt = (f"Hatches of {nw} {_and_list((CLSN[c] for c in cls_thin), 'or')} loci are too thin to see.")
    return own_txt, thin_txt


def check_claims(S, B, order, lq_rest, host_next, pl_next):
    """Exit before anything is written if a rerun makes a claim worded in the legend false; the numbers behind each
    claim go to the numbers table."""
    bad = []
    for ds in DS:
        fin = [t for t in order if t in S[ds]["ok"]]
        oth = [t for t in fin if t != "ltr_peel"]
        for t in oth:                                  # C: host share < non-nested share - 0.15 ("far less often")
            h, u = bshare(B, ds, t, "nesting", "host"), bshare(B, ds, t, "nesting", "unnested")
            if not h < u - 0.15:
                bad.append(f"{ds} {t}: host {h:.3f} vs non-nested {u:.3f}, not 'far less often'")
        lq = int(B[(ds, "ltr_peel", "nesting", "host")]["n_correct"])
        nx = int(B[(ds, host_next[ds], "nesting", "host")]["n_correct"])
        if not (lq > nx and lq_rest[(ds, "host")][0] > nx):
            bad.append(f"{ds}: LTR_peel host {lq} (without its own loci {lq_rest[(ds, 'host')][0]}) vs next {nx}")
        for t in fin:                                  # D: every caller: protein-less < coding; DANTE_LTR none
            if not bshare(B, ds, t, "protein", "proteinless") < bshare(B, ds, t, "protein", "coding"):
                bad.append(f"{ds} {t}: protein-less not below coding")
        if "dante_ltr" in fin and int(B[(ds, "dante_ltr", "protein", "proteinless")]["n_correct"]) != 0:
            bad.append(f"{ds}: DANTE_LTR genotypes some protein-less insertions")
        if not bshare(B, ds, "ltr_peel", "protein", "proteinless") > \
                bshare(B, ds, pl_next[ds], "protein", "proteinless"):
            bad.append(f"{ds}: LTR_peel not ahead on protein-less insertions")
    if bad:
        die("a legend claim no longer holds; reword the legend before rebuilding:\n  " + "\n  ".join(bad))


def _num_words(k):
    w = ("zero one two three four five six seven eight nine ten eleven twelve thirteen fourteen fifteen sixteen "
         "seventeen eighteen nineteen twenty").split()
    return w[k] if 0 <= k < len(w) else str(k)


def _and_list(xs, conj="and"):
    xs = list(xs)
    if not xs:
        return "none"
    return xs[0] if len(xs) == 1 else ", ".join(xs[:-1]) + f" {conj} " + xs[-1]


# ---------------------------------------------------------------------------------------------------- source data
# One TSV per table in source_data/figS4_ltrkin/; coordinates 0-based half-open on the named haplotype; empty cell =
# not reported (None). The figure and the numbers table are drawn from these tables only.
SCALARS = dict(example_locus=str, chrom=str, true_start=int, true_end=int, tsd_seq=str, tsd_len=int,
               ltr_identity=float, ltr5_len=int, ltr3_len=int, strand=str, classification=str, window_start=int,
               window_end=int, hap1_chrom=str, hap1_tsd_start=int, median_locus_len=None)
FEATS = ("tsd_drawn", "tsd_bp", "domain")            # A_call_features.tsv feature -> FE key tsd, tsd_bp, dom


def write_tsv(name, header, rows):
    """TSV with a header row; None -> empty cell; floats at full precision (repr)."""
    path = SRC_DIR / name
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        for r in rows:
            w.writerow(["" if x is None else x for x in r])
    return path


def read_tsv(name):
    return rd(SRC_DIR / name)


def opt_int(s):
    return None if s == "" else int(s)


def num_value(s):
    """int or float as written (repr), so a value prints the same after reading back"""
    return int(s) if re.fullmatch(r"-?\d+", s) else float(s)


def write_source(T, S, B, LC, ex, E, FE, SRC):
    """Every value the figure, legend and numbers table use, from the full dataset, to SRC_DIR."""
    SRC_DIR.mkdir(parents=True, exist_ok=True)
    for p in SRC_DIR.glob("*.tsv*"):
        p.unlink()
    ds2 = DS[0]
    r = E["r"]
    sc = dict(example_locus=r["locus"], chrom=r["chrom"], true_start=r["ts"], true_end=r["te"], tsd_seq=r["tsd"],
              tsd_len=r["d"], ltr_identity=r["ltr_id"], ltr5_len=E["l5"], ltr3_len=E["l3"], strand=E["strand"],
              classification=E["cls"], window_start=E["W0"], window_end=E["W1"], hap1_chrom=E["ochrom"],
              hap1_tsd_start=E["h_tsd"], median_locus_len=ex["median"])
    write_tsv("scalars.tsv", ["name", "value"], [(k, repr(sc[k]) if isinstance(sc[k], float) else sc[k])
                                                  for k in SCALARS])
    write_tsv("A_funnel.tsv", ["step", "loci"], ex["funnel"])
    write_tsv("A_eligible.tsv", ["locus", "length_bp", "reads_carry", "reads_skip"],
              [(x["locus"], x["len"], x["carry"], x["skip"]) for x in ex["eligible"]])
    write_tsv("A_aligned_blocks.tsv", ["hap2_start", "hap2_end", "hap1_start", "hap1_end"], E["blocks"])
    write_tsv("A_hap1_only_indels.tsv", ["hap2_pos", "length_bp"], E["gaps1"])
    write_tsv("A_reads.tsv", ["read", "kind", "start", "end"],
              [(k, kind, s, e) for k, (kind, parts) in enumerate(E["reads"].values()) for s, e in parts])
    write_tsv("callers.tsv", ["dataset", "caller", "run_status", "finished"],
              [(ds, t, S[ds]["res"][t]["status"], int(t in S[ds]["ok"])) for ds in DS for t in PRIMARY])
    rows = []
    for t in ORDER:
        info = E["calls"][t]
        if info is not None:
            err = info["err"]
            rows.append((t, info["status"], info["cat"], info["geno"][0], err and err["left"], err and err["right"]))
    write_tsv("A_caller_outcomes.tsv", ["caller", "status", "outcome", "genotype", "end_err_left", "end_err_right"],
              rows)
    calls, feats = [], []
    for t in ORDER:
        for k, f in enumerate(FE[t] or []):
            calls.append((t, k, f["s"], f["e"], int(f["matched"]), f["strand"], f["lltr"], f["rltr"]))
            for name, key in zip(FEATS, ("tsd", "tsd_bp", "dom")):
                for x in f[key]:
                    feats.append((t, k, name, x[0], x[1], x[2] if len(x) > 2 else None))
    write_tsv("A_calls.tsv", ["caller", "call", "start", "end", "at_locus", "strand", "lltr_end", "rltr_start"], calls)
    write_tsv("A_call_features.tsv", ["caller", "call", "feature", "start", "end", "name"], feats)
    write_tsv("A_feature_sources.tsv", ["caller", "feature", "source"],
              [(t, k, v) for t in ORDER if t in SRC for k, v in SRC[t].items()])
    h, d, nb = JX_HALF, r["d"], sum(JX_HALF)
    seqs = [("2", r["chrom"], a, a + nb) for a in (r["ts"] - h[0], r["te"] - h[1])]
    seqs += [("1", E["ochrom"], E["h_tsd"] - (h[0] - d), E["h_tsd"] + d),
             ("1", E["ochrom"], E["h_tsd"] + d, E["h_tsd"] + d + (h[0] - d))]
    write_tsv("B_junction_seqs.tsv", ["haplotype", "chrom", "start", "end", "seq"],
              [(hp, ch, a, b, (E["seq2"] if hp == "2" else E["seq1"])(a, b)) for hp, ch, a, b in seqs])
    cols = ["status", "n_loci", "n_correct", "n_wrong", "n_partial", "n_missed", "share_correct", "correct_lo",
            "correct_hi"]
    write_tsv("CD_breakdown.tsv", ["dataset", "caller", "dimension", "class"] + cols,
              [(ds, t, dim, c) + tuple(B[(ds, t, dim, c)][k] for k in cols)
               for ds in DS for t in ORDER for dim in DIMS for c in DIMS[dim][1]])
    write_tsv("CD_own_loci.tsv", ["dataset", "caller", "class_column", "class", "n_own", "n_own_correct"],
              [(ds, t, col or "", c or "") + own_loci_raw(LC, ds, t, col, c) for ds in DS for t in ORDER
               for col, c in [(None, None)] + [(DIMS[dim][0], c) for dim in DIMS for c in DIMS[dim][1]]])
    write_tsv("legend_yield.tsv", ["dataset", "caller", "run", "status", "share_missed"],
              [(ds2, t, "e2e", T["Y"][(ds2, t, "e2e")]["status"], T["Y"][(ds2, t, "e2e")]["share_missed"])
               for t in ORDER if (ds2, t, "e2e") in T["Y"]])


def read_source():
    """The objects the drawing code uses, rebuilt from SRC_DIR only."""
    if not SRC_DIR.is_dir():
        die(f"source data not found: {SRC_DIR}")
    sc = {x["name"]: x["value"] for x in read_tsv("scalars.tsv")}
    sc = {k: (num_value(sc[k]) if ty is None else ty(sc[k])) for k, ty in SCALARS.items()}
    funnel = [(x["step"], int(x["loci"])) for x in read_tsv("A_funnel.tsv")]
    elig = [dict(locus=x["locus"], len=int(x["length_bp"]), carry=int(x["reads_carry"]), skip=int(x["reads_skip"]))
            for x in read_tsv("A_eligible.tsv")]
    ex = dict(funnel=funnel, eligible=elig, median=sc["median_locus_len"])
    S = {ds: dict(ok=[], res={}) for ds in DS}
    for x in read_tsv("callers.tsv"):
        S[x["dataset"]]["res"][x["caller"]] = dict(status=x["run_status"])
        if x["finished"] == "1":
            S[x["dataset"]]["ok"].append(x["caller"])
    reads = {}
    for x in read_tsv("A_reads.tsv"):
        reads.setdefault(int(x["read"]), (x["kind"], []))[1].append((int(x["start"]), int(x["end"])))
    calls = {t: None for t in ORDER}
    for x in read_tsv("A_caller_outcomes.tsv"):
        err = None if x["end_err_left"] == "" else dict(left=int(x["end_err_left"]), right=int(x["end_err_right"]))
        calls[x["caller"]] = dict(status=x["status"], cat=x["outcome"], geno=(x["genotype"],), err=err)
    FE = {t: (None if t not in S[DS[0]]["ok"] else []) for t in ORDER}
    for x in read_tsv("A_calls.tsv"):
        FE[x["caller"]].append(dict(s=int(x["start"]), e=int(x["end"]), matched=x["at_locus"] == "1",
                                    strand=x["strand"] or None, lltr=opt_int(x["lltr_end"]),
                                    rltr=opt_int(x["rltr_start"]), tsd=[], tsd_bp=[], dom=[]))
    for x in read_tsv("A_call_features.tsv"):
        f = FE[x["caller"]][int(x["call"])]
        iv = (int(x["start"]), int(x["end"]))
        f[dict(zip(FEATS, ("tsd", "tsd_bp", "dom")))[x["feature"]]].append(
            iv + (x["name"],) if x["feature"] == "domain" else iv)
    FE["_src"] = {}
    for x in read_tsv("A_feature_sources.tsv"):
        FE["_src"].setdefault(x["caller"], {})[x["feature"]] = x["source"]
    seqs = {(x["haplotype"], int(x["start"]), int(x["end"])): x["seq"] for x in read_tsv("B_junction_seqs.tsv")}
    r = dict(locus=sc["example_locus"], chrom=sc["chrom"], ts=sc["true_start"], te=sc["true_end"], d=sc["tsd_len"],
             tsd=sc["tsd_seq"], ltr_id=sc["ltr_identity"])
    E = dict(r=r, W0=sc["window_start"], W1=sc["window_end"], h_tsd=sc["hap1_tsd_start"], ochrom=sc["hap1_chrom"],
             blocks=[tuple(int(x[k]) for k in ("hap2_start", "hap2_end", "hap1_start", "hap1_end"))
                     for x in read_tsv("A_aligned_blocks.tsv")],
             gaps1=[(int(x["hap2_pos"]), int(x["length_bp"])) for x in read_tsv("A_hap1_only_indels.tsv")],
             calls=calls, reads=reads, l5=sc["ltr5_len"], l3=sc["ltr3_len"], strand=sc["strand"],
             cls=sc["classification"], seq2=lambda a, b: seqs[("2", a, b)], seq1=lambda a, b: seqs[("1", a, b)])
    B = {(x["dataset"], x["caller"], x["dimension"], x["class"]): x for x in read_tsv("CD_breakdown.tsv")}
    OWN = {(x["dataset"], x["caller"], x["class_column"], x["class"]): (int(x["n_own"]), int(x["n_own_correct"]))
           for x in read_tsv("CD_own_loci.tsv")}
    T = dict(Y={(x["dataset"], x["caller"], x["run"]): x for x in read_tsv("legend_yield.tsv")})
    return T, S, B, OWN, ex, E, FE


# ---------------------------------------------------------------------------------------------------- main
GEO = dict(cgap=3.0, tw=14.6, bar_pitch=2.1, bar_h=1.6, pair_gap=0.45, hdr=5.4)


def build(H, T, S, B, OWN, E, FE, order, base, v):
    """Draw the figure page at height H (mm); returns (Page, content bottom in mm, read counts)."""
    NUMS[:] = base
    HATCH_MM.clear()
    W = 180.0
    P = Page(W, H)
    xL, xR = 21.0, 178.0
    key_A(P, 2.6, xL)
    yA = 8.4
    P.letter(0.5, yA - 1.0, "A")
    yA_bot, xts, xte, nread = draw_A(P, E, FE, S, order, yA, xL, xR, v)
    yB = yA_bot + 3.6
    P.letter(0.5, yB + 0.6, "B")
    yB_bot = draw_B(P, E, FE, order, yB + 0.6, xts, xte, yA_bot)
    yK = yB_bot + 3.4
    key_C(P, yK, xL)
    yC = yK + 2.6
    P.letter(0.5, yC - 1.0, "C")
    yC_bot = draw_breakdown(P, S, B, OWN, "nesting", "C", yC, xL, xR, GEO, False)
    yD = yC_bot + 1.6
    P.letter(0.5, yD, "D")
    yD_bot = draw_breakdown(P, S, B, OWN, "protein", "D", yD, xL, xR, GEO, True)
    log(f"layout: A {yA:.1f}-{yA_bot:.1f}, B -{yB_bot:.1f}, C {yC:.1f}-{yC_bot:.1f}, D -{yD_bot:.1f} of {H} mm", v)
    return P, yD_bot, nread


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("-v", "--verbose", action="store_true", help="per-step progress and the example funnel on stderr")
    ap.add_argument("--draft", metavar="PNG", help="write only a figure PNG at this path (inside this directory) "
                    "and stop; no PDF, PNG or numbers table")
    ap.add_argument("--from-source-data", action="store_true", help="redraw from source_data/figS4_ltrkin/ in "
                    "this directory alone (no BENCH_ROOT, no raw-data checks)")
    a = ap.parse_args()
    v = a.verbose
    if a.draft and not os.path.abspath(a.draft).startswith(str(HERE) + os.sep):
        die(f"--draft must write inside {HERE}")
    log(f"{STEM}: start", True)
    if not a.from_source_data:
        if BENCH is None:
            die("BENCH_ROOT is not set (the full dataset is needed unless --from-source-data)")
        import_raw_modules()
        T = load_tables()
        running = running_set(T)
        log(f"runs still in progress (from ltrkin_yield.tsv): {sorted(running) or 'none'}", v)
        S = score_all(T, running, v)
        caller_order(T, S)
        B, LC = load_breakdown(S, v)
        ex = select_example(S, v)
        E = example_data(S, ex, v)
        FE, SRC = reported_features(E)
        write_source(T, S, B, LC, ex, E, FE, SRC)
        log(f"{STEM}: wrote the source data to {SRC_DIR}", True)
    T, S, B, OWN, ex, E, FE = read_source()
    order = list(ORDER)
    for k, val in ex["funnel"]:
        num("A", f"example funnel: {k}", val, "selection rule in this script (memo.md)", DS[0])
    num("A", "median length of all haplotype 2 pooled loci (bp)", ex["median"], "report/tables/pool.tsv match span",
        DS[0])
    fmt_l = lambda rr: "; ".join(f"{x['locus']} ({x['len']} bp, carry/skip {x['carry']}/{x['skip']})" for x in rr)  # noqa: E731
    num("A", "eligible loci under the rule, nearest the median length first", fmt_l(ex["eligible"]),
        "selection rule (memo.md)", DS[0])
    num("A", "example locus", f"{E['r']['locus']} {E['r']['chrom']}:{E['r']['ts']}-{E['r']['te']} (0-based)",
        "selection rule (memo.md)", DS[0])
    num("A", "LTR identity, left vs right LTR (matches / longer LTR)", round(E["r"]["ltr_id"], 4),
        "parasail global alignment of LTR_peel's two LTRs (hap2 assembly)", DS[0], "ltr_peel")
    num("A", "haplotype 1 empty site: TSD start (0-based)", f"{E['ochrom']}:{E['h_tsd']}",
        "Camara_hap1.fa, junction check (30 bp each side, <= 3 mismatches, TSD exact, one hit)", DS[1])
    if E["r"]["locus"] != "camara_hap2_L00918":
        die(f"example locus changed to {E['r']['locus']} (source figure: camara_hap2_L00918); check before drawing")
    base = list(NUMS)
    P, bot, nread = build(240.0, T, S, B, OWN, E, FE, order, base, v)
    H = math.ceil(bot + 1.0)
    plt.close(P.fig)
    P, bot, nread = build(H, T, S, B, OWN, E, FE, order, base, v)
    legend = legend_text(T, S, B, OWN, E, FE, order, nread)
    if a.draft:
        P.fig.savefig(a.draft, dpi=300)
        log(f"{STEM}: draft {a.draft} ({H} mm tall)", True)
        print(legend)
        return
    out = poafig.write_figure(STEM, P.fig, legend, None, HERE, str(SRC_DIR if a.from_source_data else BENCH),
                              str(HERE / f"{STEM}.py"))
    with open(HERE / f"{STEM}_numbers.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["panel", "item", "dataset", "caller", "value", "source"], delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        for row in NUMS:
            w.writerow(row)
    log(f"{STEM}: wrote {out['pdf']}, {out['png'][0]} ({out['fig_mm']}; legend {out['short_mm']} mm, "
        f"{out['short_words']} words {out['short_par_words']}), {STEM}_numbers.tsv ({len(NUMS)} rows)", True)


if __name__ == "__main__":
    main()
