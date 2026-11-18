#!/usr/bin/env python3
"""fig_benchmark: the main benchmark figure.

Where LTR retrotransposon callers lose elements, and what recovering them costs. Panels (8):
  A  one typical nested locus of simulated genome S2 seed 1: nesting diagram
     (host in two pieces, its tenant), every tool's calls to scale (tool rows shared with B)
  B  one scorecard, S2 mean of 3 genomes: precision, recall, F1 for non-nested | tenants |
     hosts (nesting) and coding | protein-less (protein domains), one basis for both systems (sim_class_prf_joint)
  C  recall of CODING elements vs median LTR divergence (K2P) in the aging series, pooled line with a 95%
     family-resampled band, one dot at each pooled curve's recall-halved point, documented identity floors dashed
  D  the same form for the protein-less elements of the series (same axes as C; shared key between them)
  E  masking of real genomes as precision, recall, F1 for A. thaliana, rice, zebra finch
  F  share of curated LTR families recovered perfect or good (BM_RM2), 5 genomes
  G  A. thaliana: curated annotation along every protein-confirmed locus, rows by LTR divergence, beside the same
     spans moved 20-200 kb (chance), mean profiles (rendering as figS3_library_age)
  H  ltrkin bars: per caller, haplotype 2 above haplotype 1, all loci (share and 95% CI at right) and one 0-100% bar
     per class (non-nested, tenants, hosts | coding, protein-less); genotyped correctly in the caller's color, loci
     nominated by this caller alone hatched white inside it, the rest light gray
Numbers: experiments/figpub/{sim_class_prf_joint_summary,sim_class_prf_summary,ltrkin_outcomes_all,
ltrkin_locus_class}.tsv, report/tables/{masking,bm_rm2,resources,library_age_bias}.tsv, experiments/nonauto/
{div_class,div_elements,div_d50}.tsv, experiments/a4_diagnostics/annotation_gap.tsv,
ltrkin_geno/report/tables/ltrkin_yield.tsv.
Drawings: sims/S2_1 (FASTA, BED; gen500000 for the example check), sims/truth/S2_1.truth.tsv, calls/S2_1/<tool>.tsv,
eval/adj/arab/<tool>.calls.tsv, mask/arab/curated/arab.fa.out. Imported read-only: scripts/score_elements.py,
adj_summary.py, proxies.py, div_score.py; poafig (house library).
Outputs (this directory): fig_benchmark.{pdf,png}, fig_benchmark_numbers.tsv; caches and checks in tmp/.
Usage: fig_benchmark.py [-v] [--redo-prep] [--honesty] [--draft NAME]
"""
import argparse
import csv
import glob
import inspect
import math
import os
import statistics
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
os.environ["PYTHONDONTWRITEBYTECODE"] = "1"
sys.dont_write_bytecode = True                     # nothing written into bench/scripts or ltrkin_geno

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

B = Path(os.environ["BENCH_ROOT"])
HERE = Path(__file__).resolve().parent
TMP = HERE / "tmp"
STEM = "fig_benchmark"
LK = B / "ltrkin_geno"                              # another project: read only
FP = B / "experiments/figpub"
NA = B / "experiments/nonauto"
TAB = B / "report/tables"
sys.path.insert(0, str(B / "scripts"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # poafig.py sits beside this script
import score_elements as se  # noqa: E402
import adj_summary as su  # noqa: E402
import proxies  # noqa: E402
import div_score  # noqa: E402  (d50 rule)
import poafig  # noqa: E402
import matplotlib as mpl  # noqa: E402
import matplotlib.font_manager as fm  # noqa: E402

# fonts: Liberation registered in-process, poafig._advance_pt patched for matplotlib >= 3.11 (poafig.py itself is not modified)
for _f in sorted(glob.glob("/usr/share/fonts/truetype/liberation2/LiberationSans-*.ttf")
                 + glob.glob("/usr/share/fonts/truetype/liberation2/LiberationMono-*.ttf")):
    fm.fontManager.addfont(_f)
poafig.setup_mpl()
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, to_rgb  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402


def _advance_pt_mpl311(text, style):
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


from matplotlib import _text_helpers  # noqa: E402
if "kern_mode" not in inspect.signature(_text_helpers.layout).parameters:
    poafig._advance_pt = _advance_pt_mpl311

# ============================================================================ constants
TOOLS = ["ltrquest", "ltr_retriever", "ltrharvest", "ltr_finder", "look4ltrs", "hirepeat", "hite", "ltrpred",
         "dante_ltr"]
NAME = dict(ltrquest="LTRquest", ltr_retriever="LTR_retriever", ltrharvest="LTRharvest", ltr_finder="LTR_Finder",
            look4ltrs="Look4LTRs", hirepeat="HiRepEat", hite="HiTE", ltrpred="LTRpred", dante_ltr="DANTE_LTR")
COL = dict(ltrquest="#D55E00", ltr_retriever="#0072B2", ltrharvest="#56B4E9", ltr_finder="#009E73",
           look4ltrs="#CC79A7", hirepeat="#E69F00", hite="#000000", ltrpred="#A6761D", dante_ltr="#6A3D9A")
DS = ["arab", "rice", "finch", "maize", "gorilla"]
DSNAME = {"arab": r"$\it{A.\ thaliana}$", "rice": "rice", "finch": "zebra finch", "maize": "maize",
          "gorilla": "gorilla"}
S2 = ["S2_1", "S2_2", "S2_3"]
MM = poafig.MM
INK, MUTED = poafig.INK, poafig.MUTED
FS = 6.0                                        # every label >= 6 pt
LW_Q, LW_O = 1.3, 0.8                           # line weights C, D: LTRquest heavier
BAND_A = 0.22                                   # band alpha C, D
CMAP = mpl.colormaps["viridis"]                 # score map B, E, F (viridis, fixed 0-1), CMAP(v)
MISSC, HATCH, HATCHC = "#F2F2F2", "//////", "#B8B8B8"
NOCALLC = "#F2F2F2"                             # B: "–" cell, the caller makes no call of that class (light neutral)
HOSTC, HOSTC_LTR = "#5E5E5E", "#2B2B2B"         # A truth: host
TENC, TENC_LTR = "#B9B9B9", "#8A8A8A"           # A truth: tenant
NULLC = "#8C8C8C"
RESTC = "#DDDDDD"                                # H: not genotyped correctly (light gray remainder)
HATCH_LW_H, HATCH_STEP_H = 0.6, 1.4             # H: white hatch drawn as lines (pt width; mm between lines along x)
# the A example (selection rule and checks in the code below)
EX = dict(ds="S2_1", chrom="chr1", host="chr1:11858003-11868332", tenant="chr1:11861453-11866676",
          nonnested="chr1:11885081-11890067", win=(11_857_500, 11_890_600), old_gen="gen500000",
          segs=((11_857_500, 11_869_000), (11_884_600, 11_890_600)),     # A: broken axis, same bp scale
          old_host=(11_686_375, 11_691_476))
RUN = B / "runs" / "S2_1"                         # the example genome's tool outputs (sequence ids c0001.. = chr1..)
FAILC = "#5A5A5A"                                 # missing result, tool side: "tool failed"
DOMC = "#2E2E2E"                                  # A: protein domains reported by a tool
SAMTOOLS = "samtools"      # on PATH (env bench_tools)
# G (figS3_library_age C, A. thaliana only)
NB5, NBI, NB3 = 10, 30, 10
NBIN = NB5 + NBI + NB3
BINS = ["<0.01", "0.01-0.03", "0.03-0.1", ">=0.1"]
LOCIC = "#6A3D9A"
LTR5C, LTR3C, INTC = "#D55E00", "#EE9F6B", "#F6D9C7"
VERBOSE = False
NUM = []


def log(msg, force=False):
    if VERBOSE or force:
        print(f"[{STEM}] {msg}", file=sys.stderr, flush=True)


def num(panel, quantity, value, source):
    NUM.append(dict(panel=panel, quantity=quantity, value=value, source=source))
    return value


def tint(hexc, a):
    """hexc blended with white; a = share of the color"""
    r, g, b = to_rgb(hexc)
    return (1 - a + a * r, 1 - a + a * g, 1 - a + a * b)


def read_tsv(p):
    with open(p) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def pct(x):
    """share -> whole percent, rounded half up"""
    return int(100 * float(x) + 0.5)


# ============================================================================ data: A (example locus)
def fa_path(ds, gen="gen1000000"):
    """symlink a simulated FASTA into tmp/ and index it there (nothing written next to the source)"""
    TMP.mkdir(parents=True, exist_ok=True)
    link = TMP / (f"{ds}.fa" if gen == "gen1000000" else f"{ds}.{gen}.fa")
    if not link.exists():
        link.symlink_to(B / "sims" / ds / f"{gen}_final.fasta")
    if not Path(str(link) + ".fai").exists():
        subprocess.run(["nice", SAMTOOLS, "faidx", str(link)], check=True)
    return str(link)


def fetch(fa, chrom, s, e):
    """0-based half-open [s, e) -> samtools region (1-based closed) s+1..e"""
    out = subprocess.run([SAMTOOLS, "faidx", fa, f"{chrom}:{s + 1}-{e}"], capture_output=True, text=True,
                         check=True).stdout
    return "".join(out.split("\n")[1:]).upper()


def gff_rows(path, seqid, w0, w1):
    """GFF3 rows on seqid overlapping [w0, w1) (0-based); GFF3 is 1-based closed -> start - 1, end"""
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


def reported_features(calls, w0, w1):
    """per tool and call in the window: what the tool itself reports (strand, LTR bounds, TSDs, protein domains),
    from its own output. Returns {tool: [feat]} and {tool: {feature: source}}. Coordinates 0-based half-open."""
    seq = "c0001"                                  # data/S2_1/names.tsv: c0001 = chr1, same coordinates
    F, SRC = {}, {}
    span = lambda c: (c["start"], c["end"])        # noqa: E731
    flank = lambda c, n: [(c["start"] - n, c["start"]), (c["end"], c["end"] + n)]   # noqa: E731
    for t in TOOLS:
        F[t] = []
        for c in calls[t]:
            F[t].append(dict(s=c["start"], e=c["end"], matched=c["matched"], id=c["id"],
                             strand=c["strand"] if c["strand"] in ("+", "-") else None,
                             lltr=c["lltr"], rltr=c["rltr"], tsd=None, dom=[]))
    # LTRquest: strand + LTRs (calls tsv <- clean_ltr.tsv); TSD sequence, offsets 0,0 (clean_ltr.tsv tsd); domains
    # protein_match children of its LTR_retrotransposon in the GFF3
    g = gff_rows(RUN / "ltrquest/S2_1_LTRs_all_depth_LTR_cleaned.gff3", seq, w0, w1)
    ids = {(r["s"], r["e"]): r["att"]["ID"] for r in g if r["type"] == "LTR_retrotransposon"}
    for f in F["ltrquest"]:
        c = next(x for x in calls["ltrquest"] if x["id"] == f["id"])
        if c["tsd"] not in (".", ""):
            f["tsd"] = flank(c, len(c["tsd"]))
        pid = ids[span(c)]
        f["dom"] = [(r["s"], r["e"], r["att"]["Name"]) for r in g
                    if r["type"] == "protein_match" and r["att"].get("Parent") == pid]
    SRC["ltrquest"] = dict(strand="calls tsv strand (clean_ltr.tsv)", LTR="calls tsv lltr_end/rltr_start",
                           TSD="clean_ltr.tsv tsd sequence, tsd_offset 0,0 (flanks drawn adjacent)",
                           domains="S2_1_LTRs_all_depth_LTR_cleaned.gff3 protein_match (1-based)")
    # LTR_retriever: TSD coordinates from the pass list (5_TSD, 3_TSD; 1-based closed)
    with open(RUN / "ltr_retriever/S2_1.fa.mod.pass.list") as fh:
        pl = [x.rstrip("\n").split("\t") for x in fh if not x.startswith("#")]
    for f in F["ltr_retriever"]:
        for r in pl:
            ch, rng = r[0].split(":")
            a, b = sorted(int(v) for v in rng.split(".."))
            if ch == seq and (a - 1, b) == (f["s"], f["e"]):
                t5, t3 = ([int(v) for v in x.split("..")] for x in (r[4], r[5]))
                f["tsd"] = [(t5[0] - 1, t5[1]), (t3[0] - 1, t3[1])]
    SRC["ltr_retriever"] = dict(strand="pass.list Strand (calls tsv)", LTR="pass.list (calls tsv)",
                                TSD="S2_1.fa.mod.pass.list 5_TSD / 3_TSD (1-based)", domains="none")
    # LTRharvest: strand '?', TSD features in its GFF3 (seq0 = first FASTA record = c0001)
    g = gff_rows(RUN / "ltrharvest/ltrharvest.gff3", "seq0", w0, w1)
    for f in F["ltrharvest"]:
        rr = next(r for r in g if r["type"] == "LTR_retrotransposon" and (r["s"], r["e"]) == (f["s"], f["e"]))
        par = rr["att"]["Parent"]
        f["tsd"] = [(r["s"], r["e"]) for r in g if r["type"] == "target_site_duplication" and r["att"]["Parent"] == par]
    SRC["ltrharvest"] = dict(strand="not reported ('?')", LTR="ltrharvest.tabout (calls tsv)",
                             TSD="ltrharvest.gff3 target_site_duplication (1-based)", domains="none")
    # LTR_Finder: strand, TSR sequence and domain columns (RT, IN(core), IN(c-term), RH; 1-based, N-N = none) in scn
    rows = {}
    with open(RUN / "ltr_finder/S2_1.fa.finder.combine.scn") as fh:
        for x in fh:
            c = x.rstrip("\n").split("\t")
            if len(c) > 12 and c[1] == seq:
                a, b = (int(v) for v in c[2].split("-"))
                rows[(a - 1, b)] = c
    for f in F["ltr_finder"]:
        c = rows[(f["s"], f["e"])]
        tsr = c[5]
        if tsr not in ("N", "NA", "") and set(tsr) <= set("ACGTN"):
            f["tsd"] = flank(dict(start=f["s"], end=f["e"]), len(tsr))
        for k, nm in ((8, "RT"), (9, "IN core"), (10, "IN c-term"), (11, "RH")):
            if c[k] != "N-N":
                a, b = (int(v) for v in c[k].split("-"))
                f["dom"].append((a - 1, b, nm))
    SRC["ltr_finder"] = dict(strand="scn Strand", LTR="scn Location, LTR len",
                             TSD="scn TSR sequence (flanks drawn adjacent)",
                             domains="scn RT / IN(core) / IN(c-term) / RH columns: N-N for all three calls here")
    # Look4LTRs: strand from RC (README: 1 = reverse complemented), TSD coordinates (0-based) in the rtr table
    with open(RUN / "look4ltrs/out/Rtr/c0001.rtr") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        rt = [dict(zip(hdr, x.rstrip("\n").split("\t"))) for x in fh]
    for f in F["look4ltrs"]:
        r = next(x for x in rt if x["RightEnd"] != "NA" and (int(x["LeftStart"]), int(x["RightEnd"])) == (f["s"], f["e"]))
        if r["TSDStart"] != "NA":
            f["tsd"] = [(int(r["TSDStart"]), f["s"]), (f["e"], int(r["TSDEnd"]))]
    SRC["look4ltrs"] = dict(strand="rtr RC (calls tsv)", LTR="rtr LeftEnd / RightStart (0-based)",
                            TSD="out/Rtr/c0001.rtr TSDStart / TSDEnd (0-based; NA for the host)", domains="none")
    SRC["hirepeat"] = dict(strand="not reported", LTR="not reported", TSD="not reported", domains="none",
                           note="S2_1_all_rounds_original.bed: span only")
    # HiTE: TSD sequence only (intact_LTR.list TSD; coordinates NA..NA), strand '.'
    for f in F["hite"]:
        c = next(x for x in calls["hite"] if x["id"] == f["id"])
        if c["tsd"] not in (".", ""):
            f["tsd"] = flank(c, len(c["tsd"]))
    SRC["hite"] = dict(strand="not reported ('.')", LTR="intact_LTR.list IN (calls tsv)",
                       TSD="intact_LTR.list TSD sequence (flanks drawn adjacent)", domains="none")
    # LTRpred: TSD coordinates in the DataSheet (1-based), domains = LTRdigest protein_match of the same element
    with open(RUN / "ltrpred/S2_1_ltrpred/S2_1_LTRpred_DataSheet.tsv") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ds_ = [dict(zip(hdr, x.rstrip("\n").split("\t"))) for x in fh]
    g = gff_rows(RUN / "ltrpred/S2_1_ltrpred/S2_1_ltrdigest/S2_1_LTRdigestPrediction.gff", "seq0", w0, w1)
    for f in F["ltrpred"]:
        r = next(x for x in ds_ if x["chromosome"] == seq and (int(x["start"]) - 1, int(x["end"])) == (f["s"], f["e"]))
        f["tsd"] = [(int(r["lTSD_start"]) - 1, int(r["lTSD_end"])), (int(r["rTSD_start"]) - 1, int(r["rTSD_end"]))]
        pid = r["ID"].replace("S2_1_", "")
        f["dom"] = [(x["s"], x["e"], x["att"]["name"]) for x in g
                    if x["type"] == "protein_match" and x["att"].get("Parent") == pid]
    SRC["ltrpred"] = dict(strand="DataSheet strand", LTR="DataSheet lLTR / rLTR (1-based)",
                          TSD="S2_1_LTRpred_DataSheet.tsv lTSD / rTSD (1-based)",
                          domains="S2_1_LTRdigestPrediction.gff protein_match of the element (1-based)")
    # DANTE_LTR: TSD and protein_domain children of the call in dante_ltr.gff3 (1-based)
    g = gff_rows(RUN / "dante_ltr/dante_ltr.gff3", seq, w0, w1)
    for f in F["dante_ltr"]:
        f["tsd"] = [(r["s"], r["e"]) for r in g if r["type"] == "target_site_duplication"
                    and r["att"].get("Parent") == f["id"]]
        f["dom"] = [(r["s"], r["e"], r["att"]["Name"]) for r in g if r["type"] == "protein_domain"
                    and r["att"].get("Parent") == f["id"]]
    SRC["dante_ltr"] = dict(strand="dante_ltr.gff3 transposable_element", LTR="dante_ltr.gff3 long_terminal_repeat",
                            TSD="dante_ltr.gff3 target_site_duplication (1-based)",
                            domains="dante_ltr.gff3 protein_domain children (1-based)")
    for t in TOOLS:                                # sanity: features lie on or beside their call
        for f in F[t]:
            for a, b in (f["tsd"] or []):
                assert b <= f["s"] + 1 or a >= f["e"] - 1, (t, f["id"], a, b)
                assert 2 <= b - a <= 8, (t, f["id"], a, b)
            for a, b, _ in f["dom"]:
                assert f["s"] <= a < b <= f["e"], (t, f["id"], a, b)
    return F, SRC


def load_example():
    """A: three true elements of S2 seed 1 (host, its tenant, a non-nested element 16.7 kb away with no other truth
    LTR-RT between), every tool's calls matched (score_elements.match, RO 0.9, one-to-one), the structural checks of
    fig1_nested_sims extended to the non-nested element (asserted), and each tool's reported features"""
    ds, chrom = EX["ds"], EX["chrom"]
    tr = se.read_truth(str(B / "sims/truth" / f"{ds}.truth.tsv"))
    calls, c2t, found = {}, {}, {}
    for t in TOOLS:
        c = se.read_calls(str(B / "calls" / ds / f"{t}.tsv"))
        m, t2c = se.match(c, tr, 0.9)
        calls[t], c2t[t] = c, m
        found[t] = t2c
    tid = {x["id"]: i for i, x in enumerate(tr)}
    hi, ti, ni = tid[EX["host"]], tid[EX["tenant"]], tid[EX["nonnested"]]
    h, u, nn = tr[hi], tr[ti], tr[ni]
    assert nn["class"] == "intact" and nn["chrom"] == chrom and 0 < nn["start"] - h["end"] <= 20_000
    between = [x for x in tr if x["chrom"] == chrom and x["start"] < nn["start"] and x["end"] > h["end"]]
    assert not between, between                   # no other truth LTR-RT between host and non-nested element
    pieces = [tuple(map(int, p.split("-"))) for p in h["pieces"].split(",")]
    (p1s, p1e), (p2s, p2e) = pieces
    tsdlen = len(u["tsd"])
    fa = fa_path(ds)
    chk = {}
    chk["piece1_end_eq_tenant_start"] = p1e == u["start"]
    chk["piece2_start_eq_tenant_end_plus_tsd"] = p2s == u["end"] + tsdlen
    for nm, x in (("host", h), ("tenant", u), ("nonnested", nn)):
        L = int(x["ltr_len"])
        tl = fetch(fa, chrom, x["start"] - len(x["tsd"]), x["start"])
        tr_ = fetch(fa, chrom, x["end"], x["end"] + len(x["tsd"]))
        chk[f"{nm}_tsd_identical"] = tl == tr_ == x["tsd"]
        L1, L2 = fetch(fa, chrom, x["start"], x["start"] + L), fetch(fa, chrom, x["end"] - L, x["end"])
        chk[f"{nm}_ltrs_TG_CA"] = L1[:2] == "TG" and L1[-2:] == "CA" and L2[:2] == "TG" and L2[-2:] == "CA"
    ofa = fa_path(ds, EX["old_gen"])
    old = fetch(ofa, chrom, *EX["old_host"])
    now = fetch(fa, chrom, p1s, p1e) + fetch(fa, chrom, p2s, p2e)
    chk["old_host_len_eq_pieces"] = len(old) == len(now) == int(h["seq_len"])
    chk["old_single_tsd_at_site"] = old[p1e - p1s - tsdlen:p1e - p1s] == u["tsd"]
    for k, x in chk.items():
        log(f"  example check {k}: {x}")
    assert all(chk.values()), chk
    w0, w1 = EX["win"]
    wc = {}
    for t in TOOLS:
        raw = {r["id"]: r for r in read_tsv(B / "calls" / ds / f"{t}.tsv")
               if r["chrom"] == chrom and int(r["end"]) > w0 and int(r["start"]) < w1}
        out = []
        for k, c in enumerate(calls[t]):
            if c["chrom"] == chrom and c["end"] > w0 and c["start"] < w1:
                m = c2t[t].get(k)
                r = raw[c["id"]]
                out.append(dict(c, matched=None if m is None else tr[m]["id"],
                                strand=r["strand"], tsd=r["tsd"],
                                lltr=int(r["lltr_end"]) if r["lltr_end"] not in (".", "") else None,
                                rltr=int(r["rltr_start"]) if r["rltr_start"] not in (".", "") else None))
        wc[t] = out
    for c in (x for t in TOOLS for x in wc[t]):
        assert c["matched"] is not None, "every call in the window matches a true element"
    found_by = [t for t in TOOLS if hi in found[t]]
    tfound_by = [t for t in TOOLS if ti in found[t]]
    nfound_by = [t for t in TOOLS if ni in found[t]]
    assert "ltrquest" in found_by and "ltrquest" in tfound_by and "ltrquest" in nfound_by, "no LTRquest miss shown"
    # how typical: other tools finding each host / each non-nested element of this genome
    hosts = [i for i, x in enumerate(tr) if x["class"] == "host"]
    intact = [i for i, x in enumerate(tr) if x["class"] == "intact"]
    oth = lambda i: sum(i in found[t] for t in TOOLS[1:])     # noqa: E731
    med_h, med_n = float(np.median([oth(i) for i in hosts])), float(np.median([oth(i) for i in intact]))
    F, SRC = reported_features(wc, w0, w1)
    num("A", "example host (S2 seed 1, 0-based half-open)", f"{h['id']} {h['family']} {h['superfamily']} "
        f"{h['strand']}", "sims/truth/S2_1.truth.tsv")
    num("A", "example tenant", f"{u['id']} {u['family']} {u['superfamily']} {u['strand']}", "sims/truth/S2_1.truth.tsv")
    num("A", "example non-nested element", f"{nn['id']} {nn['family']} {nn['superfamily']} {nn['strand']} "
        f"(LTR {nn['ltr_len']} bp, TSD {nn['tsd']}, {nn['start'] - h['end']} bp from the host)",
        "sims/truth/S2_1.truth.tsv")
    num("A", "tools matching host / tenant / non-nested element",
        f"{','.join(NAME[t] for t in found_by)} / {','.join(NAME[t] for t in tfound_by)} / "
        f"{','.join(NAME[t] for t in nfound_by)}", "score_elements.match RO 0.9")
    num("A", "other tools finding the host / non-nested element here; median over all S2_1 hosts / non-nested",
        f"{len(found_by) - 1} / {len(nfound_by) - 1}; {med_h:g} / {med_n:g}", "score_elements.match RO 0.9")
    for t in TOOLS:
        num("A", f"{NAME[t]} reported features (source)", "; ".join(f"{k}: {v}" for k, v in SRC[t].items()),
            f"runs/S2_1/{t}/ and calls/S2_1/{t}.tsv")
        for f in F[t]:
            num("A", f"{NAME[t]} call {f['s']}-{f['e']} (0-based): strand / LTR ends / TSDs / domains",
                f"{f['strand'] or 'not reported'} / {f['lltr']},{f['rltr']} / {f['tsd'] or 'not reported'} / "
                f"{', '.join(f'{d[2]} {d[0]}-{d[1]}' for d in f['dom']) or 'none'}", f"runs/S2_1/{t}/")
    (a1, b1), (a2, b2) = EX["segs"]
    om = (b1, a2)                                  # omitted stretch of A's broken axis: no truth element, no call
    assert not [x for x in tr if x["chrom"] == chrom and x["end"] > om[0] and x["start"] < om[1]]
    assert not [(t, c["id"]) for t in TOOLS for c in wc[t] if c["end"] > om[0] and c["start"] < om[1]]
    for x in (h, u, nn):
        assert a1 <= x["start"] - 300 and x["end"] + 300 <= b1 or a2 <= x["start"] - 300 and x["end"] + 300 <= b2
    num("A", "broken axis: segments drawn (0-based) and omitted stretch (no truth element, no call of any tool)",
        f"chr1:{a1}-{b1} | chr1:{a2}-{b2}; omitted {om[1] - om[0]} bp", "this script (asserted)")
    return dict(h=h, u=u, nn=nn, pieces=pieces, calls=wc, F=F, SRC=SRC, found_by=found_by, tfound_by=tfound_by,
                nfound_by=nfound_by, med_h=med_h, med_n=med_n)


# ============================================================================ data: B (one scorecard, two class systems)
BSYS = [("nesting", "nesting", [("unnested", "non-nested"), ("tenant", "tenants"), ("host", "hosts")]),
        ("protein", "protein domains", [("coding", "coding"), ("proteinless", "protein-less")])]
BCLS = [(s, c, lab) for s, _, cl in BSYS for c, lab in cl]
CMET = [("precision", "prec."), ("recall", "rec."), ("f1", "F1")]


def load_B():
    """S2 (3 deeply nested genomes), mean of 3 genomes: precision, recall, F1 per class of both class systems on one
    basis (experiments/figpub/sim_class_prf_joint_summary.tsv, ext truth = every truth.tsv element,
    one-to-one RO >= 0.9; nesting classes as sim_class_prf; protein class of a true element = its library entry's
    class, of a call = the domains inside it). NA = the caller makes no call of that class (precision, F1 undefined;
    drawn "–"). Checks: nesting rows equal sim_class_prf_summary.tsv (printed strings); every finite value comes from
    3 genomes; NA only for precision / F1 with 0 genomes; n_truth per class equal for every tool."""
    want = {(s, c) for s, c, _ in BCLS}
    rows = [r for r in read_tsv(FP / "sim_class_prf_joint_summary.tsv")
            if r["set"] == "S2" and r["tool"] in TOOLS and (r["system"], r["cls"]) in want]
    assert len(rows) == len(TOOLS) * len(want), len(rows)
    old = {(r["tool"], r["cls"]): r for r in read_tsv(FP / "sim_class_prf_summary.tsv") if r["set"] == "S2"}
    src = "experiments/figpub/sim_class_prf_joint_summary.tsv set=S2"
    V, SD, n = {}, {}, {}
    for r in rows:
        t, c = r["tool"], r["cls"]
        for m, _ in CMET:
            v = r[f"{m}_mean"]
            if v == "NA":
                assert m != "recall" and r[f"{m}_n_genomes"] == "0" and int(r["n_calls"]) == 0, (t, c, m)
                V[(c, m, t)] = None
                num("B", f"S2 {NAME[t]} {c} {m}", "NA (no calls of this class; drawn '–')",
                    f"{src} tool={t} cls={c} {m}_mean (n_calls 0)")
                continue
            assert r[f"{m}_n_genomes"] == "3", (t, c, m)
            V[(c, m, t)], SD[(c, m, t)] = float(v), float(r[f"{m}_sd"])
            num("B", f"S2 {NAME[t]} {c} {m} (mean of 3 genomes; SD {float(r[m + '_sd']):.4f})", f"{float(v):.4f}",
                f"{src} tool={t} system={r['system']} cls={c} {m}_mean, {m}_sd")
        if r["system"] == "nesting":                 # the nesting half = as in sim_class_prf
            o = old[(t, c)]
            for m, _ in CMET:
                assert (r[f"{m}_mean"], r[f"{m}_sd"]) == (o[f"{m}_mean"], o[f"{m}_sd"]), (t, c, m)
        assert n.setdefault(c, int(r["n_truth"])) == int(r["n_truth"]), (t, c)
    for s, c, lab in BCLS:
        num("B", f"true elements in class {lab} (S2, 3 genomes summed; header n)", n[c],
            f"{src} n_truth (system {s}, cls {c})")
    # spread across the three genomes (one legend clause): largest SD per tool over the drawn classes and metrics
    sd_max = {t: max(v for (c, m, tt), v in SD.items() if tt == t) for t in TOOLS}
    for t in TOOLS:
        num("B", f"S2 {NAME[t]} largest across-genome SD over the drawn classes and metrics", f"{sd_max[t]:.4f}",
            f"{src} *_sd")
    return dict(V=V, SD=SD, n=n, sd_max=sd_max)

# ============================================================================ data: C, D (aging series)
def k2p_of_identity(p100):
    """LTR-LTR identity (%) -> K2P under the simulation's model (substitutions only, transitions = transversions)"""
    p = 1 - p100 / 100
    P = Q = p / 2
    return -0.5 * math.log(1 - 2 * P - Q) - 0.25 * math.log(1 - 2 * Q)


FLOOR = {"ltrquest": None, "ltr_retriever": 85, "ltrharvest": 85, "ltr_finder": 85, "look4ltrs": None,
         "hirepeat": None, "hite": 85, "ltrpred": 70, "dante_ltr": 70}


def pooled_d50(xs, rs):
    """K2P at which the pooled recall falls to half its timepoint-0 value: scripts/div_score.d50 (linear between
    timepoints), on the pooled curve"""
    return div_score.d50([(t, float(x), float(r)) for t, (x, r) in enumerate(zip(xs, rs))])


def load_curve(cls, panel):
    """one protein class of the aging series (C: coding, D: protein-less), 3 genomes pooled: experiments/figpub/
    div_all_ci.py's procedure on one class (as for the protein-less panel). Per tool x timepoint recall over the
    pooled elements of the class; 95% CI from 1,000 resamples of their library families (numpy default_rng(1), every
    copy of a drawn family carried along, percentiles 2.5 / 97.5 as nonauto_stats.ci); x = median K2P of the pooled
    elements of the class at the timepoint (each element once); one dot at the pooled recall-halved point
    (scripts/div_score.d50, linear between timepoints). Checks: pooled n / tp = div_class.tsv sums over seeds; the
    pooled d50 agrees with div_d50.tsv (same pool, x = class median over seeds) within 0.002 and lies inside its CI."""
    de = pd.read_csv(NA / "div_elements.tsv", sep="\t",
                     usecols=["ds", "timepoint", "tool", "element", "family", "cls", "k2p", "hit"])
    s_ = de[(de.cls == cls) & de.tool.isin(TOOLS)]
    x = np.array([statistics.median({(r.ds, r.element): r.k2p
                                     for r in s_[s_.timepoint == t].itertuples()}.values()) for t in range(11)])
    dc = pd.read_csv(NA / "div_class.tsv", sep="\t")
    dc = dc[(dc.cls == cls) & dc.tool.isin(TOOLS)]
    ref = pd.read_csv(NA / "div_d50.tsv", sep="\t", dtype=str, keep_default_na=False)
    fams = sorted(s_.family.unique())
    fi = {f: i for i, f in enumerate(fams)}
    nF = len(fams)
    draws = np.random.default_rng(1).integers(0, nF, size=(1000, nF))
    M = np.stack([np.bincount(z, minlength=nF) for z in draws]).astype(float)
    lab = {"coding": "coding", "proteinless": "protein-less"}[cls]
    S = {}
    for t in TOOLS:
        st = s_[s_.tool == t]
        N, Hh = np.zeros((nF, 11)), np.zeros((nF, 11))
        f_idx = st.family.map(fi).to_numpy()
        np.add.at(N, (f_idx, st.timepoint.to_numpy()), 1)
        np.add.at(Hh, (f_idx, st.timepoint.to_numpy()), st.hit.to_numpy())
        n, h = N.sum(0), Hh.sum(0)
        g = dc[dc.tool == t].groupby("timepoint")[["n", "tp"]].sum()
        if not (np.array_equal(g.n.to_numpy(), n) and np.array_equal(g.tp.to_numpy(), h)):
            raise SystemExit(f"{STEM}: {panel} pooled {cls} counts {t} differ from div_class.tsv")
        r = h / n
        with np.errstate(invalid="ignore", divide="ignore"):
            bs = (M @ Hh) / (M @ N)
        lo, hi = np.nanpercentile(bs, 2.5, axis=0), np.nanpercentile(bs, 97.5, axis=0)
        v = pooled_d50(x, r)
        rr = ref[(ref.tool == t) & (ref.cls == cls)].iloc[0]
        if v[0].isdigit() and rr.d50[:1].isdigit():
            assert abs(float(v) - float(rr.d50)) <= 0.002, (t, cls, v, rr.d50)
            assert float(rr.d50_lo) <= float(v) <= float(rr.d50_hi), (t, cls, v, rr.d50_lo, rr.d50_hi)
        S[t] = dict(x=x, r=r, lo=lo, hi=hi, d50=v, r0=r[0], ref=rr.d50, ref_lo=rr.d50_lo, ref_hi=rr.d50_hi)
        num(panel, f"{NAME[t]} {lab} pooled d50 (K2P at which pooled recall halves; dot)", v,
            f"computed: scripts/div_score.d50 on pooled {cls} counts (experiments/nonauto/div_elements.tsv)")
        num(panel, f"{NAME[t]} {lab} d50 in div_d50.tsv [95% CI] (x = class median over seeds), for comparison",
            f"{rr.d50} [{rr.d50_lo}, {rr.d50_hi}]", "experiments/nonauto/div_d50.tsv")
        for tt in range(11):
            num(panel, f"{NAME[t]} timepoint {tt} {lab} pooled recall [95% CI] at median K2P {x[tt]:.4f}",
                f"{r[tt]:.4f} [{lo[tt]:.4f}, {hi[tt]:.4f}] (n {int(n[tt])})",
                f"computed from experiments/nonauto/div_elements.tsv cls={cls} (= div_class.tsv sums)")
    num(panel, f"{lab} element families (unit of the resampling)", nF, "experiments/nonauto/div_elements.tsv")
    num(panel, f"{lab} elements per timepoint (3 genomes pooled)", int(s_[(s_.tool == TOOLS[0]) &
                                                                           (s_.timepoint == 0)].shape[0]),
        "experiments/nonauto/div_elements.tsv")
    if panel == "C":
        for fl in (85, 70):
            num("C", f"documented identity floor {fl}% as K2P (simulation model)", f"{k2p_of_identity(fl):.3f}",
                "computed (substitutions only, transitions = transversions); fig2_divergence memo")
    return dict(S=S, n_fam=nF)


# ============================================================================ data: E, F (real genomes)
E_DS = ["arab", "rice", "finch"]                  # E: masking P / R / F1 (maize and gorilla are in Fig. S2C)
EMET = [("prec", "precision", "prec."), ("sens", "recall", "rec."), ("F1", "f1", "F1")]
# missing cells as drawn, asserted per panel (reasons: resources.tsv status / masking.tsv 'incomplete')
E_MISSING = {"failed": {("finch", "ltr_retriever"), ("finch", "dante_ltr")},
             "noteval": {("finch", "hirepeat")}}
F_MISSING = {"failed": {("finch", "ltr_retriever"), ("maize", "ltr_retriever"), ("gorilla", "ltr_retriever"),
                        ("gorilla", "ltrharvest"), ("gorilla", "ltrpred"), ("gorilla", "look4ltrs"),
                        ("finch", "dante_ltr"), ("gorilla", "dante_ltr")},
             "noteval": {("gorilla", "hirepeat")}}


def missing_reasons():
    """(ds, tool) -> ("failed" | "noteval", source). failed = the tool did not finish within 36 h (timeout), exceeded
    200 GB (mem_guard) or crashed with no calls (failed): resources.tsv status. noteval = the tool finished but our
    masking or library build with its library would have run past 36 h or was stopped at 36 h: masking.tsv ci
    'incomplete: ...'. Every (ds, tool) gets at most one class (asserted)."""
    res = pd.read_csv(TAB / "resources.tsv", sep="\t")
    m = pd.read_csv(TAB / "masking.tsv", sep="\t")
    out = {}
    for ds in DS:
        for t in TOOLS:
            r = res[(res.ds == ds) & (res.tool == t)]
            st = r.status.iloc[0] if len(r) else "?"
            row = m[(m.ds == ds) & (m.lib == f"{t}.intact") & (m.role == "primary")]
            ci = str(row.ci.iloc[0]) if len(row) else ""
            failed = st in ("timeout", "mem_guard", "failed")
            noteval = len(row) and row.F1.iloc[0] == "." and ci.startswith("incomplete")
            assert not (failed and noteval), (ds, t)
            if failed:
                out[ds, t] = ("failed", f"resources.tsv status={st}")
            elif noteval:
                out[ds, t] = ("noteval", f"masking.tsv ci: {ci}")
    return out


def check_missing(panel, dss, have, reasons, expected):
    """every missing cell of a panel has a class, no failed tool has a value, and the drawn sets equal the lists"""
    got = {"failed": set(), "noteval": set()}
    for ds in dss:
        for t in TOOLS:
            if (ds, t) in have:
                if (ds, t) in reasons and reasons[(ds, t)][0] == "failed":
                    raise SystemExit(f"{STEM}: {panel} result present for a failed tool {ds} {t}")
                continue
            if (ds, t) not in reasons:
                raise SystemExit(f"{STEM}: no class for the missing {panel} result {ds} {t}")
            got[reasons[(ds, t)][0]].add((ds, t))
    assert got == expected, (panel, got)
    for kind in ("failed", "noteval"):
        for ds, t in sorted(got[kind]):
            num(panel, f"{ds} {NAME[t]} missing result class (cell text)",
                "tool failed" if kind == "failed" else "not evaluated", reasons[(ds, t)][1])


def load_EF():
    """E: masking precision, recall (lib-test sensitivity), F1 of each tool's intact library against curated
    annotation (report/tables/masking.tsv role primary, <tool>.intact: prec, sens, F1) for A. thaliana, rice, zebra
    finch. F: share of curated LTR entries each library recovers perfect or good (report/tables/bm_rm2.tsv), 5
    genomes. Missing cells classed and asserted per panel."""
    m = pd.read_csv(TAB / "masking.tsv", sep="\t")
    m = m[(m.role == "primary") & (m.F1 != ".")].copy()
    m["tool"] = m.lib.str.replace(".intact", "", regex=False)
    m = m[m.tool.isin(TOOLS)]
    mk = {(r.ds, r.tool): dict(precision=float(r.prec), recall=float(r.sens), f1=float(r.F1), f1_lo=r.F1_lo,
                               f1_hi=r.F1_hi) for r in m.itertuples() if r.ds in E_DS}
    r = pd.read_csv(TAB / "bm_rm2.tsv", sep="\t")
    r = r[(r.role == "primary") & r.lib.str.endswith(".intact")].copy()
    r["tool"] = r.lib.str.replace(".intact", "", regex=False)
    r = r[r.tool.isin(TOOLS)]
    pg = {(x.ds, x.tool): float(x.frac_perfect_or_good) for x in r.itertuples()}
    ncur = r.groupby("ds").n_curated_LTR.first().to_dict()
    reasons = missing_reasons()
    check_missing("E", E_DS, mk, reasons, E_MISSING)
    check_missing("F", DS, pg, reasons, F_MISSING)
    for (ds, t), v in sorted(mk.items()):
        for col, key, _ in EMET:
            extra = f" [95% CI {v['f1_lo']}, {v['f1_hi']}]" if key == "f1" else ""
            num("E", f"{ds} {NAME[t]} masking {key}{extra}", f"{v[key]:.4f}",
                f"report/tables/masking.tsv {col} (role primary, {t}.intact)")
    for (ds, t), v in pg.items():
        num("F", f"{ds} {NAME[t]} share of curated LTR entries perfect or good", f"{v:.4f}",
            "report/tables/bm_rm2.tsv frac_perfect_or_good")
    for ds in DS:
        num("F", f"{ds} curated LTR entries (n)", int(ncur[ds]), "report/tables/bm_rm2.tsv n_curated_LTR")
    return dict(mk=mk, pg=pg, ncur=ncur, reasons=reasons, m=m)


# ============================================================================ data: G (A. thaliana loci)
def load_coding(ds, tool):
    out = []
    with open(B / f"eval/adj/{ds}/{tool}.calls.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["coding"] != "1":
                continue
            out.append((r["chrom"], int(r["start"]), int(r["end"]), r["k2p"], int(r["int_start"]), int(r["int_end"]),
                        r["supported"] == "1", r["null_shift"], r["null_supported"], r["tsd"] == "1",
                        r["tgca"] == "1"))
    return out


class Cover:
    """curated LTR coverage from merged intervals (proxies.ref_ltr: lib-test LTR pattern, 0-based half-open)"""

    def __init__(self, out_path):
        self.iv = {}
        for c, v in proxies.ref_ltr(str(out_path)).items():
            a = np.array(v, dtype=np.int64)
            ln = a[:, 1] - a[:, 0]
            self.iv[c] = (a[:, 0], ln, np.concatenate([[0], np.cumsum(ln)]))

    def cum(self, c, x):
        x = np.asarray(x, dtype=np.int64)
        if c not in self.iv:
            return np.zeros(x.shape)
        s, ln, cs = self.iv[c]
        i = np.searchsorted(s, x, side="right") - 1
        ok = i >= 0
        ii = np.where(ok, i, 0)
        v = cs[ii] + np.clip(x - s[ii], 0, ln[ii])
        return np.where(ok, v, 0).astype(float)


def edges(s, e, a, b):
    return np.round(np.concatenate([np.linspace(s, a, NB5 + 1)[:-1], np.linspace(a, b, NBI + 1)[:-1],
                                    np.linspace(b, e, NB3 + 1)])).astype(np.int64)


def profile(cov, c, ed):
    v = cov.cum(c, ed)
    w = np.diff(ed).astype(float)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(w > 0, np.diff(v) / np.maximum(w, 1), np.nan)


def prep_G(ds="arab"):
    """copied from figS3_library_age.prep (then figS5_library_age): coding-unit loci pooled over callers (adj_summary rule), a coverage
    profile of a representative member and of its moved span; checked against library_age_bias.tsv"""
    pool = [t for t in TOOLS if (B / f"eval/adj/{ds}/{t}.calls.tsv").exists()]
    rows = {t: load_coding(ds, t) for t in pool}
    mem = [(r[0], r[1], r[2], (t, k)) for t in pool for k, r in enumerate(rows[t])]
    lab, n = su.loci(mem)
    by = defaultdict(list)
    for c, s, e, key in mem:
        by[lab[key]].append(key)
    cov = Cover(B / f"mask/{ds}/curated/{ds}.fa.out")
    k2p, sup, nul, prof, nprof = [], [], [], [], []
    for lid in sorted(by):
        rs = [rows[t][k] for t, k in by[lid]]
        s_ = sum(r[6] for r in rs) >= len(rs) / 2
        ks = [float(r[3]) for r in rs if r[3] != "NA"]
        km = statistics.median(ks) if ks else None
        nv = [r[8] == "1" for r in rs if r[7] != "NA" and r[8] in ("0", "1")]
        nl = (sum(nv) >= len(nv) / 2) if nv else None
        r = rs[0] if km is None else min(rs, key=lambda x: abs(float(x[3]) - km) if x[3] != "NA" else math.inf)
        c, s, e, a, b = r[0], r[1], r[2], r[4], r[5]
        if not (s < a < b < e):
            a, b = s + (e - s) // 10, e - (e - s) // 10
        ed = edges(s, e, a, b)
        p = profile(cov, c, ed)
        q = profile(cov, c, ed + int(r[7])) if r[7] != "NA" else np.full(NBIN, np.nan)
        k2p.append(np.nan if km is None else km)
        sup.append(s_)
        nul.append(np.nan if nl is None else float(nl))
        prof.append(p)
        nprof.append(q)
    k2p, sup = np.array(k2p), np.array(sup)
    np.savez_compressed(TMP / f"prep_{ds}.npz", k2p=k2p, sup=sup, null=np.array(nul),
                        prof=np.array(prof, np.float32), nprof=np.array(nprof, np.float32))


def load_G(redo=False, ds="arab"):
    if redo or not (TMP / f"prep_{ds}.npz").exists():
        log("G: building A. thaliana loci and coverage profiles (cached in tmp/)", True)
        prep_G(ds)
    z = np.load(TMP / f"prep_{ds}.npz")
    d = {k: z[k] for k in z.files}
    tab = {r["age_bin"]: r for r in read_tsv(TAB / "library_age_bias.tsv") if r["ds"] == ds}
    k2p, sup = d["k2p"], d["sup"]
    for i, (lo, hi) in enumerate([(0, .01), (.01, .03), (.03, .1), (.1, math.inf)]):
        m = (k2p >= lo) & (k2p < hi)
        t = tab[BINS[i]]
        if int(t["n_loci"]) != int(m.sum()) or abs(float(t["supported_share"]) - sup[m].mean()) > 5e-4:
            raise SystemExit(f"{STEM}: G {BINS[i]}: {m.sum()} / {sup[m].mean():.4f} do not reproduce "
                             f"library_age_bias.tsv ({t['n_loci']} / {t['supported_share']})")
        num("G", f"A. thaliana protein-confirmed loci, K2P {BINS[i]}: n, library-supported share",
            f"{t['n_loci']}, {t['supported_share']} [{t['wilson95']}]", "report/tables/library_age_bias.tsv")
    gap = {(r["ds"], r["set"]): r for r in read_tsv(B / "experiments/a4_diagnostics/annotation_gap.tsv")}
    d["gap_mb"] = float(gap[(ds, "all")]["gap_mb"])
    num("G", "moved spans: each locus shifted 20-200 kb (chance)", "20-200 kb",
        "eval/adj/arab/<tool>.calls.tsv null_shift (approach 4 null; figS3_library_age memo)")
    d["tab"] = tab
    num("G", "A. thaliana protein-confirmed bases outside the curated annotation (Mb, all calls)", d["gap_mb"],
        "experiments/a4_diagnostics/annotation_gap.tsv gap_mb (arab, all)")
    return d


# ============================================================================ data: H (ltrkin)
HDS = ("camara_hap2", "camara_hap1")
HAPN = {"camara_hap2": "2", "camara_hap1": "1"}
HCOLS = [("all", None, None, "all"),
         ("nesting", "nesting", "unnested", "non-nested"), ("nesting", "nesting", "tenant", "tenants"),
         ("nesting", "nesting", "host", "hosts"),
         ("protein", "protein_class", "coding", "coding"), ("protein", "protein_class", "proteinless", "protein-less")]


def load_H():
    """ltrkin, real insertion polymorphisms of the two Cardamine amara haplotypes. All loci: share genotyped correctly
    and its block-bootstrap 95% CI (experiments/figpub/ltrkin_outcomes_all.tsv; = ltrkin_yield.tsv yield, asserted).
    Per class (non-nested, tenants, hosts | coding, protein-less; partly coding loci only in Fig. S4) and caller, from
    the per-locus table experiments/figpub/ltrkin_locus_class.tsv: n loci, genotyped correctly (out_<tool> = correct)
    split into loci also nominated by another caller (solid) and loci this caller alone nominated (nominators = this
    caller only; hatched). Checks: per caller the locus table's all-loci n and correct = ltrkin_outcomes_all.tsv;
    LTRquest's alone-nominated loci and their correct calls = n_ltrquest_only_loci / n_correct_ltrquest_only (and the
    table's ltrquest_only flag); a caller that did not finish has '.' at every locus. Also every caller's share on the
    loci two or more callers nominated (legend)."""
    rows = {(r["ds"], r["tool"]): r for r in read_tsv(FP / "ltrkin_outcomes_all.tsv")
            if r["role"] == "primary" and r["tool"] in TOOLS}
    Y = {(r["ds"], r["tool"]): r for r in read_tsv(LK / "report/tables/ltrkin_yield.tsv") if r["run"] == "e2e"}
    d = pd.read_csv(FP / "ltrkin_locus_class.tsv", sep="\t", dtype=str, keep_default_na=False)
    src = "experiments/figpub/ltrkin_locus_class.tsv"
    out, ncls, core = {}, {}, {}
    for ds in HDS:
        g = d[d.ds == ds].reset_index(drop=True)
        noms = g.nominators.str.split(",")
        nnom = noms.str.len()
        assert ((noms.apply(lambda z: z == ["ltrquest"])) == (g.ltrquest_only == "1")).all(), ds
        masks = {lab: (g.index >= 0) if col is None else (g[col] == val) for _, col, val, lab in HCOLS}
        for _, _, _, lab in HCOLS:
            ncls[(ds, lab)] = int(masks[lab].sum())
        core[ds] = dict(n=int((nnom >= 2).sum()), n_single=int((nnom == 1).sum()))
        num("H", f"{ds} loci nominated by two or more callers / by one caller only",
            f"{core[ds]['n']} / {core[ds]['n_single']}", f"{src} nominators")
        for t in TOOLS:
            r = rows[(ds, t)]
            o = g[f"out_{t}"]
            if r["status"] != "ok":
                assert r["status"].startswith("DNF") and (o == ".").all(), (ds, t, r["status"])
                out[(ds, t)] = dict(status=r["status"])
                num("H", f"{ds} {NAME[t]} status", r["status"], "experiments/figpub/ltrkin_outcomes_all.tsv")
                continue
            n, nc = int(r["n_loci"]), int(r["n_correct"])
            assert abs(float(r["share_correct"]) - float(Y[(ds, t)]["yield"])) < 1e-9, (ds, t)   # = ltrkin_yield
            corr = (o == "correct").to_numpy()
            alone = noms.apply(lambda z, t=t: z == [t]).to_numpy()
            assert len(g) == n and int(corr.sum()) == nc, (ds, t, len(g), n, int(corr.sum()), nc)
            if t == "ltrquest":
                assert int(alone.sum()) == int(r["n_ltrquest_only_loci"]), ds
                assert int((alone & corr).sum()) == int(r["n_correct_ltrquest_only"]), ds
            cls = {}
            for _, _, _, lab in HCOLS:
                mk = masks[lab].to_numpy() if hasattr(masks[lab], "to_numpy") else masks[lab]
                cls[lab] = dict(n=int(mk.sum()), solid=int((mk & corr & ~alone).sum()),
                                alone=int((mk & corr & alone).sum()))
            cc = (nnom >= 2).to_numpy()
            out[(ds, t)] = dict(status="ok", n=n, correct=nc, share=float(r["share_correct"]),
                                lo=float(r["correct_lo"]), hi=float(r["correct_hi"]), cls=cls,
                                n_alone=int(alone.sum()), alone_correct=int((alone & corr).sum()),
                                core_share=float((corr & cc).sum() / cc.sum()))
            num("H", f"{ds} {NAME[t]} all loci genotyped correctly [95% CI] (n of {n})",
                f"{float(r['share_correct']):.4f} [{float(r['correct_lo']):.4f}, {float(r['correct_hi']):.4f}] ({nc})",
                "experiments/figpub/ltrkin_outcomes_all.tsv share_correct, correct_lo, correct_hi (= ltrkin_yield.tsv)")
            for _, _, _, lab in HCOLS:
                c = cls[lab]
                num("H", f"{ds} {NAME[t]} {lab}: genotyped correctly, of them nominated by this caller alone (n loci)",
                    f"{(c['solid'] + c['alone']) / c['n']:.4f}: {c['solid'] + c['alone']}, {c['alone']} "
                    f"(n {c['n']})", f"{src} out_{t} = correct; nominators = {t} only")
            num("H", f"{ds} {NAME[t]} loci nominated by this caller alone (genotyped correctly)",
                f"{int(alone.sum())} ({int((alone & corr).sum())})", f"{src} nominators")
            num("H", f"{ds} {NAME[t]} share genotyped correctly on loci two or more callers nominated",
                f"{float((corr & cc).sum() / cc.sum()):.4f} ({int((corr & cc).sum())} of {int(cc.sum())})",
                f"{src} nominators, out_{t}")
        for _, _, _, lab in HCOLS:
            num("H", f"{ds} loci in class {lab} (header n)", ncls[(ds, lab)], f"{src} nesting / protein_class")
    return dict(o=out, n=ncls, core=core)


# ============================================================================ drawing helpers
def text_color(v):
    """white on dark cells, ink on light ones: relative luminance below 0.20 (equal contrast of white and ink)"""
    r, g, b = (c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4 for c in to_rgb(CMAP(v)))
    return "white" if 0.2126 * r + 0.7152 * g + 0.0722 * b < 0.20 else INK


def score_cell(ax, x, y, w, h, v, bold, pad=0.15, fs=FS):
    ax.add_patch(Rectangle((x + pad, y + pad), w - 2 * pad, h - 2 * pad, fc=CMAP(v), ec="none", lw=0))
    ax.text(x + w / 2, y + h / 2, f"{v:.2f}", ha="center", va="center", fontsize=fs, color=text_color(v),
            fontweight="bold" if bold else "normal")


def missing_cell(ax, x, y, w, h, kind, pad=0.15, fs=FS):
    """two styles: tool side "tool failed" (dark gray, white text); evaluation side "not evaluated" (light hatch)"""
    if kind == "failed":
        ax.add_patch(Rectangle((x + pad, y + pad), w - 2 * pad, h - 2 * pad, fc=FAILC, ec="none", lw=0))
        ax.text(x + w / 2, y + h / 2, "tool failed", ha="center", va="center", fontsize=fs, color="white")
    else:
        ax.add_patch(Rectangle((x + pad, y + pad), w - 2 * pad, h - 2 * pad, fc=MISSC, ec="none", lw=0, hatch=HATCH,
                               hatchcolor=HATCHC))
        ax.text(x + w / 2, y + h / 2, "not evaluated", ha="center", va="center", fontsize=fs, color=MUTED,
                bbox=dict(boxstyle="square,pad=0.05", fc=MISSC, ec="none"))


def best_set(vals):
    """tools whose printed (2-decimal) value equals the column's best printed value (ties all bold); no bold when
    fewer than two tools have a value in the column (same rule in B, E and F)"""
    pv = {t: float(f"{v:.2f}") for t, v in vals.items()}
    if len(pv) < 2:
        return set()
    b = max(pv.values())
    return {t for t, v in pv.items() if v == b}


def names(ax, x_right, ys, rows=TOOLS, fs=6.5):
    """tool names (LTRquest bold) and a square of the tool color, right-aligned to x_right (page mm)"""
    for y, t in zip(ys, rows):
        ax.text(x_right - 2.1, y, NAME[t], ha="right", va="center", fontsize=fs,
                fontweight="bold" if t == "ltrquest" else "normal")
        ax.add_patch(Rectangle((x_right - 1.5, y - 0.65), 1.3, 1.3, fc=COL[t], ec="none"))


def bp_axes(fig, x, y, w, h, xlim):
    """axes at (x, y, w, h) mm; data x = xlim, data y = page mm (downward)"""
    ax = poafig.mm_axes(fig, x, y, w, h)
    ax.set_xlim(*xlim)
    ax.set_ylim(y + h, y)
    ax.axis("off")
    return ax


def element(ax, s, e, y, h, body, ltr=None, lltr=None, rltr=None, ec="none", lw=0.0, z=2):
    ax.add_patch(Rectangle((s, y - h / 2), e - s, h, fc=body, ec=ec, lw=lw, zorder=z))
    if ltr is not None and lltr is not None:
        ax.add_patch(Rectangle((s, y - h / 2), lltr - s, h, fc=ltr, ec="none", zorder=z + 0.1))
    if ltr is not None and rltr is not None:
        ax.add_patch(Rectangle((rltr, y - h / 2), e - rltr, h, fc=ltr, ec="none", zorder=z + 0.1))


def chevron(ax, x, y, h, strand, col, bp_per_mm):
    d = 0.55 * bp_per_mm * (1 if strand == "+" else -1)
    ax.plot([x - d / 2, x + d / 2, x - d / 2], [y - h * 0.28, y, y + h * 0.28], color=col, lw=0.5,
            solid_capstyle="butt", solid_joinstyle="miter", zorder=5)


def clean(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="both", length=1.8, pad=1.2, labelsize=FS)


# ============================================================================ layout (mm from the top-left)
W = 180.0
G = dict(
    # band 1: nested simulations (A, B share the tool rows)
    XN=19.5,                       # right edge of the tool-name column (+ color square)
    XA0=20.5, XA1=76.0, XB0=80.5, XB1=179.5,
    TOP=6.5, TRUTH=13.5, RH=4.3,
    BCW=6.0, BGAP=1.6, BSYS=4.2,   # B: column width, gap between class groups, gap between the two class systems
    # band 2: aging series (C coding, key, D protein-less)
    Y2=71.0, XC0=20.5, XC1=82.5, XK0=86.0, XD0=115.5, XD1=177.5, HD=37.0,
    # band 3: real genomes (E, F share the tool rows)
    Y3=122.5, XE0=20.5, XE1=97.5, XF0=102.5, XF1=179.5, RH3=4.0, EGAP=1.6,
    # band 4: G (A. thaliana loci) and H (ltrkin)
    Y4=171.5, XG0=20.5, GBW=14.0, XH=53.0, XHN=73.0, XHB0=77.0, RH4=4.4,
    GHT=30.0,                      # G: heatmap height
    HALLW=20.0, HCLW=10.0, HGAP=1.2, HSYS=3.0,   # H: all-loci bar, class bars, gaps within / between systems
)
G["ROW0"] = G["TOP"] + G["TRUTH"]


def row_y(k):
    return G["ROW0"] + (k + 0.5) * G["RH"]


# ============================================================================ panels
def free_center(s0, e0, blocks):
    """center of the widest stretch of [s0, e0) not covered by blocks (for the strand chevron)"""
    cuts = sorted((max(a, s0), min(b, e0)) for a, b in blocks if b > s0 and a < e0)
    best, x = (0, (s0 + e0) / 2), s0
    for a, b in cuts + [(e0, e0)]:
        if a - x > best[0]:
            best = (a - x, (x + a) / 2)
        x = max(x, b)
    return best[1]


def panel_A(fig, can, ex):
    (a1, b1), (a2, b2) = EX["segs"]
    XL, XR = G["XA0"], G["XA1"]
    GAP = 3.0                                      # mm of the axis break
    bppm = ((b1 - a1) + (b2 - a2)) / (XR - XL - GAP)          # one bp scale for both segments
    gap_bp = GAP * bppm

    def X(v):
        """genome bp -> axis x (bp units, the omitted stretch replaced by the break)"""
        if a1 <= v <= b1:
            return v - a1
        assert a2 <= v <= b2, v
        return (b1 - a1) + gap_bp + (v - a2)
    TOP, ROW0, RH = G["TOP"], G["ROW0"], G["RH"]
    xmax = X(b2)
    ax = bp_axes(fig, XL, TOP, XR - XL, ROW0 + 9 * RH - TOP + 0.5, (0, xmax))
    h, u, nn = ex["h"], ex["u"], ex["nn"]
    (p1s, p1e), (p2s, p2e) = ex["pieces"]
    yh, yt = TOP + 6.0, TOP + 10.4
    eh = 2.0

    def el(s0, e0, y, hh, body, ltr=None, lltr=None, rltr=None, **kw):
        element(ax, X(s0), X(e0), y, hh, body, ltr, lltr=None if lltr is None else X(lltr),
                rltr=None if rltr is None else X(rltr), **kw)
    # truth: nesting diagram (host in two pieces, its tenant) and a non-nested element, real coordinates
    L = int(h["ltr_len"])
    el(p1s, p1e, yh, eh, HOSTC, HOSTC_LTR, lltr=p1s + L)
    el(p2s, p2e, yh, eh, HOSTC, HOSTC_LTR, rltr=p2e - L)
    yb = yh - eh / 2 - 1.1
    ax.plot([X(p1e) - 0.25 * bppm, X(p1e) - 0.25 * bppm, X(p2s) + 0.25 * bppm, X(p2s) + 0.25 * bppm],
            [yh - eh / 2, yb, yb, yh - eh / 2], color=HOSTC, lw=0.6, ls=(0, (2.2, 1.4)), zorder=3)
    ax.text(X(p1s), yb - 0.3, "host", ha="left", va="bottom", fontsize=FS, color=INK)
    chevron(ax, X((p1s + L + p1e) / 2), yh, eh, h["strand"], "white", bppm)
    Lt = int(u["ltr_len"])
    el(u["start"], u["end"], yt, eh, TENC, TENC_LTR, lltr=u["start"] + Lt, rltr=u["end"] - Lt)
    chevron(ax, X((u["start"] + u["end"]) / 2), yt, eh, u["strand"], INK, bppm)
    ax.text(X(u["end"]) + 0.6 * bppm, yt, "tenant", ha="left", va="center", fontsize=FS, color=INK)
    for x in (u["start"], u["end"]):
        ax.plot([X(x), X(x)], [yh + eh / 2, yt - eh / 2], color="#9A9A9A", lw=0.35, zorder=0)
    Ln = int(nn["ltr_len"])
    el(nn["start"], nn["end"], yh, eh, HOSTC, HOSTC_LTR, lltr=nn["start"] + Ln, rltr=nn["end"] - Ln)
    chevron(ax, X((nn["start"] + nn["end"]) / 2), yh, eh, nn["strand"], "white", bppm)
    ax.text(X((nn["start"] + nn["end"]) / 2), yh - eh / 2 - 0.4, "non-nested", ha="center", va="bottom",
            fontsize=FS, color=INK)
    # each tool: what it reports (LTR end blocks, TSD boxes, strand chevron, protein domains); host-level calls above,
    # tenant calls below the row line
    bh = 1.4
    tsd_w = 0.5 * bppm                             # TSDs (4-6 bp) drawn 0.5 mm wide so they can be seen
    for k, t in enumerate(TOOLS):
        y = row_y(k)
        for sa, sb in EX["segs"]:
            ax.plot([X(sa), X(sb)], [y, y], color="#E4E4E4", lw=0.4, zorder=0)
        for f in ex["F"][t]:
            yy = y + 0.9 if f["matched"] == u["id"] else y - 0.9
            body = tint(COL[t], 0.38) if f["lltr"] is not None else tint(COL[t], 0.75)
            el(f["s"], f["e"], yy, bh, body, COL[t], lltr=f["lltr"], rltr=f["rltr"])
            for a, b, _ in f["dom"]:
                ax.add_patch(Rectangle((X(a), yy - bh * 0.28), b - a, bh * 0.56, fc=DOMC, ec="none", zorder=3))
            for a, b in (f["tsd"] or []):
                xa = X(f["s"]) - tsd_w if b <= f["s"] + 1 else X(f["e"])
                ax.add_patch(Rectangle((xa, yy - bh * 0.36), tsd_w, bh * 0.72, fc="white", ec=INK, lw=0.4,
                                       zorder=4))
            if f["strand"]:
                lo = f["lltr"] if f["lltr"] is not None else f["s"]
                hi = f["rltr"] if f["rltr"] is not None else f["e"]
                xc = free_center(lo, hi, [(a, b) for a, b, _ in f["dom"]])
                chevron(ax, X(xc), yy, bh, f["strand"], INK, bppm)
    # broken axis: two segments, "//" at the break, a faint break across the tracks, the omitted length labeled
    yax = ROW0 + 9 * RH + 0.5
    xg0, xg1 = X(b1), X(a2)
    for sa, sb in EX["segs"]:
        ax.plot([X(sa), X(sb)], [yax, yax], color=INK, lw=0.5, clip_on=False)
        for x in np.arange(np.ceil(sa / 5000) * 5000, sb + 1, 5000):
            ax.plot([X(x), X(x)], [yax, yax + 0.8], color=INK, lw=0.5, clip_on=False)
            if abs(X(x) - (xg0 + xg1) / 2) > 7.0 * bppm:          # leave room for the omitted-length label
                ax.text(X(x), yax + 1.2, f"{x / 1e6:.3f}", ha="center", va="top", fontsize=FS)
    for xb in (xg0 + 0.22 * gap_bp, xg1 - 0.22 * gap_bp):
        ax.plot([xb - 0.35 * bppm, xb + 0.35 * bppm], [yax + 0.75, yax - 0.75], color=INK, lw=0.6, clip_on=False)
        ax.plot([xb, xb], [TOP + 1.0, yax - 0.9], color="#C8C8C8", lw=0.4, ls=(0, (1.2, 1.2)), zorder=0)
    ax.text((xg0 + xg1) / 2, yax + 1.2, f"{(a2 - b1) / 1000:.1f} kb", ha="center", va="top", fontsize=FS,
            color=MUTED)
    ax.text(xmax / 2, yax + 4.0, "simulated chromosome 1 (Mb)", ha="center", va="top", fontsize=6.5)
    can.text(G["XN"] - 2.1, yh, "true elements", ha="right", va="center", fontsize=6.5)
    # compact glyph key above the panel
    ky, x = 3.4, XL
    rend = fig.canvas.get_renderer()
    inv = can.transData.inverted()

    def lab(xx, txt):
        tt = can.text(xx, ky, txt, ha="left", va="center", fontsize=FS)
        bb_ = tt.get_window_extent(renderer=rend)
        (a_, _), (b_, _) = inv.transform([(bb_.x0, bb_.y0), (bb_.x1, bb_.y0)])
        return xx + abs(b_ - a_) + 2.4
    can.add_patch(Rectangle((x, ky - 0.7), 3.2, 1.4, fc=tint("#7A7A7A", 0.38), ec="none"))
    can.add_patch(Rectangle((x, ky - 0.7), 0.8, 1.4, fc="#7A7A7A", ec="none"))
    can.add_patch(Rectangle((x + 2.4, ky - 0.7), 0.8, 1.4, fc="#7A7A7A", ec="none"))
    x = lab(x + 3.8, "LTRs")
    can.add_patch(Rectangle((x, ky - 0.5), 0.5, 1.0, fc="white", ec=INK, lw=0.4))
    x = lab(x + 1.1, "TSD")
    can.plot([x, x + 0.6, x], [ky - 0.45, ky, ky + 0.45], color=INK, lw=0.5, solid_joinstyle="miter")
    x = lab(x + 1.2, "strand")
    can.add_patch(Rectangle((x, ky - 0.4), 1.6, 0.8, fc=DOMC, ec="none"))
    x = lab(x + 2.2, "protein domain")
    assert x - 2.4 <= G["XB0"] - 4.5, "A key runs into B"


def nocall_cell(ax, x, y, w, h, pad=0.15, fs=FS):
    """B: the caller makes no call of the class, so precision and F1 are undefined: light neutral cell with '–'"""
    ax.add_patch(Rectangle((x + pad, y + pad), w - 2 * pad, h - 2 * pad, fc=NOCALLC, ec="none", lw=0))
    ax.text(x + w / 2, y + h / 2, "–", ha="center", va="center", fontsize=fs, color=MUTED)


def panel_B(fig, can, Bd):
    """one scorecard on A's nine tool rows: precision, recall, F1 for the nesting classes (non-nested, tenants,
    hosts) and the protein classes (coding, protein-less); system super-header with a rule over its groups, class
    name, n (true elements, 3 genomes summed), group rule, metric labels; viridis 0-1, bold = column best"""
    TOP, ROW0, RH = G["TOP"], G["ROW0"], G["RH"]
    w, gap, sysgap = G["BCW"], G["BGAP"], G["BSYS"]
    x = G["XB0"]
    Bd["drawn_nocall"] = []
    for si, (sysk, syslab, classes) in enumerate(BSYS):
        if si:
            x += sysgap - gap
        sx0 = x
        for ci, (cls, clab) in enumerate(classes):
            gx0 = x
            for m, mlab in CMET:
                can.text(x + w / 2, TOP + 10.6, mlab, ha="center", va="center", fontsize=FS)
                vals = {t: Bd["V"][(cls, m, t)] for t in TOOLS if Bd["V"][(cls, m, t)] is not None}
                best = best_set(vals)
                for k, t in enumerate(TOOLS):
                    if t in vals:
                        score_cell(can, x, ROW0 + k * RH, w, RH, vals[t], t in best)
                    else:
                        nocall_cell(can, x, ROW0 + k * RH, w, RH)
                        Bd["drawn_nocall"].append((cls, m, t))
                x += w
            gw = x - gx0
            can.text(gx0 + gw / 2, TOP + 1.4, clab, ha="center", va="center", fontsize=6.5)
            can.text(gx0 + gw / 2, TOP + 4.3, f"n = {Bd['n'][cls]:,}", ha="center", va="center", fontsize=FS,
                     color=MUTED)
            can.plot([gx0 + 0.4, x - 0.4], [TOP + 6.6, TOP + 6.6], color=INK, lw=0.4)
            x += gap
        sx1 = x - gap
        can.text((sx0 + sx1) / 2, TOP - 4.6, syslab, ha="center", va="center", fontsize=6.5, color=MUTED)
        can.plot([sx0 + 0.4, sx1 - 0.4], [TOP - 3.0, TOP - 3.0], color=INK, lw=0.4)
    assert abs(x - gap - G["XB1"]) < 1e-6, x
    assert Bd["drawn_nocall"] == [("proteinless", "precision", "dante_ltr"), ("proteinless", "f1", "dante_ltr")], \
        Bd["drawn_nocall"]


def recall_plot(fig, x0, x1, y0, h, S, ylab):
    """C and D: pooled recall vs median K2P, one line + 95% band per tool, one dot at the pooled recall-halved point,
    documented identity floors dashed; same axes, sizes, line weights and styling for both"""
    ax = poafig.mm_axes(fig, x0, y0, x1 - x0, h)
    for fl, lab, ha, dx in ((85, "85% identity", "left", 0.005), (70, "70%", "right", -0.005)):
        xf = k2p_of_identity(fl)
        ax.axvline(xf, color=INK, lw=0.5, ls=(0, (2.2, 1.6)), zorder=1)
        ax.text(xf + dx, 1.02, lab, ha=ha, va="top", fontsize=FS)
    for t in TOOLS[::-1]:
        v = S[t]
        ax.fill_between(v["x"], v["lo"], v["hi"], color=COL[t], alpha=BAND_A, lw=0, zorder=2)
        ax.plot(v["x"], v["r"], color=COL[t], lw=LW_Q if t == "ltrquest" else LW_O, zorder=3 + (t == "ltrquest"),
                solid_capstyle="round")
        if v["d50"] not in ("NA", "incomplete") and not v["d50"].startswith(">"):
            ax.plot(float(v["d50"]), v["r0"] / 2, "o", ms=2.6, mfc=COL[t], mec="white", mew=0.4, zorder=7)
    ax.set_xlim(0, 0.40)
    ax.set_ylim(-0.03, 1.04)                 # lines at 0 stay visible above the axis
    ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4])
    ax.set_xticklabels(["0", "0.1", "0.2", "0.3", "0.4"])
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_yticklabels(["0", "0.25", "0.5", "0.75", "1"])
    clean(ax)
    ax.set_xlabel("median LTR divergence (K2P)", fontsize=6.5, labelpad=1.5)
    ax.set_ylabel(ylab, fontsize=6.5, labelpad=1.5)
    return ax


def panel_C(fig, can, Cd):
    return recall_plot(fig, G["XC0"], G["XC1"], G["Y2"] + 3.0, G["HD"], Cd["S"], "coding elements found (recall)")


def key_CD(fig, can):
    """one tool key for C and D, between them"""
    x0, y = G["XK0"], G["Y2"] + 5.0
    dy = 2.7
    for i, t in enumerate(TOOLS):
        yy = y + i * dy
        can.plot([x0, x0 + 3.4], [yy, yy], color=COL[t], lw=LW_Q if t == "ltrquest" else LW_O,
                 solid_capstyle="butt")
        can.text(x0 + 4.1, yy, NAME[t], ha="left", va="center", fontsize=FS,
                 fontweight="bold" if t == "ltrquest" else "normal")
    yy = y + 9 * dy + 1.2
    can.add_patch(Rectangle((x0, yy - 0.7), 3.4, 1.4, fc="#8C8C8C", alpha=BAND_A * 1.6, ec="none"))
    can.text(x0 + 4.1, yy, "95% CI", ha="left", va="center", fontsize=FS)
    yy += dy
    can.plot([x0 + 1.7], [yy], "o", ms=2.6, mfc=INK, mec="white", mew=0.4)
    can.text(x0 + 4.1, yy, "recall halved", ha="left", va="center", fontsize=FS)


def panel_D(fig, can, Dd):
    return recall_plot(fig, G["XD0"], G["XD1"], G["Y2"] + 3.0, G["HD"], Dd["S"], "protein-less found (recall)")


def panel_E_F(fig, can, Fd):
    """band 3 on shared tool rows: E masking precision / recall / F1 for three genomes (blocks of three columns, laid
    out like B; a missing result spans its genome's block), F family share for five genomes; the one color bar"""
    y0 = G["Y3"]
    RH = G["RH3"]
    yR = y0 + 7.4
    ys = [yR + (k + 0.5) * RH for k in range(len(TOOLS))]
    names(can, G["XN"], ys)
    # the one score scale of B, E and F (viridis, 0-1)
    xb, yb, bw = 6.4, y0 + 5.2, 10.0
    for i_ in range(50):
        can.add_patch(Rectangle((xb + i_ * bw / 50, yb - 0.65), bw / 50 + 0.02, 1.3, fc=CMAP((i_ + 0.5) / 50),
                                ec="none"))
    can.add_patch(Rectangle((xb, yb - 0.65), bw, 1.3, fc="none", ec=INK, lw=0.3))
    can.text(xb - 0.6, yb, "0", ha="right", va="center", fontsize=FS)
    can.text(xb + bw + 0.6, yb, "1", ha="left", va="center", fontsize=FS)
    can.text(xb + bw / 2, yb - 2.3, "score (B, E, F)", ha="center", va="center", fontsize=FS)
    drawn = []
    # E
    gap = G["EGAP"]
    cw = (G["XE1"] - G["XE0"] - (len(E_DS) - 1) * gap) / (len(E_DS) * len(EMET))
    for j, ds in enumerate(E_DS):
        gx = G["XE0"] + j * (len(EMET) * cw + gap)
        gw = len(EMET) * cw
        can.text(gx + gw / 2, y0 + 1.6, DSNAME[ds], ha="center", va="center", fontsize=6.5)
        can.plot([gx + 0.4, gx + gw - 0.4], [y0 + 3.4, y0 + 3.4], color=INK, lw=0.4)
        for i, (_, key, mlab) in enumerate(EMET):
            x = gx + i * cw
            can.text(x + cw / 2, y0 + 5.4, mlab, ha="center", va="center", fontsize=FS)
            col = {t: Fd["mk"][(ds, t)][key] for t in TOOLS if (ds, t) in Fd["mk"]}
            best = best_set(col)
            for k, t in enumerate(TOOLS):
                if t in col:
                    score_cell(can, x, yR + k * RH, cw, RH, col[t], t in best)
        for k, t in enumerate(TOOLS):
            if (ds, t) not in Fd["mk"]:
                kind = Fd["reasons"][(ds, t)][0]
                missing_cell(can, gx, yR + k * RH, gw, RH, kind)
                drawn.append(("E", ds, t, kind))
    # F
    cw = (G["XF1"] - G["XF0"]) / len(DS)
    for j, ds in enumerate(DS):
        x = G["XF0"] + j * cw
        can.text(x + cw / 2, y0 + 1.6, DSNAME[ds], ha="center", va="center", fontsize=6.5)
        can.text(x + cw / 2, y0 + 5.4, f"n = {Fd['ncur'][ds]}", ha="center", va="center", fontsize=FS, color=MUTED)
        col = {t: Fd["pg"][(ds, t)] for t in TOOLS if (ds, t) in Fd["pg"]}
        best = best_set(col)
        for k, t in enumerate(TOOLS):
            if t in col:
                score_cell(can, x, yR + k * RH, cw, RH, col[t], t in best)
            else:
                kind = Fd["reasons"][(ds, t)][0]
                missing_cell(can, x, yR + k * RH, cw, RH, kind)
                drawn.append(("F", ds, t, kind))
    for panel, exp in (("E", E_MISSING), ("F", F_MISSING)):
        got = {k: {(ds, t) for p_, ds, t, kk in drawn if p_ == panel and kk == k} for k in ("failed", "noteval")}
        assert got == exp, (panel, got)
    return yR + 9 * RH


def panel_G(fig, can, Gd):
    """A. thaliana loci (same rendering as figS3_library_age; the heatmap height
    ends level with H's last bar row)"""
    y0 = G["Y4"]
    x0, BW = G["XG0"], G["GBW"]
    gap = 1.4
    cmap = LinearSegmentedColormap.from_list("cov", ["#FFFFFF", "#262626"])
    k = Gd["k2p"]
    ok = ~np.isnan(k)
    bins = np.digitize(k[ok], [0.01, 0.03, 0.1])
    order = np.lexsort((k[ok],))
    P, Q = Gd["prof"][ok][order], Gd["nprof"][ok][order]
    bb = bins[order]
    n_tot = int(ok.sum())
    can.text(x0 + BW + gap / 2, y0 + 1.0, r"$\it{A.\ thaliana}$" + f", n = {n_tot:,}", ha="center", va="center",
             fontsize=6.5)
    num("G", "A. thaliana protein-confirmed loci drawn (K2P measured)", n_tot,
        "report/tables/library_age_bias.tsv n_loci (four bins); rows rebuilt from eval/adj")
    yP, hP = y0 + 7.6, 7.0
    yHt, hHt = yP + hP + 1.3, G["GHT"]
    for j, (M, title) in enumerate(((P, "loci"), (Q, "moved"))):
        bx = x0 + j * (BW + gap)
        can.text(bx + BW / 2, y0 + 3.6, title, ha="center", va="center", fontsize=FS, color=INK if j == 0 else MUTED)
        gy = y0 + 5.9
        sc = BW / NBIN
        can.add_patch(Rectangle((bx + NB5 * sc, gy - 0.35), NBI * sc, 0.7, fc=INTC, ec="none"))
        can.add_patch(Rectangle((bx, gy - 0.9), NB5 * sc, 1.8, fc=LTR5C, ec="none"))
        can.add_patch(Rectangle((bx + (NB5 + NBI) * sc, gy - 0.9), NB3 * sc, 1.8, fc=LTR3C, ec="none"))
        axp = poafig.mm_axes(fig, bx, yP, BW, hP)
        xs = np.arange(NBIN) + 0.5
        if j == 0:
            for b_, shade in ((0, 0.0), (3, 0.55)):
                prof_ = np.nanmean(M[bb == b_], axis=0)
                axp.plot(xs, prof_, color=tint(LOCIC, 1 - shade), lw=0.9)
                num("G", f"mean curated coverage along loci, K2P {BINS[b_]} (mean over bins)",
                    round(float(np.nanmean(prof_)), 3), "computed from eval/adj + mask/arab/curated (tmp cache)")
        else:
            prof_ = np.nanmean(M, axis=0)
            axp.plot(xs, prof_, color=NULLC, lw=0.9)
            num("G", "mean curated coverage along moved loci (mean over bins)", round(float(np.nanmean(prof_)), 3),
                "computed from eval/adj null_shift + mask/arab/curated (tmp cache)")
        axp.set_xlim(0, NBIN)
        axp.set_ylim(0, 1.05)
        clean(axp)
        axp.set_xticks([])
        axp.set_yticks([0, 1])
        if j == 0:
            axp.set_yticklabels(["0", "1"])
        else:
            axp.set_yticklabels([])
        axh = poafig.mm_axes(fig, bx, yHt, BW, hHt)
        for b_ in range(4):
            img = np.nan_to_num(M[bb == b_], nan=0.0)
            if img.shape[0]:
                axh.imshow(img, cmap=cmap, vmin=0, vmax=1, aspect="auto", interpolation="antialiased",
                           extent=(0, NBIN, b_ + 1, b_), origin="upper")
        for b_ in (1, 2, 3):
            axh.axhline(b_, color="white", lw=0.8)
        axh.set_xlim(0, NBIN)
        axh.set_ylim(4, 0)
        axh.set_xticks([])
        for s_ in axh.spines.values():
            s_.set_linewidth(0.4)
        axh.set_yticks([1, 2, 3])
        axh.tick_params(axis="y", length=1.8, pad=1.2, labelsize=FS)
        axh.set_yticklabels(["0.01", "0.03", "0.1"] if j == 0 else [])
    for b_ in range(4):
        num("G", f"loci in K2P class {BINS[b_]} (drawn)", int((bb == b_).sum()), "report/tables/library_age_bias.tsv")
    can.text(x0 - 6.4, yP + hP / 2, "mean\ncoverage", ha="center", va="center", rotation=90, fontsize=FS,
             linespacing=0.95)
    can.text(x0 - 8.6, yHt + hHt / 2, "LTR divergence (K2P)", ha="center", va="center", rotation=90, fontsize=6.5)
    # key: profile lines, then the shading scale
    yk = yHt + hHt + 2.4
    items = [(LOCIC, "K2P < 0.01"), (tint(LOCIC, 0.45), "K2P ≥ 0.1"), (NULLC, "moved")]
    xk = x0 - 9.0
    for i, (c, lab) in enumerate(items):
        xx = xk + (0 if i != 1 else 18.5)
        yy = yk + (0 if i < 2 else 2.6)
        can.plot([xx, xx + 3.0], [yy, yy], color=c, lw=1.0, solid_capstyle="butt")
        can.text(xx + 3.7, yy, lab, ha="left", va="center", fontsize=FS)
    xb, yb, bw = xk + 18.5 + 7.6, yk + 2.6, 9.0
    for i_ in range(40):
        can.add_patch(Rectangle((xb + i_ * bw / 40, yb - 0.65), bw / 40 + 0.02, 1.3, fc=cmap((i_ + 0.5) / 40),
                                ec="none"))
    can.add_patch(Rectangle((xb, yb - 0.65), bw, 1.3, fc="none", ec=INK, lw=0.3))
    can.text(xb - 0.6, yb, "0", ha="right", va="center", fontsize=FS)
    can.text(xb + bw + 0.6, yb, "100%", ha="left", va="center", fontsize=FS)
    can.text(xb + bw / 2, yb + 2.4, "curated annotation", ha="center", va="center", fontsize=FS)
    return dict(bb=bb, P=P, Q=Q)


def h_columns():
    """H: x position and width of each bar column (all loci, then the classes), and the right end"""
    xs, ws, x, prev = [], [], G["XHB0"], None
    for sysk, _, _, _ in HCOLS:
        if prev is not None:
            x += G["HGAP"] if sysk == prev else G["HSYS"]
        w = G["HALLW"] if sysk == "all" else G["HCLW"]
        xs.append(x)
        ws.append(w)
        x += w
        prev = sysk
    return xs, ws, x


def h_rows():
    """H: y of the header block and of the first bar row"""
    yH = G["Y4"] + 6.0
    return yH, yH + 6.2


def hatch_rect(can, x, y, w, h, fc, lc="white", lw=HATCH_LW_H, step=HATCH_STEP_H):
    """rectangle (x, top y, w, h in page mm, y down) in fc (None: nothing drawn, only the clip for lines laid over
    what is already there) with 45-degree lines in lc clipped to it: explicit lines, not a PDF tiling pattern
    (poppler renders those faint), so the PNG and every PDF viewer show the same marks; the lines sit on one
    page-wide grid, so neighbouring slices continue the same pattern"""
    r = Rectangle((x, y), w, h, fc="none" if fc is None else fc, ec="none", lw=0)
    can.add_patch(r)
    if w <= 0:
        return
    s0 = math.floor((x - h) / step) * step
    segs = [[(s, y + h), (s + h, y)] for s in np.arange(s0, x + w + step, step)]
    lc_ = LineCollection(segs, colors=lc, linewidths=lw, capstyle="butt")
    can.add_collection(lc_)
    lc_.set_clip_path(r)


def fill_bar(can, x, yc, w, h, n, solid, alone, t):
    """one 0-100% bar: light gray ground; genotyped correctly in the caller's color, first the loci another caller
    also nominated (solid), then the loci this caller alone nominated (white hatch over its color); true widths"""
    can.add_patch(Rectangle((x, yc - h / 2), w, h, fc=RESTC, ec="none", lw=0))
    ws, wa = w * solid / n, w * alone / n
    if solid:
        can.add_patch(Rectangle((x, yc - h / 2), ws, h, fc=COL[t], ec="none", lw=0))
    if alone:
        hatch_rect(can, x + ws, yc - h / 2, wa, h, COL[t])


def panel_H(fig, can, Hd):
    """ltrkin: per caller, haplotype 2 above haplotype 1: one 0-100% bar for all loci (share and 95% CI printed at
    the right) and one per class (non-nested, tenants, hosts | coding, protein-less; no intervals); colored = genotyped
    correctly (hatched: loci this caller alone nominated), light gray = the rest; tool failed = dark gray"""
    y0 = G["Y4"]
    RH = G["RH4"]
    xs, ws, X1 = h_columns()
    yH, yR = h_rows()
    ys = [yR + (k + 0.5) * RH for k in range(len(TOOLS))]
    names(can, G["XHN"], ys)
    hb, hg = 1.75, 0.45
    can.text(G["XHN"] + 2.1, yR - 1.2, "hap.", ha="center", va="center", fontsize=FS, color=MUTED)
    # headers: class name, n "hap. 2 / hap. 1"; system super-headers with a rule over their classes
    for (_, _, _, lab), x, w in zip(HCOLS, xs, ws):
        can.text(x + w / 2, yH + 2.2, lab, ha="center", va="center", fontsize=FS)
        can.text(x + w / 2, yH + 4.4, f"{Hd['n'][('camara_hap2', lab)]} / {Hd['n'][('camara_hap1', lab)]}",
                 ha="center", va="center", fontsize=FS, color=MUTED)
    for sysk, slab in (("nesting", "nesting"), ("protein", "protein domains")):
        idx = [i for i, c in enumerate(HCOLS) if c[0] == sysk]
        xa, xb = xs[idx[0]], xs[idx[-1]] + ws[idx[-1]]
        can.plot([xa + 0.2, xb - 0.2], [yH + 0.6, yH + 0.6], color=INK, lw=0.4)
        can.text((xa + xb) / 2, yH - 0.4, slab, ha="center", va="bottom", fontsize=FS, color=MUTED)
    can.text(X1 + 1.2, yR - 1.2, "all: correct (95% CI)", ha="left", va="center", fontsize=FS, color=MUTED)
    for k, t in enumerate(TOOLS):
        for j, ds in enumerate(HDS):
            yc = ys[k] + (-1 if j == 0 else 1) * (hb + hg) / 2
            can.text(G["XHN"] + 2.1, yc, HAPN[ds], ha="center", va="center", fontsize=FS, color=MUTED)
            o = Hd["o"][(ds, t)]
            if o["status"] != "ok":                  # ltrkin_yield.tsv: the caller did not finish
                for x, w in zip(xs, ws):
                    can.add_patch(Rectangle((x, yc - hb / 2), w, hb, fc=FAILC, ec="none", lw=0))
                # baseline placement: 6-pt ascenders (~1.54 mm) fit inside the 1.75-mm bar (va="center" let them
                # rise above it, white on white)
                can.text(xs[0] + ws[0] / 2, yc + 0.77, "tool failed", ha="center", va="baseline", fontsize=FS,
                         color="white")
                continue
            for (_, _, _, lab), x, w in zip(HCOLS, xs, ws):
                c = o["cls"][lab]
                fill_bar(can, x, yc, w, hb, c["n"], c["solid"], c["alone"], t)
            can.text(X1 + 1.2, yc, f"{pct(o['share'])}% ({pct(o['lo'])}–{pct(o['hi'])})", ha="left", va="center",
                     fontsize=FS, fontweight="bold" if t == "ltrquest" else "normal")
    yax = yR + 9 * RH + 0.4
    for i, (x, w) in enumerate(zip(xs, ws)):
        can.plot([x, x + w], [yax, yax], color=INK, lw=0.5)
        for v in ((0, 50, 100) if i == 0 else (0, 100)):
            can.plot([x + w * v / 100] * 2, [yax, yax + 0.8], color=INK, lw=0.5)
            if i == 0:
                can.text(x + w * v / 100, yax + 1.2, f"{v}", ha="center", va="top", fontsize=FS)
    can.text((xs[0] + X1) / 2, yax + 4.0, "share genotyped correctly (%)", ha="center", va="top", fontsize=6.5)
    # key, one row (as Fig. S4's): genotyped correctly (three caller colors as a sample strip), nominated by this
    # caller alone (the same strip under the white hatch, so it does not read as a third gray category), not
    # genotyped correctly, tool failed
    ky = y0 + 1.0
    x = G["XH"] + 4.5
    items = [("strip", "genotyped correctly (caller's color)"),
             ("strip_hatch", "nominated by this caller alone"),
             (dict(fc=RESTC), "not genotyped correctly"),
             (dict(fc=FAILC), "tool failed")]
    rend = fig.canvas.get_renderer()
    inv = can.transData.inverted()
    for st, lab in items:
        if st in ("strip", "strip_hatch"):
            for i_, t in enumerate(("ltrquest", "ltr_retriever", "dante_ltr")):
                can.add_patch(Rectangle((x + i_ * 1.0, ky - 0.75), 1.0, 1.5, fc=COL[t], ec="none", lw=0))
            if st == "strip_hatch":
                hatch_rect(can, x, ky - 0.75, 3.0, 1.5, None)
        else:
            can.add_patch(Rectangle((x, ky - 0.75), 3.0, 1.5, ec="none", lw=0, **st))
        tt = can.text(x + 3.7, ky, lab, fontsize=FS, ha="left", va="center")
        bb_ = tt.get_window_extent(renderer=rend)
        (a_, _), (b_, _) = inv.transform([(bb_.x0, bb_.y0), (bb_.x1, bb_.y0)])
        x += 3.7 + abs(b_ - a_) + 3.0
    log(f"H key ends at {x - 3.0:.1f} mm; label column from {X1 + 1.2:.1f} mm")
    assert x - 3.0 <= W - 0.5, f"H key overflows the page ({x - 3.0:.1f} mm)"
    return yax


# ============================================================================ legend
def legend_text(ex, Bd, Cd, Dd, Fd, Gd, Hd):
    found_o = [NAME[t] for t in ex["found_by"] if t != "ltrquest"]
    miss_t = [NAME[t] for t in TOOLS if t not in ex["tfound_by"]]
    assert len(found_o) == 3 and len(miss_t) == 1 and len(ex["nfound_by"]) == len(TOOLS)
    # B claims
    V = Bd["V"]
    val = lambda c, m, ts: [V[(c, m, t)] for t in ts if V[(c, m, t)] is not None]     # noqa: E731
    for _, c, _ in BCLS:
        for m in ("recall", "f1"):
            assert V[(c, m, "ltrquest")] > max(val(c, m, TOOLS[1:])), (c, m)
        assert V[(c, "precision", "ltrquest")] < max(val(c, "precision", TOOLS[1:])), c
        assert V[(c, "precision", "hite")] >= max(val(c, "precision", TOOLS)), c
    rh_q, rp_q = V[("host", "recall", "ltrquest")], V[("proteinless", "recall", "ltrquest")]
    rh_o, rp_o = val("host", "recall", TOOLS[1:]), val("proteinless", "recall", TOOLS[1:])
    nocall = sorted({(c, t) for (c, m, t), v in V.items() if v is None})
    assert nocall == [("proteinless", "dante_ltr")], nocall
    sd = Bd["sd_max"]
    sd_o = max(v for t, v in sd.items() if t != "ltrpred")
    assert sd_o <= 0.023 and sd["ltrpred"] == max(sd.values()) and sd["ltrpred"] > 0.05, sd
    num("memo", "B: LTRquest host / protein-less recall; other tools' ranges (printed in B's cells; legend says only "
        "that LTRquest leads on recall and F1 in every class)",
        f"{rh_q:.4f} / {rp_q:.4f}; {min(rh_o):.4f}-{max(rh_o):.4f} / {min(rp_o):.4f}-{max(rp_o):.4f}",
        "experiments/figpub/sim_class_prf_joint_summary.tsv set=S2 recall_mean")
    num("legend", "B: largest across-genome SD, every tool but LTRpred / LTRpred (drawn classes, P/R/F1)",
        f"{sd_o:.4f} / {sd['ltrpred']:.4f}", "experiments/figpub/sim_class_prf_joint_summary.tsv *_sd")
    num("legend", "B: LTRpred's S2_3 run incomplete (no chr6 calls, chr5 partly; gt ltrdigest HMMER error, recorded ok)",
        "one run of 60", "experiments/figpub/memo.md section e (runs/S2_3/ltrpred/run.log)")
    # C claims (coding, pooled recall-halved points)
    CS = Cd["S"]
    d50c = {t: float(CS[t]["d50"]) for t in TOOLS if CS[t]["d50"][0].isdigit()}
    seven = [t for t in TOOLS if t not in ("ltrquest", "dante_ltr")]
    lo7, hi7 = min(d50c[t] for t in seven), max(d50c[t] for t in seven)
    q50 = d50c["ltrquest"]
    assert q50 > 0.3 and hi7 < 0.2 and CS["dante_ltr"]["d50"].startswith(">"), (q50, hi7)
    for t in seven:
        if FLOOR[t] is not None:
            assert d50c[t] < k2p_of_identity(FLOOR[t]), t
    dq0 = CS["dante_ltr"]["r0"]
    num("legend", "C: coding pooled d50 range of the seven other tools (K2P)", f"{lo7:.4f}-{hi7:.4f}", "computed (C dots)")
    num("memo", "C: DANTE_LTR coding pooled recall at timepoint 0 (legend: never halves)", f"{dq0:.4f}",
        "computed (C)")
    # D claims (protein-less)
    DS_ = Dd["S"]
    eq, eh = float(DS_["ltrquest"]["d50"]), float(DS_["hirepeat"]["d50"])
    dh = d50c["hirepeat"]
    assert eq < q50 and eh > dh
    assert DS_["dante_ltr"]["r"].max() == 0
    lp0 = DS_["ltrpred"]["r0"]
    assert lp0 < 0.25
    num("legend", "C/D: pooled d50 LTRquest coding / protein-less; HiRepEat coding / protein-less",
        f"{q50:.4f} / {eq:.4f}; {dh:.4f} / {eh:.4f}", "computed (C and D dots)")
    num("legend", "D: LTRpred protein-less pooled recall at timepoint 0", f"{lp0:.4f}", "computed (D)")
    # E claims
    mk = Fd["mk"]
    rec = [v["recall"] for v in mk.values()]
    prc = [v["precision"] for v in mk.values()]
    aq, ar, ah = (mk[("arab", t)]["f1"] for t in ("ltrquest", "ltr_retriever", "hite"))
    assert aq < min(ar, ah)
    for ds in E_DS:                                 # "LTRquest trails the strictest in F1", every E genome
        assert mk[(ds, "ltrquest")]["f1"] < max(v["f1"] for (d, t), v in mk.items() if d == ds), ds
    sq = [mk[(ds, "ltrquest")]["recall"] for ds in E_DS]
    pq = [mk[(ds, "ltrquest")]["precision"] for ds in E_DS]
    assert min(sq) >= 0.985 and min(rec) > max(prc) - 0.1 and min(rec) >= 0.75, (sq, min(rec))
    num("legend", "E: recall range / precision range over the drawn cells (3 genomes)",
        f"{min(rec):.4f}-{max(rec):.4f} / {min(prc):.4f}-{max(prc):.4f}", "report/tables/masking.tsv sens, prec")
    num("legend", "E: LTRquest recall / precision (A. thaliana, rice, zebra finch)",
        " / ".join(f"{a:.4f}" for a in sq) + "; " + " / ".join(f"{a:.4f}" for a in pq), "report/tables/masking.tsv")
    num("legend", "E: 'LTRquest trails the strictest in F1', A. thaliana F1 LTRquest / LTR_retriever / HiTE (also in "
        "E's cells)", f"{aq:.4f} / {ar:.4f} / {ah:.4f}", "report/tables/masking.tsv F1")
    # F claims
    pg = Fd["pg"]
    lead_pg = [ds for ds in DS if max(v for (d, t), v in pg.items() if d == ds) == pg.get((ds, "ltrquest"), -1)]
    assert lead_pg == ["arab", "rice", "maize", "gorilla"], lead_pg
    fin = {t: v for (d, t), v in pg.items() if d == "finch"}
    assert max(fin, key=fin.get) == "look4ltrs"
    assert max(v for (d, t), v in pg.items() if d == "gorilla") < 0.1
    # G claims
    ht = Gd["tab"]
    s_y, s_o = float(ht["<0.01"]["supported_share"]), float(ht[">=0.1"]["supported_share"])
    fpq = float(Fd["m"][(Fd["m"].ds == "arab") & (Fd["m"].tool == "ltrquest")].FP_mb.iloc[0])
    frac = Gd["gap_mb"] / fpq
    assert frac <= 0.205
    num("legend", "G: A. thaliana library-supported share K2P < 0.01 / >= 0.1", f"{s_y} / {s_o}",
        "report/tables/library_age_bias.tsv supported_share")
    num("legend", "G: protein-confirmed bases outside curated annotation / LTRquest masking false-positive bases",
        f"{Gd['gap_mb']} / {fpq} Mb = {frac:.3f}", "annotation_gap.tsv gap_mb / masking.tsv FP_mb (arab)")
    # H claims
    o = Hd["o"]
    okt = {ds: [t for t in TOOLS if o[(ds, t)]["status"] == "ok"] for ds in HDS}
    for ds in HDS:
        srt = sorted(okt[ds], key=lambda t: -o[(ds, t)]["share"])
        assert srt[:2] == ["ltrquest", "hirepeat"], (ds, srt)
        srt = sorted(okt[ds], key=lambda t: -o[(ds, t)]["core_share"])
        assert srt[:2] == ["ltrquest", "hirepeat"], (ds, srt)
    assert o[("camara_hap2", "look4ltrs")]["status"].startswith("DNF") and okt["camara_hap1"] == TOOLS

    def csh(ds, t, lab):
        c = o[(ds, t)]["cls"][lab]
        return (c["solid"] + c["alone"]) / c["n"]
    hq = [csh(ds, "ltrquest", "hosts") for ds in HDS]
    pq_ = [csh(ds, "ltrquest", "protein-less") for ds in HDS]
    ho = [csh(ds, t, "hosts") for ds in HDS for t in okt[ds] if t != "ltrquest"]
    po = [csh(ds, t, "protein-less") for ds in HDS for t in okt[ds] if t != "ltrquest"]
    assert min(hq) > max(ho) and min(pq_) > max(po) and min(ho) == 0 and min(po) == 0
    # "partly on loci only its calls nominated": LTRquest's hatched part in hosts and protein-less loci, both haplotypes
    for lab in ("hosts", "protein-less"):
        c = [o[(ds, "ltrquest")]["cls"][lab] for ds in HDS]
        assert all(x["alone"] > 0 for x in c), lab
        oth = [csh(ds, t, lab) for ds in HDS for t in okt[ds] if t != "ltrquest"]
        num("memo", f"H: LTRquest {lab} correct share without the loci only its calls nominated (hap2, hap1) / "
            f"best other caller", f"{c[0]['solid'] / c[0]['n']:.4f}, {c[1]['solid'] / c[1]['n']:.4f} / "
            + ", ".join(f"{max(csh(ds, t, lab) for t in okt[ds] if t != 'ltrquest'):.4f}" for ds in HDS),
            "experiments/figpub/ltrkin_locus_class.tsv (solid part of the bar)")
        assert min(oth) == 0
    alone = {(ds, t): (o[(ds, t)]["n_alone"], o[(ds, t)]["alone_correct"]) for ds in HDS for t in okt[ds]}
    exp_alone = {("camara_hap2", "ltrquest"): (39, 38), ("camara_hap1", "ltrquest"): (26, 23),
                 ("camara_hap2", "hirepeat"): (5, 3), ("camara_hap1", "hirepeat"): (4, 1),
                 ("camara_hap2", "ltrharvest"): (2, 1), ("camara_hap2", "ltr_finder"): (1, 1),
                 ("camara_hap2", "hite"): (1, 1), ("camara_hap1", "look4ltrs"): (1, 0)}
    assert {k: v for k, v in alone.items() if v[0]} == exp_alone, {k: v for k, v in alone.items() if v[0]}
    assert alone[("camara_hap2", "ltr_finder")][0] == alone[("camara_hap2", "hite")][0]   # "1 each"
    num("legend", "H: loci nominated by one caller alone, other callers (hatching too thin to see), hap2 / hap1",
        "; ".join(f"{NAME[t]} {alone.get(('camara_hap2', t), ('failed',))[0]} / {alone[('camara_hap1', t)][0]}"
                  for t in TOOLS[1:] if any(alone.get((ds, t), (0,))[0] for ds in HDS)),
        "experiments/figpub/ltrkin_locus_class.tsv nominators")
    cq = [o[(ds, "ltrquest")]["core_share"] for ds in HDS]
    ch = [o[(ds, "hirepeat")]["core_share"] for ds in HDS]
    lead = [100 * (a - b) for a, b in zip(cq, ch)]
    num("legend", "H: LTRquest correct share of hosts (hap2, hap1) / other finished callers' range",
        f"{hq[0]:.4f}, {hq[1]:.4f} / {min(ho):.4f}-{max(ho):.4f}", "experiments/figpub/ltrkin_locus_class.tsv")
    num("legend", "H: LTRquest correct share of protein-less loci (hap2, hap1) / other finished callers' range",
        f"{pq_[0]:.4f}, {pq_[1]:.4f} / {min(po):.4f}-{max(po):.4f}", "experiments/figpub/ltrkin_locus_class.tsv")
    num("legend", "H: on loci two or more callers nominated, LTRquest / HiRepEat correct share (hap2, hap1; shares "
        "memo only); lead (points, legend)",
        f"{cq[0]:.4f} / {ch[0]:.4f}, {cq[1]:.4f} / {ch[1]:.4f}; {lead[0]:.1f}, {lead[1]:.1f}",
        "experiments/figpub/ltrkin_locus_class.tsv nominators, out_<tool>")

    txt = (
        f"**Fig. 1. Where LTR retrotransposon callers lose elements, and what recovering them costs.** Callers lose "
        f"elements split by a later insertion, old copies and protein-less elements, which carry no LTR "
        f"retrotransposon protein domain and so cannot transpose on their own. LTRquest recovers most and genotypes "
        f"the most real polymorphisms but masks much sequence curated libraries leave unannotated. Panels A–D use "
        f"genomes simulated with PrinTE, which shares an author with LTRquest, and LTRquest was tuned on PrinTE "
        f"simulations and on rice.\n\n"
        f"**A**, A host, its tenant and a non-nested element in a simulated genome (broken axis). A host is an element "
        f"split in two by a later insertion, its tenant, and may itself be a tenant. A non-nested element is neither. Each "
        f"tool's calls follow to scale with the features it reports (TSDs enlarged). LTRquest, {found_o[0]}, "
        f"{found_o[1]} and {found_o[2]} rebuild the host, all but {miss_t[0]} the tenant, and every tool the "
        f"non-nested element.\n\n"
        f"**B**, Precision, recall and F1 by class, averaged over three deeply nested simulated genomes. Bold marks "
        f"each column's best, n the true elements, and a dash a class the caller never calls. Coding elements carry "
        f"at least two kinds of LTR retrotransposon domain from one superfamily (partly coding ones in Figs. S1 and "
        f"S4). Recall counts true elements by their family's class, precision counts calls by their own domains. "
        f"LTRquest leads on recall and F1 in every class, HiTE on precision. The genomes agree (SD at most "
        f"{sd_o:.2f}) except for LTRpred, one of whose runs lost calls to a pipeline error. Precision runs higher in "
        f"simulations than on real genomes, which carry tandem duplicated genes and incomplete curated "
        f"libraries.\n\n"
        f"**C**, Recall of coding elements as one simulated burst's LTRs diverge, three genomes pooled (bands, 95% "
        f"intervals over families). Seven tools lose half (dots) by a divergence of {lo7:.2f}–{hi7:.2f}, those with a "
        f"documented identity floor (dashed) well before it. LTRquest keeps half to {q50:.2f}, and DANTE_LTR, which "
        f"reads proteins, never halves.\n\n"
        f"**D**, Protein-less elements of the same series. LTRquest loses half by {eq:.2f}, while HiRepEat keeps them "
        f"longer than its coding elements ({eh:.2f} against {dh:.2f}). DANTE_LTR finds none and LTRpred few. Details "
        f"in Fig. S1.\n\n"
        f"**E**, Masking of three real genomes by each tool's library, against curated annotation. Recall is high "
        f"for every tool ({min(rec):.2f}–{max(rec):.2f}), so tools differ in precision, and LTRquest trails the "
        f"strictest in F1. Maize and gorilla are in Fig. S2C. Dark gray marks tools that failed (no result within 36 h or "
        f"200 GB of memory, or a crash without calls), and hatching marks libraries whose masking or building would have run far past 36 h or was "
        f"stopped there.\n\n"
        f"**F**, Share of curated LTR families each library recovers perfectly or well (n, curated families). "
        f"LTRquest leads in four genomes and Look4LTRs in zebra finch, though every finished library recovers few "
        f"gorilla families. Missing cells as in E.\n\n"
        f"**G**, Every protein-confirmed *A. thaliana* locus any tool called, sorted by LTR divergence and shaded by "
        f"curated annotation, beside the same spans moved 20–200 kb (chance). Curated libraries support fewer old "
        f"loci ({s_o:.2f} against {s_y:.2f} for young), but such sequence is at most {pct(frac)}% of LTRquest's "
        f"false-positive masking. Fig. S3 shows four more genomes.\n\n"
        f"**H**, Real insertion polymorphisms of a *Cardamine amara* plant's two haplotypes, genotyped from each "
        f"caller's annotation (haplotype 2 first). Color marks the share genotyped correctly, with a 95% interval "
        f"for all loci (per class in Fig. S4). True loci pool all callers' nominations, each confirmed by an empty "
        f"site on the other haplotype, and hatching marks those only that caller nominated. LTRquest genotypes "
        f"{hq[0]:.2f} and {hq[1]:.2f} of hosts and {pq_[0]:.2f} and {pq_[1]:.2f} of protein-less loci, partly its "
        f"own, against {min(ho):.0f}–{max(ho):.2f} and {min(po):.0f}–{max(po):.2f} for the others. Other callers' "
        f"hatching is too thin to see (HiRepEat {alone[('camara_hap2', 'hirepeat')][0]} and "
        f"{alone[('camara_hap1', 'hirepeat')][0]} loci, LTRharvest {alone[('camara_hap2', 'ltrharvest')][0]}, "
        f"LTR_Finder and HiTE {alone[('camara_hap2', 'ltr_finder')][0]} each on haplotype 2, Look4LTRs "
        f"{alone[('camara_hap1', 'look4ltrs')][0]} on haplotype 1). On loci two or more callers nominated, LTRquest "
        f"still leads HiRepEat by {lead[0]:.1f} and {lead[1]:.1f} points. The genotyper was developed with LTRquest "
        f"annotations, and Look4LTRs failed on haplotype 2.")
    for bad in (":", ";", "—"):                    # house style: no colons, semicolons or dash asides
        assert bad not in txt, bad
    return txt


# ============================================================================ honesty check (plain bars)
def honesty_plots(Bd, Fd, Hd):
    """plain bar charts of the numbers behind B, E and H's class bars (tmp/), to compare with the drawn forms"""
    fig, axs = plt.subplots(4, 5, figsize=(180 * MM, 200 * MM))
    xs = np.arange(len(TOOLS))
    for j, (_, c, lab) in enumerate(BCLS):
        ax = axs[0, j]
        for i, (m, _) in enumerate(CMET):
            ax.bar(xs + (i - 1) * 0.27, [Bd["V"][(c, m, t)] or 0 for t in TOOLS], 0.27, label=m,
                   color=["#999999", "#555555", "#111111"][i])
        ax.set_title(f"B {lab}", fontsize=6)
        ax.set_xticks(xs)
        ax.set_xticklabels([NAME[t] for t in TOOLS], rotation=90, fontsize=4.5)
        ax.set_ylim(0, 1)
        ax.tick_params(labelsize=5)
    axs[0, 0].legend(fontsize=5)
    for j, ds in enumerate(E_DS):
        ax = axs[1, j]
        for i, (_, key, _) in enumerate(EMET):
            ax.bar(xs + (i - 1) * 0.27, [Fd["mk"].get((ds, t), {}).get(key, 0) for t in TOOLS], 0.27,
                   color=["#999999", "#555555", "#111111"][i], label=key)
        ax.set_title(f"E {ds} (P, R, F1)", fontsize=6)
        ax.set_xticks(xs)
        ax.set_xticklabels([NAME[t] for t in TOOLS], rotation=90, fontsize=4.5)
        ax.set_ylim(0, 1)
        ax.tick_params(labelsize=5)
    axs[1, 3].axis("off")
    axs[1, 4].axis("off")
    for r, ds in ((2, "camara_hap2"), (3, "camara_hap1")):
        for j, (_, _, _, lab) in enumerate(HCOLS[1:]):
            ax = axs[r, j]
            v = [0 if Hd["o"][(ds, t)]["status"] != "ok" else
                 (Hd["o"][(ds, t)]["cls"][lab]["solid"] + Hd["o"][(ds, t)]["cls"][lab]["alone"]) /
                 Hd["o"][(ds, t)]["cls"][lab]["n"] for t in TOOLS]
            ax.bar(xs, v, color=[COL[t] for t in TOOLS])
            vs = [0 if Hd["o"][(ds, t)]["status"] != "ok" else
                  Hd["o"][(ds, t)]["cls"][lab]["solid"] / Hd["o"][(ds, t)]["cls"][lab]["n"] for t in TOOLS]
            ax.bar(xs, vs, color="none", edgecolor="k", lw=0.4)
            ax.set_title(f"H {ds[-4:]} {lab} (outline: without own-nominated)", fontsize=5)
            ax.set_xticks(xs)
            ax.set_xticklabels([NAME[t] for t in TOOLS], rotation=90, fontsize=4.5)
            ax.set_ylim(0, 1)
            ax.tick_params(labelsize=5)
    fig.tight_layout()
    fig.savefig(TMP / "honesty_plain_bars.png", dpi=200)
    plt.close(fig)


# ============================================================================ main
def main(argv=None):
    global VERBOSE
    ap = argparse.ArgumentParser(description="Build fig_benchmark (Fig. 1): PDF (figure + legend), PNG, numbers.")
    ap.add_argument("-v", "--verbose", action="store_true", help="progress and checks")
    ap.add_argument("--redo-prep", action="store_true", help="rebuild the A. thaliana locus cache (panel G)")
    ap.add_argument("--honesty", action="store_true", help="also write plain bar charts of B, E, H to tmp/")
    ap.add_argument("--draft", default=None, help="write only a PNG draft to this name inside tmp/")
    a = ap.parse_args(argv)
    VERBOSE = a.verbose
    TMP.mkdir(exist_ok=True)
    log("start", True)
    ex = load_example()
    log("A: example loaded and checked")
    Bd = load_B()
    Cd, Dd = load_curve("coding", "C"), load_curve("proteinless", "D")
    log("B-D: tables loaded and checked")
    Fd = load_EF()
    Gd = load_G(a.redo_prep)
    Hd = load_H()
    log("E-H: tables loaded and checked")
    if a.honesty:
        honesty_plots(Bd, Fd, Hd)
        log(f"honesty plots: {TMP / 'honesty_plain_bars.png'}", True)

    _, yR = h_rows()
    G["GHT"] = yR + 9 * G["RH4"] - (G["Y4"] + 7.6 + 7.0 + 1.3)     # G's heatmap ends level with H's last bar row
    Hfig = G["Y4"] + 61.0
    fig = plt.figure(figsize=(W * MM, Hfig * MM), dpi=300)      # text measured at the output resolution
    can = poafig.key_axes(fig, 0, 0, W, Hfig)
    # band 1
    names(can, G["XN"], [row_y(k) for k in range(len(TOOLS))])
    panel_A(fig, can, ex)
    panel_B(fig, can, Bd)
    # band 2
    panel_C(fig, can, Cd)
    key_CD(fig, can)
    panel_D(fig, can, Dd)
    # band 3
    panel_E_F(fig, can, Fd)
    # band 4
    panel_G(fig, can, Gd)
    panel_H(fig, can, Hd)
    for x, y, L in ((0, 0.6, "A"), (G["XB0"] - 4.3, 0.6, "B"),
                    (0, G["Y2"], "C"), (G["XD0"] - 5.5, G["Y2"], "D"),
                    (0, G["Y3"], "E"), (G["XF0"] - 3.8, G["Y3"], "F"),
                    (0, G["Y4"], "G"), (G["XH"], G["Y4"], "H")):
        poafig.panel_letter(fig, x, y, L)
    if a.draft:
        p = (TMP / Path(a.draft).name)
        fig.savefig(p, dpi=300)
        log(f"draft {p}", True)
        return
    txt = legend_text(ex, Bd, Cd, Dd, Fd, Gd, Hd)
    out = poafig.write_figure(STEM, fig, txt, None, HERE, str(B), str(HERE / f"{STEM}.py"))
    pd.DataFrame(NUM).to_csv(HERE / f"{STEM}_numbers.tsv", sep="\t", index=False)
    log(f"done {out['pdf']} fig {out['fig_mm']} legend {out['short_mm']} mm, {out['short_words']} words "
        f"{out['short_par_words']}; {len(NUM)} numbers", True)


if __name__ == "__main__":
    main()
