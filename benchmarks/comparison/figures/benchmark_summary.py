#!/usr/bin/env python3
"""benchmark_summary: the benchmark figure of the LTR_peel README ("How good is it").

Four panels of the published benchmark figure (fig_benchmark, panels A, B, C and E), re-laid out at README size
(GitHub shows it about 880 px wide; text >= 8 pt, labels a reader must read >= 9 pt, lines >= 0.6 pt).
  A  one nested host and its tenant plus a non-nested element in simulated genome S2 seed 1, every tool's calls to
     scale with the features that tool reports (fig_benchmark A: same example, same structural checks)
  B  precision / recall / F1 by nesting class and protein class, S2 mean of 3 genomes (fig_benchmark B)
  C  recall of coding elements vs median LTR divergence (K2P) in the aging series, every tool, pooled line with a 95%
     family-resampled band, one dot at each recall-halved point, documented identity floors dashed (fig_benchmark C)
  D  masking of three real genomes against curated annotation, precision / recall / F1 (fig_benchmark E)
Data loading, selection rules and checks follow fig_benchmark.py. Every number recorded for panels A-D is
cross-checked against fig_benchmark_numbers.tsv (beside this script; same quantity, same value) before writing.
Inputs under $BENCH_ROOT: experiments/figpub/{sim_class_prf_joint_summary,sim_class_prf_summary}.tsv,
experiments/nonauto/{div_elements,div_class,div_d50}.tsv, report/tables/{masking,resources}.tsv; sims/S2_1 (FASTA;
gen500000 for the example check), sims/S2_{1,2,3}/gen1000000_final.bed, sims/truth/S2_1.truth.tsv,
calls/S2_1/<tool>.tsv, runs/S2_1/<tool>/ (each tool's own output). Imported read-only: scripts/score_elements.py,
div_score.py; poafig (beside this script).
Full mode reads $BENCH_ROOT, writes the tables the figure is drawn from to source_data/benchmark_summary/, reads them
back and draws from them; --from-source-data skips $BENCH_ROOT and redraws from those tables (checks that need the
raw data run in full mode only).
Outputs: ../../../docs/images/benchmark_summary.png (300 dpi), benchmark_summary.pdf and benchmark_summary_numbers.tsv
(this directory); scratch and checks in tmp/benchmark_summary/.
Usage: benchmark_summary.py [-v] [--from-source-data] [--honesty] [--draft NAME]
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
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
os.environ["PYTHONDONTWRITEBYTECODE"] = "1"
sys.dont_write_bytecode = True                     # nothing written beside the imported benchmark scripts

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = Path(__file__).resolve().parent
STEM = "benchmark_summary"
TMP = HERE / "tmp" / STEM
SD = HERE / "source_data" / STEM                    # the tables the figure is drawn from
PNG = HERE.parents[2] / "docs" / "images" / f"{STEM}.png"     # the README figure
OLD = HERE / "fig_benchmark_numbers.tsv"            # the published figure's numbers (read for the cross-check)
B = FP = NA = TAB = RUN = se = div_score = None     # full mode only: set by bench_root()
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # poafig.py sits beside this script
import poafig  # noqa: E402
import matplotlib as mpl  # noqa: E402
import matplotlib.font_manager as fm  # noqa: E402

# fonts: Liberation registered in-process, poafig._advance_pt patched for matplotlib >= 3.11 (poafig.py itself is not
# modified)
for _f in sorted(glob.glob("/usr/share/fonts/truetype/liberation2/LiberationSans-*.ttf")
                 + glob.glob("/usr/share/fonts/truetype/liberation2/LiberationMono-*.ttf")):
    fm.fontManager.addfont(_f)
poafig.setup_mpl()
mpl.rcParams.update({"axes.linewidth": 0.7, "xtick.major.width": 0.7, "ytick.major.width": 0.7,
                     "xtick.major.size": 2.2, "ytick.major.size": 2.2, "hatch.linewidth": 0.6})   # README: >= 0.6 pt
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import to_rgb  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402


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
TOOLS = ["ltr_peel", "ltr_retriever", "ltrharvest", "ltr_finder", "look4ltrs", "hirepeat", "hite", "ltrpred",
         "dante_ltr"]
NAME = dict(ltr_peel="LTR_peel", ltr_retriever="LTR_retriever", ltrharvest="LTRharvest", ltr_finder="LTR_Finder",
            look4ltrs="Look4LTRs", hirepeat="HiRepEat", hite="HiTE", ltrpred="LTRpred", dante_ltr="DANTE_LTR")
COL = dict(ltr_peel="#D55E00", ltr_retriever="#0072B2", ltrharvest="#56B4E9", ltr_finder="#009E73",
           look4ltrs="#CC79A7", hirepeat="#E69F00", hite="#000000", ltrpred="#A6761D", dante_ltr="#6A3D9A")
DSNAME = {"arab": "A. thaliana", "rice": "rice", "finch": "zebra finch"}
DSITAL = {"arab"}                               # species names set in italic (plain text, normal word space)
MM = poafig.MM
INK, MUTED = poafig.INK, poafig.MUTED
FT = 8.0                                        # ticks, cell values, n, metric labels, small annotations
FL = 9.0                                        # tool names, headers, axis labels, key entries (must-read labels)
FLET = 12.0                                     # panel letters
LW_Q, LW_O = 1.8, 1.1                           # line weights C: LTR_peel heavier
BAND_A = 0.22                                   # band alpha C (as published)
CMAP = mpl.colormaps["viridis"]                 # score map B, D (viridis, fixed 0-1)
MISSC, HATCHC = "#F2F2F2", "#B8B8B8"              # "not evaluated" cell: light fill, gray hatch lines
NOCALLC = "#F2F2F2"                             # B: "–" cell, the caller makes no call of that class
HOSTC, HOSTC_LTR = "#5E5E5E", "#2B2B2B"         # A truth: host
TENC, TENC_LTR = "#B9B9B9", "#8A8A8A"           # A truth: tenant
FAILC = "#5A5A5A"                               # missing result, tool side: "tool failed"
DOMC = "#2E2E2E"                                # A: protein domains reported by a tool
SAMTOOLS = "/home/chris/bin/mambaforge/envs/bench_tools/bin/samtools"
# the A example (fig_benchmark A: same selection rule and checks)
EX = dict(ds="S2_1", chrom="chr1", host="chr1:11858003-11868332", tenant="chr1:11861453-11866676",
          nonnested="chr1:11885081-11890067", win=(11_857_500, 11_890_600), old_gen="gen500000",
          segs=((11_857_500, 11_869_000), (11_884_600, 11_890_600)),     # A: broken axis, same bp scale
          old_host=(11_686_375, 11_691_476))
OLDPANEL = {"A": "A", "B": "B", "C": "C", "D": "E"}   # this figure's panel -> fig_benchmark's panel (cross-check)
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


def bench_root():
    """full mode: the benchmark dataset ($BENCH_ROOT) and the benchmark's scoring modules"""
    global B, FP, NA, TAB, RUN, se, div_score
    if not os.environ.get("BENCH_ROOT"):
        raise SystemExit(f"{STEM}: set BENCH_ROOT to the benchmark dataset, or redraw with --from-source-data")
    B = Path(os.environ["BENCH_ROOT"])
    FP, NA, TAB = B / "experiments/figpub", B / "experiments/nonauto", B / "report/tables"
    RUN = B / "runs" / "S2_1"                     # the example genome's tool outputs (sequence ids c0001.. = chr1..)
    sys.path.insert(0, str(B / "scripts"))
    import score_elements as se  # noqa: F811
    import div_score  # noqa: F811  (d50 rule)


# ============================================================================ source data (the drawn tables)
def write_source_data(tables):
    """one TSV per table in SD (gzip, mtime 0, only above 1 MB uncompressed); full float precision"""
    SD.mkdir(parents=True, exist_ok=True)
    for name, df in tables.items():
        text = df.to_csv(sep="\t", index=False)
        big = len(text.encode()) > 1_000_000
        p, other = SD / f"{name}.tsv", SD / f"{name}.tsv.gz"
        if big:
            p, other = other, p
            df.to_csv(p, sep="\t", index=False, compression={"method": "gzip", "mtime": 0})
        else:
            p.write_text(text)
        if other.exists():
            other.unlink()
        log(f"source data {p.name}: {len(df)} rows")


def read_sd(name):
    """one source-data table, every column as text (callers convert)"""
    p = SD / f"{name}.tsv"
    if not p.exists():
        p = SD / f"{name}.tsv.gz"
    if not p.exists():
        raise SystemExit(f"{STEM}: source-data table {name}.tsv missing in {SD}")
    return pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)


def opt_int(v):
    return int(v) if v != "" else None


def scalars():
    """scalars.tsv as {name: value text}"""
    return {r.name: r.value for r in read_sd("scalars").itertuples()}


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
    # LTR_peel: strand + LTRs (calls tsv <- clean_ltr.tsv); TSD sequence, offsets 0,0 (clean_ltr.tsv tsd); domains
    # protein_match children of its LTR_retrotransposon in the GFF3
    g = gff_rows(RUN / "ltr_peel/S2_1_LTRs_all_depth_LTR_cleaned.gff3", seq, w0, w1)
    ids = {(r["s"], r["e"]): r["att"]["ID"] for r in g if r["type"] == "LTR_retrotransposon"}
    for f in F["ltr_peel"]:
        c = next(x for x in calls["ltr_peel"] if x["id"] == f["id"])
        if c["tsd"] not in (".", ""):
            f["tsd"] = flank(c, len(c["tsd"]))
        pid = ids[span(c)]
        f["dom"] = [(r["s"], r["e"], r["att"]["Name"]) for r in g
                    if r["type"] == "protein_match" and r["att"].get("Parent") == pid]
    SRC["ltr_peel"] = dict(strand="calls tsv strand (clean_ltr.tsv)", LTR="calls tsv lltr_end/rltr_start",
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
    fig_benchmark A extended to the non-nested element (asserted), and each tool's reported features"""
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
    for t in TOOLS:
        assert len({c["id"] for c in wc[t]}) == len(wc[t]), (t, "call ids in the window are unique")
    found_by = [t for t in TOOLS if hi in found[t]]
    tfound_by = [t for t in TOOLS if ti in found[t]]
    nfound_by = [t for t in TOOLS if ni in found[t]]
    assert "ltr_peel" in found_by and "ltr_peel" in tfound_by and "ltr_peel" in nfound_by, "no LTR_peel miss shown"
    # how typical: other tools finding each host / each non-nested element of this genome
    hosts = [i for i, x in enumerate(tr) if x["class"] == "host"]
    intact = [i for i, x in enumerate(tr) if x["class"] == "intact"]
    oth = lambda i: sum(i in found[t] for t in TOOLS[1:])     # noqa: E731
    med_h, med_n = float(np.median([oth(i) for i in hosts])), float(np.median([oth(i) for i in intact]))
    F, SRC = reported_features(wc, w0, w1)
    (a1, b1), (a2, b2) = EX["segs"]
    om = (b1, a2)                                  # omitted stretch of A's broken axis: no truth element, no call
    assert not [x for x in tr if x["chrom"] == chrom and x["end"] > om[0] and x["start"] < om[1]]
    assert not [(t, c["id"]) for t in TOOLS for c in wc[t] if c["end"] > om[0] and c["start"] < om[1]]
    for x in (h, u, nn):
        assert a1 <= x["start"] - 300 and x["end"] + 300 <= b1 or a2 <= x["start"] - 300 and x["end"] + 300 <= b2
    blank = lambda v: "" if v is None else v       # noqa: E731
    T = {}
    T["A_truth"] = pd.DataFrame([dict(role=role, id=x["id"], chrom=x["chrom"], start=x["start"], end=x["end"],
                                      strand=x["strand"], family=x["family"], superfamily=x["superfamily"],
                                      ltr_len=int(x["ltr_len"]), tsd=x["tsd"])
                                 for role, x in (("host", h), ("tenant", u), ("nonnested", nn))])
    T["A_host_pieces"] = pd.DataFrame([dict(piece=k + 1, start=a, end=b) for k, (a, b) in enumerate(pieces)])
    T["A_calls"] = pd.DataFrame([dict(tool=t, call_id=f["id"], start=f["s"], end=f["e"], matched=f["matched"],
                                      strand=blank(f["strand"]), lltr_end=blank(f["lltr"]),
                                      rltr_start=blank(f["rltr"])) for t in TOOLS for f in F[t]])
    T["A_tsd"] = pd.DataFrame([dict(tool=t, call_id=f["id"], start=a, end=b)
                               for t in TOOLS for f in F[t] for a, b in (f["tsd"] or [])],
                              columns=["tool", "call_id", "start", "end"])
    T["A_domains"] = pd.DataFrame([dict(tool=t, call_id=f["id"], start=a, end=b, name=nm)
                                   for t in TOOLS for f in F[t] for a, b, nm in f["dom"]])
    T["A_feature_sources"] = pd.DataFrame([dict(tool=t, feature=k, source=v) for t in TOOLS for k, v in SRC[t].items()])
    T["A_tool_matches"] = pd.DataFrame([dict(tool=t, host=int(t in found_by), tenant=int(t in tfound_by),
                                             nonnested=int(t in nfound_by)) for t in TOOLS])
    scal = [dict(panel="A", name="other_tools_finding_host_median_over_S2_1_hosts", value=med_h),
            dict(panel="A", name="other_tools_finding_nonnested_median_over_S2_1_nonnested", value=med_n)]
    return T, scal


def example_from_tables():
    """A as drawn, from the source-data tables; records A's numbers"""
    tr = {r.role: dict(id=r.id, chrom=r.chrom, start=int(r.start), end=int(r.end), strand=r.strand, family=r.family,
                       superfamily=r.superfamily, ltr_len=int(r.ltr_len), tsd=r.tsd)
          for r in read_sd("A_truth").itertuples()}
    h, u, nn = tr["host"], tr["tenant"], tr["nonnested"]
    pieces = [(int(r.start), int(r.end)) for r in read_sd("A_host_pieces").itertuples()]
    tsd, dom = {}, {}
    for r in read_sd("A_tsd").itertuples():
        tsd.setdefault((r.tool, r.call_id), []).append((int(r.start), int(r.end)))
    for r in read_sd("A_domains").itertuples():
        dom.setdefault((r.tool, r.call_id), []).append((int(r.start), int(r.end), r.name))
    F = {t: [] for t in TOOLS}
    for r in read_sd("A_calls").itertuples():
        F[r.tool].append(dict(s=int(r.start), e=int(r.end), matched=r.matched, id=r.call_id, strand=r.strand or None,
                              lltr=opt_int(r.lltr_end), rltr=opt_int(r.rltr_start),
                              tsd=tsd.get((r.tool, r.call_id)), dom=dom.get((r.tool, r.call_id), [])))
    SRC = {t: {} for t in TOOLS}
    for r in read_sd("A_feature_sources").itertuples():
        SRC[r.tool][r.feature] = r.source
    tm = read_sd("A_tool_matches")
    found_by, tfound_by, nfound_by = ([t for t in TOOLS if t in set(tm.tool[tm[c] == "1"])]
                                      for c in ("host", "tenant", "nonnested"))
    sc = scalars()
    med_h = float(sc["other_tools_finding_host_median_over_S2_1_hosts"])
    med_n = float(sc["other_tools_finding_nonnested_median_over_S2_1_nonnested"])
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
    num("A", "broken axis: segments drawn (0-based) and omitted stretch (no truth element, no call of any tool)",
        f"chr1:{a1}-{b1} | chr1:{a2}-{b2}; omitted {a2 - b1} bp", "this script (asserted)")
    return dict(h=h, u=u, nn=nn, pieces=pieces, F=F, SRC=SRC, found_by=found_by, tfound_by=tfound_by,
                nfound_by=nfound_by, med_h=med_h, med_n=med_n)


# ======================================================================= data: B (one scorecard, two class systems)
BSYS = [("nesting", "nesting", [("unnested", "non-nested"), ("tenant", "tenants"), ("host", "hosts")]),
        ("protein", "protein domains", [("coding", "coding"), ("proteinless", "protein-less")])]
BCLS = [(s, c, lab) for s, _, cl in BSYS for c, lab in cl]
CMET = [("precision", "prec."), ("recall", "rec."), ("f1", "F1")]


def load_B():
    """S2 (3 deeply nested genomes), mean of 3 genomes: precision, recall, F1 per class of both class systems on one
    basis (experiments/figpub/sim_class_prf_joint_summary.tsv: ext truth = every truth.tsv element,
    one-to-one RO >= 0.9; nesting classes as sim_class_prf; protein class of a true element = its library entry's
    class, of a call = the domains inside it). NA = the caller makes no call of that class (precision, F1 undefined;
    drawn "–"). Checks: nesting rows equal sim_class_prf_summary.tsv (printed strings); every finite value comes from
    3 genomes; NA only for precision / F1 with 0 genomes; n_truth per class equal for every tool."""
    want = {(s, c) for s, c, _ in BCLS}
    rows = [r for r in read_tsv(FP / "sim_class_prf_joint_summary.tsv")
            if r["set"] == "S2" and r["tool"] in TOOLS and (r["system"], r["cls"]) in want]
    assert len(rows) == len(TOOLS) * len(want), len(rows)
    old = {(r["tool"], r["cls"]): r for r in read_tsv(FP / "sim_class_prf_summary.tsv") if r["set"] == "S2"}
    rec, n = [], {}
    for r in rows:
        t, c = r["tool"], r["cls"]
        for m, _ in CMET:
            v = r[f"{m}_mean"]
            if v == "NA":
                assert m != "recall" and r[f"{m}_n_genomes"] == "0" and int(r["n_calls"]) == 0, (t, c, m)
                rec.append(dict(tool=t, system=r["system"], cls=c, metric=m, mean="", sd="", n_true=int(r["n_truth"])))
                continue
            assert r[f"{m}_n_genomes"] == "3", (t, c, m)
            rec.append(dict(tool=t, system=r["system"], cls=c, metric=m, mean=float(v), sd=float(r[f"{m}_sd"]),
                            n_true=int(r["n_truth"])))
        if r["system"] == "nesting":                 # the nesting half = fig_benchmark's nesting scorecard, unchanged
            o = old[(t, c)]
            for m, _ in CMET:
                assert (r[f"{m}_mean"], r[f"{m}_sd"]) == (o[f"{m}_mean"], o[f"{m}_sd"]), (t, c, m)
        assert n.setdefault(c, int(r["n_truth"])) == int(r["n_truth"]), (t, c)
    return {"B_scorecard": pd.DataFrame(rec)}


def B_from_tables():
    """B as drawn, from the source-data table (rows in the order of sim_class_prf_joint_summary.tsv); records B's
    numbers"""
    src = "experiments/figpub/sim_class_prf_joint_summary.tsv set=S2"
    V, SD_, n = {}, {}, {}
    for r in read_sd("B_scorecard").itertuples():
        t, c, m = r.tool, r.cls, r.metric
        n[c] = int(r.n_true)
        if r.mean == "":
            V[(c, m, t)] = None
            num("B", f"S2 {NAME[t]} {c} {m}", "NA (no calls of this class; drawn '–')",
                f"{src} tool={t} cls={c} {m}_mean (n_calls 0)")
            continue
        V[(c, m, t)], SD_[(c, m, t)] = float(r.mean), float(r.sd)
        num("B", f"S2 {NAME[t]} {c} {m} (mean of 3 genomes; SD {float(r.sd):.4f})", f"{float(r.mean):.4f}",
            f"{src} tool={t} system={r.system} cls={c} {m}_mean, {m}_sd")
    for s, c, lab in BCLS:
        num("B", f"true elements in class {lab} (S2, 3 genomes summed; header n)", n[c],
            f"{src} n_truth (system {s}, cls {c})")
    # spread across the three genomes (one legend clause): largest SD per tool over the drawn classes and metrics
    sd_max = {t: max(v for (c, m, tt), v in SD_.items() if tt == t) for t in TOOLS}
    for t in TOOLS:
        num("B", f"S2 {NAME[t]} largest across-genome SD over the drawn classes and metrics", f"{sd_max[t]:.4f}",
            f"{src} *_sd")
    return dict(V=V, SD=SD_, n=n, sd_max=sd_max)


# ============================================================================ data: C (aging series)
def k2p_of_identity(p100):
    """LTR-LTR identity (%) -> K2P under the simulation's model (substitutions only, transitions = transversions)"""
    p = 1 - p100 / 100
    P = Q = p / 2
    return -0.5 * math.log(1 - 2 * P - Q) - 0.25 * math.log(1 - 2 * Q)


FLOOR = {"ltr_peel": None, "ltr_retriever": 85, "ltrharvest": 85, "ltr_finder": 85, "look4ltrs": None,
         "hirepeat": None, "hite": 85, "ltrpred": 70, "dante_ltr": 70}


def pooled_d50(xs, rs):
    """K2P at which the pooled recall falls to half its timepoint-0 value: scripts/div_score.d50 (linear between
    timepoints), on the pooled curve"""
    return div_score.d50([(t, float(x), float(r)) for t, (x, r) in enumerate(zip(xs, rs))])


def load_curve(cls, tools=TOOLS):
    """one protein class of the aging series (coding for C; protein-less only for one legend clause), 3 genomes
    pooled: experiments/figpub/div_all_ci.py's procedure on one class (as the published panels). Per tool x timepoint
    recall over the pooled elements of the class; 95% CI from 1,000 resamples of their library families (numpy
    default_rng(1), every copy of a drawn family carried along, percentiles 2.5 / 97.5 as nonauto_stats.ci); x = median
    K2P of the pooled elements of the class at the timepoint (each element once); one dot at the pooled recall-halved
    point (scripts/div_score.d50, linear between timepoints). Checks: pooled n / tp = div_class.tsv sums over seeds;
    the pooled d50 agrees with div_d50.tsv within 0.002 and lies inside its CI. Returns the curve rows (tools), the
    d50 rows (tools) and the class's family / element counts."""
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
    curve, d50 = [], []
    for t in TOOLS:
        st = s_[s_.tool == t]
        N, Hh = np.zeros((nF, 11)), np.zeros((nF, 11))
        f_idx = st.family.map(fi).to_numpy()
        np.add.at(N, (f_idx, st.timepoint.to_numpy()), 1)
        np.add.at(Hh, (f_idx, st.timepoint.to_numpy()), st.hit.to_numpy())
        n, h = N.sum(0), Hh.sum(0)
        g = dc[dc.tool == t].groupby("timepoint")[["n", "tp"]].sum()
        if not (np.array_equal(g.n.to_numpy(), n) and np.array_equal(g.tp.to_numpy(), h)):
            raise SystemExit(f"{STEM}: pooled {cls} counts {t} differ from div_class.tsv")
        r = h / n
        with np.errstate(invalid="ignore", divide="ignore"):
            bs = (M @ Hh) / (M @ N)
        lo, hi = np.nanpercentile(bs, 2.5, axis=0), np.nanpercentile(bs, 97.5, axis=0)
        v = pooled_d50(x, r)
        rr = ref[(ref.tool == t) & (ref.cls == cls)].iloc[0]
        if v[0].isdigit() and rr.d50[:1].isdigit():
            assert abs(float(v) - float(rr.d50)) <= 0.002, (t, cls, v, rr.d50)
            assert float(rr.d50_lo) <= float(v) <= float(rr.d50_hi), (t, cls, v, rr.d50_lo, rr.d50_hi)
        if t not in tools:
            continue
        d50.append(dict(cls=cls, tool=t, d50_pooled=v, d50_table=rr.d50, d50_table_lo=rr.d50_lo,
                        d50_table_hi=rr.d50_hi))
        curve += [dict(cls=cls, tool=t, timepoint=tt, median_k2p=float(x[tt]), recall=float(r[tt]),
                       ci_lo=float(lo[tt]), ci_hi=float(hi[tt]), n=int(n[tt])) for tt in range(11)]
    n_el = int(s_[(s_.tool == TOOLS[0]) & (s_.timepoint == 0)].shape[0])
    return curve, d50, nF, n_el


def C_from_tables():
    """C as drawn (coding curves, dots) and the protein-less d50 of the legend, from the source-data tables;
    records C's numbers"""
    cv = read_sd("C_curves")
    dd = {(r.cls, r.tool): r for r in read_sd("C_d50").itertuples()}
    sc = scalars()
    S = {}
    for t in TOOLS:
        c = cv[(cv.cls == "coding") & (cv.tool == t)]
        assert list(c.timepoint) == [str(k) for k in range(11)], t
        x, r, lo, hi = (c[k].astype(float).to_numpy() for k in ("median_k2p", "recall", "ci_lo", "ci_hi"))
        n = c.n.astype(int).to_numpy()
        rr = dd[("coding", t)]
        S[t] = dict(x=x, r=r, lo=lo, hi=hi, d50=rr.d50_pooled, r0=r[0])
        num("C", f"{NAME[t]} coding pooled d50 (K2P at which pooled recall halves; dot)", rr.d50_pooled,
            "computed: scripts/div_score.d50 on pooled coding counts (experiments/nonauto/div_elements.tsv)")
        num("C", f"{NAME[t]} coding d50 in div_d50.tsv [95% CI] (x = class median over seeds), for comparison",
            f"{rr.d50_table} [{rr.d50_table_lo}, {rr.d50_table_hi}]", "experiments/nonauto/div_d50.tsv")
        for tt in range(11):
            num("C", f"{NAME[t]} timepoint {tt} coding pooled recall [95% CI] at median K2P {x[tt]:.4f}",
                f"{r[tt]:.4f} [{lo[tt]:.4f}, {hi[tt]:.4f}] (n {int(n[tt])})",
                "computed from experiments/nonauto/div_elements.tsv cls=coding (= div_class.tsv sums)")
    num("C", "coding element families (unit of the resampling)", int(sc["coding_element_families"]),
        "experiments/nonauto/div_elements.tsv")
    num("C", "coding elements per timepoint (3 genomes pooled)", int(sc["coding_elements_per_timepoint"]),
        "experiments/nonauto/div_elements.tsv")
    for fl in (85, 70):
        num("C", f"documented identity floor {fl}% as K2P (simulation model)", f"{k2p_of_identity(fl):.3f}",
            "computed (substitutions only, transitions = transversions); fig2_divergence memo")
    Cp = {"S": {t: dict(d50=r.d50_pooled) for (c, t), r in dd.items() if c == "proteinless"}}
    return dict(S=S, n_fam=int(sc["coding_element_families"])), Cp


# ============================================================================ data: D (real genomes; fig_benchmark E)
E_DS = ["arab", "rice", "finch"]                  # masking P / R / F1 (maize and gorilla are in the paper's Fig. S2C)
EMET = [("prec", "precision", "prec."), ("sens", "recall", "rec."), ("F1", "f1", "F1")]
# missing cells as drawn, asserted (reasons: resources.tsv status / masking.tsv 'incomplete')
E_MISSING = {"failed": {("finch", "ltr_retriever"), ("finch", "dante_ltr")},
             "noteval": {("finch", "hirepeat")}}


def missing_reasons():
    """(ds, tool) -> ("failed" | "noteval", source). failed = the tool did not finish within 36 h (timeout), exceeded
    200 GB (mem_guard) or crashed with no calls (failed): resources.tsv status. noteval = the tool finished but our
    masking or library build with its library would have run past 36 h or was stopped at 36 h: masking.tsv ci
    'incomplete: ...'. Every (ds, tool) gets at most one class (asserted)."""
    res = pd.read_csv(TAB / "resources.tsv", sep="\t")
    m = pd.read_csv(TAB / "masking.tsv", sep="\t")
    out = {}
    for ds in E_DS:
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


def load_D():
    """masking precision, recall (lib-test sensitivity), F1 of each tool's intact library against curated annotation
    (report/tables/masking.tsv role primary, <tool>.intact: prec, sens, F1) for A. thaliana, rice, zebra finch
    (fig_benchmark E). Missing cells classed and asserted."""
    m = pd.read_csv(TAB / "masking.tsv", sep="\t")
    m = m[(m.role == "primary") & (m.F1 != ".")].copy()
    m["tool"] = m.lib.str.replace(".intact", "", regex=False)
    m = m[m.tool.isin(TOOLS)]
    mk = {(r.ds, r.tool): dict(precision=float(r.prec), recall=float(r.sens), f1=float(r.F1), f1_lo=r.F1_lo,
                               f1_hi=r.F1_hi) for r in m.itertuples() if r.ds in E_DS}
    reasons = missing_reasons()
    check_missing("D", E_DS, mk, reasons, E_MISSING)
    return {"D_masking": pd.DataFrame([dict(ds=ds, tool=t, **v) for (ds, t), v in mk.items()]),
            "D_missing": pd.DataFrame([dict(ds=ds, tool=t, kind=k, reason=why)
                                       for (ds, t), (k, why) in reasons.items()])}


def D_from_tables():
    """D as drawn, from the source-data tables; records D's numbers"""
    mk = {(r.ds, r.tool): dict(precision=float(r.precision), recall=float(r.recall), f1=float(r.f1), f1_lo=r.f1_lo,
                               f1_hi=r.f1_hi) for r in read_sd("D_masking").itertuples()}
    reasons = {(r.ds, r.tool): (r.kind, r.reason) for r in read_sd("D_missing").itertuples()}
    for kind in ("failed", "noteval"):
        for ds, t in sorted(k for k, v in reasons.items() if v[0] == kind and k not in mk):
            num("D", f"{ds} {NAME[t]} missing result class (cell text)",
                "tool failed" if kind == "failed" else "not evaluated", reasons[(ds, t)][1])
    for (ds, t), v in sorted(mk.items()):
        for col, key, _ in EMET:
            extra = f" [95% CI {v['f1_lo']}, {v['f1_hi']}]" if key == "f1" else ""
            num("D", f"{ds} {NAME[t]} masking {key}{extra}", f"{v[key]:.4f}",
                f"report/tables/masking.tsv {col} (role primary, {t}.intact)")
    return dict(mk=mk, reasons=reasons)


# ============================================================================ drawing helpers
def text_color(v):
    """white on dark cells, ink on light ones: relative luminance below 0.20 (equal contrast of white and ink)"""
    r, g, b = (c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4 for c in to_rgb(CMAP(v)))
    return "white" if 0.2126 * r + 0.7152 * g + 0.0722 * b < 0.20 else INK


def score_cell(ax, x, y, w, h, v, bold, pad=0.15, fs=FT):
    ax.add_patch(Rectangle((x + pad, y + pad), w - 2 * pad, h - 2 * pad, fc=CMAP(v), ec="none", lw=0))
    ax.text(x + w / 2, y + h / 2, f"{v:.2f}", ha="center", va="center", fontsize=fs, color=text_color(v),
            fontweight="bold" if bold else "normal")


def missing_cell(ax, x, y, w, h, kind, pad=0.15, fs=FT):
    """two styles: tool side "tool failed" (dark gray, white text); evaluation side "not evaluated" (light hatch)"""
    if kind == "failed":
        ax.add_patch(Rectangle((x + pad, y + pad), w - 2 * pad, h - 2 * pad, fc=FAILC, ec="none", lw=0))
        ax.text(x + w / 2, y + h / 2, "tool failed", ha="center", va="center", fontsize=fs, color="white")
    else:
        # hatch drawn as explicit 45-degree lines clipped to the cell: a matplotlib hatch pattern renders much fainter
        # in the PDF than in the PNG (checked with pdftoppm), so PNG and PDF match this way (ax: 1 unit = 1 mm)
        x0, y0, w0, h0 = x + pad, y + pad, w - 2 * pad, h - 2 * pad
        cell = Rectangle((x0, y0), w0, h0, fc=MISSC, ec="none", lw=0)
        ax.add_patch(cell)
        for xa in np.arange(x0 - h0, x0 + w0, 0.9):
            ln, = ax.plot([xa, xa + h0], [y0 + h0, y0], color=HATCHC, lw=0.6, solid_capstyle="butt")
            ln.set_clip_path(cell)
        ax.text(x + w / 2, y + h / 2, "not evaluated", ha="center", va="center", fontsize=fs, color=MUTED,
                bbox=dict(boxstyle="square,pad=0.08", fc=MISSC, ec="none"))


def nocall_cell(ax, x, y, w, h, pad=0.15, fs=FT):
    """B: the caller makes no call of the class, so precision and F1 are undefined: light neutral cell with '–'"""
    ax.add_patch(Rectangle((x + pad, y + pad), w - 2 * pad, h - 2 * pad, fc=NOCALLC, ec="none", lw=0))
    ax.text(x + w / 2, y + h / 2, "–", ha="center", va="center", fontsize=fs, color=MUTED)


def best_set(vals):
    """tools whose printed (2-decimal) value equals the column's best printed value (ties all bold); no bold when
    fewer than two tools have a value in the column"""
    pv = {t: float(f"{v:.2f}") for t, v in vals.items()}
    if len(pv) < 2:
        return set()
    b = max(pv.values())
    return {t for t, v in pv.items() if v == b}


def names(ax, x_right, ys, rows=TOOLS, fs=FL):
    """tool names (LTR_peel bold) and a square of the tool color, right-aligned to x_right (page mm)"""
    for y, t in zip(ys, rows):
        ax.text(x_right - 2.7, y, NAME[t], ha="right", va="center", fontsize=fs,
                fontweight="bold" if t == "ltr_peel" else "normal")
        ax.add_patch(Rectangle((x_right - 2.0, y - 0.85), 1.7, 1.7, fc=COL[t], ec="none"))


def text_mm(fig, s, fs, bold=False):
    """rendered width of a text in mm (figure made at its output dpi)"""
    t = fig.text(0, 0, s, fontsize=fs, fontweight="bold" if bold else "normal")
    w = t.get_window_extent(renderer=fig.canvas.get_renderer()).width / fig.dpi * 25.4
    t.remove()
    return w


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


def chevron(ax, x, y, h, strand, col, bp_per_mm, size=0.8, lw=0.8):
    d = size * bp_per_mm * (1 if strand == "+" else -1)
    ax.plot([x - d / 2, x + d / 2, x - d / 2], [y - h * 0.30, y, y + h * 0.30], color=col, lw=lw,
            solid_capstyle="butt", solid_joinstyle="miter", zorder=5)


def clean(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="both", length=2.2, pad=1.5, labelsize=FT)


def letter(fig, x, y, L):
    W, H = fig.get_size_inches() / MM
    fig.text(x / W, 1 - y / H, L, fontsize=FLET, fontweight="bold", va="top", ha="left")


# ============================================================================ layout (mm from the top-left)
W = 180.0
G = dict(
    XN=23.5,                       # right edge of the tool-name column (A, B) incl. the color square
    X0=25.0, X1=177.5,             # data left edge (A tracks, B cells, C plot) and right edge (A tracks, D cells)
    # band 1: A
    KEYA=2.6,                      # glyph key row (center)
    YH=11.0, YT=15.4, EH=2.4,      # truth: host / non-nested level, tenant level, element height
    ROWA0=18.8, RHA=4.6,           # A tool rows (two sub-rows each: host-level calls above, tenant calls below)
    # band 2: B
    Y2=73.2, RH=4.5,               # B and D share row height and column width (D as plain as B)
    CW=7.4, GGAP=1.6, BSYS=4.2,    # column width, gap between class / genome groups, gap between B's two systems
    # band 3: C | D
    Y3=0.0,                       # set in layout(): B rows end + 5.4
    KEYP=3.6,                      # C key row pitch
)


def layout(fig):
    """derived positions: B rows, D columns and names, C plot edges, band 3 header from the measured C key"""
    G["ROWB0"] = G["Y2"] + 15.4
    G["Y3"] = G["ROWB0"] + 9 * G["RH"] + 5.4
    G["XD0"] = G["X1"] - (9 * G["CW"] + 2 * G["GGAP"])
    G["XDN"] = G["XD0"] - 0.8
    longest = max(text_mm(fig, NAME[t], FL, bold=t == "ltr_peel") for t in TOOLS)
    G["XDL"] = G["XDN"] - 2.7 - longest - 4.2          # letter D: left of D's longest name
    G["XC1"] = G["XDL"] + 0.6                           # C plot right edge (D's names start 3.6 mm further right)
    rows = key_C(fig, None, dry=True)
    G["ROWD0"] = G["Y3"] + max(12.8, 1.0 + rows * G["KEYP"] + 1.2)
    G["HFIG"] = G["ROWD0"] + 9 * G["RH"] + 10.2


def row_y(k):
    return G["ROWA0"] + (k + 0.5) * G["RHA"]


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
    """fig_benchmark A at README scale: nesting diagram (host in two pieces, its tenant) and a non-nested element on a
    broken axis (one bp scale), then each tool's calls with what that tool reports; host-level calls above the row
    line, tenant calls below"""
    (a1, b1), (a2, b2) = EX["segs"]
    XL, XR = G["X0"], G["X1"]
    GAP = 3.5                                      # mm of the axis break
    bppm = ((b1 - a1) + (b2 - a2)) / (XR - XL - GAP)          # one bp scale for both segments
    gap_bp = GAP * bppm

    def X(v):
        """genome bp -> axis x (bp units, the omitted stretch replaced by the break)"""
        if a1 <= v <= b1:
            return v - a1
        assert a2 <= v <= b2, v
        return (b1 - a1) + gap_bp + (v - a2)
    ROW0, RH = G["ROWA0"], G["RHA"]
    TOP = 4.8
    xmax = X(b2)
    ax = bp_axes(fig, XL, TOP, XR - XL, ROW0 + 9 * RH - TOP + 0.6, (0, xmax))
    h, u, nn = ex["h"], ex["u"], ex["nn"]
    (p1s, p1e), (p2s, p2e) = ex["pieces"]
    yh, yt, eh = G["YH"], G["YT"], G["EH"]

    def el(s0, e0, y, hh, body, ltr=None, lltr=None, rltr=None, **kw):
        element(ax, X(s0), X(e0), y, hh, body, ltr, lltr=None if lltr is None else X(lltr),
                rltr=None if rltr is None else X(rltr), **kw)
    # truth: nesting diagram (host in two pieces, its tenant) and a non-nested element, real coordinates
    L = int(h["ltr_len"])
    el(p1s, p1e, yh, eh, HOSTC, HOSTC_LTR, lltr=p1s + L)
    el(p2s, p2e, yh, eh, HOSTC, HOSTC_LTR, rltr=p2e - L)
    yb = yh - eh / 2 - 1.4
    ax.plot([X(p1e) - 0.3 * bppm, X(p1e) - 0.3 * bppm, X(p2s) + 0.3 * bppm, X(p2s) + 0.3 * bppm],
            [yh - eh / 2, yb, yb, yh - eh / 2], color=HOSTC, lw=0.8, ls=(0, (2.6, 1.6)), zorder=3)
    ax.text(X(p1s), yb - 0.5, "host", ha="left", va="baseline", fontsize=FL, color=INK)
    chevron(ax, X((p1s + L + p1e) / 2), yh, eh, h["strand"], "white", bppm, size=0.9, lw=0.9)
    Lt = int(u["ltr_len"])
    el(u["start"], u["end"], yt, eh, TENC, TENC_LTR, lltr=u["start"] + Lt, rltr=u["end"] - Lt)
    chevron(ax, X((u["start"] + u["end"]) / 2), yt, eh, u["strand"], INK, bppm, size=0.9, lw=0.9)
    ax.text(X(u["end"]) + 1.0 * bppm, yt, "tenant", ha="left", va="center", fontsize=FL, color=INK)
    for x in (u["start"], u["end"]):
        ax.plot([X(x), X(x)], [yh + eh / 2, yt - eh / 2], color="#9A9A9A", lw=0.6, zorder=0)
    Ln = int(nn["ltr_len"])
    el(nn["start"], nn["end"], yh, eh, HOSTC, HOSTC_LTR, lltr=nn["start"] + Ln, rltr=nn["end"] - Ln)
    chevron(ax, X((nn["start"] + nn["end"]) / 2), yh, eh, nn["strand"], "white", bppm, size=0.9, lw=0.9)
    ax.text(X((nn["start"] + nn["end"]) / 2), yh - eh / 2 - 0.9, "non-nested", ha="center", va="baseline",
            fontsize=FL, color=INK)
    # each tool: what it reports (LTR end blocks, TSD boxes, strand chevron, protein domains)
    bh = 1.6
    off = 1.05
    tsd_w = 0.8 * bppm                             # TSDs (4-6 bp) drawn 0.8 mm wide so they can be seen
    for k, t in enumerate(TOOLS):
        y = row_y(k)
        for sa, sb in EX["segs"]:
            ax.plot([X(sa), X(sb)], [y, y], color="#E2E2E2", lw=0.6, zorder=0)
        for f in ex["F"][t]:
            yy = y + off if f["matched"] == u["id"] else y - off
            body = tint(COL[t], 0.38) if f["lltr"] is not None else tint(COL[t], 0.75)
            el(f["s"], f["e"], yy, bh, body, COL[t], lltr=f["lltr"], rltr=f["rltr"])
            for a, b, _ in f["dom"]:
                ax.add_patch(Rectangle((X(a), yy - bh * 0.28), b - a, bh * 0.56, fc=DOMC, ec="none", zorder=3))
            for a, b in (f["tsd"] or []):
                xa = X(f["s"]) - tsd_w if b <= f["s"] + 1 else X(f["e"])
                ax.add_patch(Rectangle((xa, yy - bh * 0.38), tsd_w, bh * 0.76, fc="white", ec=INK, lw=0.6,
                                       zorder=4))
            if f["strand"]:
                lo = f["lltr"] if f["lltr"] is not None else f["s"]
                hi = f["rltr"] if f["rltr"] is not None else f["e"]
                xc = free_center(lo, hi, [(a, b) for a, b, _ in f["dom"]])
                chevron(ax, X(xc), yy, bh, f["strand"], INK, bppm)
    # broken axis: two segments, "//" at the break, a faint break across the tracks, the omitted length labeled
    yax = ROW0 + 9 * RH + 0.6
    xg0, xg1 = X(b1), X(a2)
    for sa, sb in EX["segs"]:
        ax.plot([X(sa), X(sb)], [yax, yax], color=INK, lw=0.7, clip_on=False)
        for x in np.arange(np.ceil(sa / 5000) * 5000, sb + 1, 5000):
            ax.plot([X(x), X(x)], [yax, yax + 1.0], color=INK, lw=0.7, clip_on=False)
            if abs(X(x) - (xg0 + xg1) / 2) > 11.0 * bppm:         # leave room for the omitted-length label
                ax.text(X(x), yax + 1.5, f"{x / 1e6:.3f}", ha="center", va="top", fontsize=FT)
    for xb in (xg0 + 0.22 * gap_bp, xg1 - 0.22 * gap_bp):
        ax.plot([xb - 0.4 * bppm, xb + 0.4 * bppm], [yax + 0.9, yax - 0.9], color=INK, lw=0.8, clip_on=False)
        ax.plot([xb, xb], [TOP + 1.2, yax - 1.0], color="#C8C8C8", lw=0.6, ls=(0, (1.4, 1.4)), zorder=0)
    ax.text((xg0 + xg1) / 2, yax + 1.5, f"{(a2 - b1) / 1000:.1f} kb", ha="center", va="top", fontsize=FT,
            color=MUTED)
    ax.text(xmax / 2, yax + 5.0, "simulated chromosome 1 (Mb)", ha="center", va="top", fontsize=FL)
    can.text(G["XN"] - 2.7, (yh + yt) / 2, "true elements", ha="right", va="center", fontsize=FL)
    # compact glyph key above the panel
    ky, x = G["KEYA"], XL
    can.add_patch(Rectangle((x, ky - 0.8), 4.2, 1.6, fc=tint("#7A7A7A", 0.38), ec="none"))
    can.add_patch(Rectangle((x, ky - 0.8), 1.0, 1.6, fc="#7A7A7A", ec="none"))
    can.add_patch(Rectangle((x + 3.2, ky - 0.8), 1.0, 1.6, fc="#7A7A7A", ec="none"))
    can.text(x + 5.0, ky, "LTRs", ha="left", va="center", fontsize=FL)
    x += 5.0 + text_mm(fig, "LTRs", FL) + 3.6
    can.add_patch(Rectangle((x, ky - 0.6), 0.8, 1.2, fc="white", ec=INK, lw=0.6))
    can.text(x + 1.6, ky, "TSD", ha="left", va="center", fontsize=FL)
    x += 1.6 + text_mm(fig, "TSD", FL) + 3.6
    can.plot([x, x + 0.8, x], [ky - 0.55, ky, ky + 0.55], color=INK, lw=0.8, solid_joinstyle="miter")
    can.text(x + 1.6, ky, "strand", ha="left", va="center", fontsize=FL)
    x += 1.6 + text_mm(fig, "strand", FL) + 3.6
    can.add_patch(Rectangle((x, ky - 0.45), 2.2, 0.9, fc=DOMC, ec="none"))
    can.text(x + 3.0, ky, "protein domain", ha="left", va="center", fontsize=FL)


def panel_B(fig, can, Bd):
    """one scorecard: precision, recall, F1 for the nesting classes (non-nested, tenants, hosts) and the protein
    classes (coding, protein-less); system super-header with a rule over its groups, class name, n (true elements,
    3 genomes summed), group rule, metric labels; viridis 0-1, bold = column best"""
    Y, R0, RH = G["Y2"], G["ROWB0"], G["RH"]
    w, gap, sysgap = G["CW"], G["GGAP"], G["BSYS"]
    names(can, G["XN"], [R0 + (k + 0.5) * RH for k in range(len(TOOLS))])
    x = G["X0"]
    Bd["drawn_nocall"] = []
    for si, (sysk, syslab, classes) in enumerate(BSYS):
        if si:
            x += sysgap - gap
        sx0 = x
        for ci, (cls, clab) in enumerate(classes):
            gx0 = x
            for m, mlab in CMET:
                can.text(x + w / 2, R0 - 2.0, mlab, ha="center", va="center", fontsize=FT)
                vals = {t: Bd["V"][(cls, m, t)] for t in TOOLS if Bd["V"][(cls, m, t)] is not None}
                best = best_set(vals)
                for k, t in enumerate(TOOLS):
                    if t in vals:
                        score_cell(can, x, R0 + k * RH, w, RH, vals[t], t in best)
                    else:
                        nocall_cell(can, x, R0 + k * RH, w, RH)
                        Bd["drawn_nocall"].append((cls, m, t))
                x += w
            gw = x - gx0
            can.text(gx0 + gw / 2, Y + 5.9, clab, ha="center", va="center", fontsize=FL)
            can.text(gx0 + gw / 2, Y + 9.3, f"n = {Bd['n'][cls]:,}", ha="center", va="center", fontsize=FT,
                     color=MUTED)
            can.plot([gx0 + 0.5, x - 0.5], [Y + 11.3, Y + 11.3], color=INK, lw=0.6)
            x += gap
        sx1 = x - gap
        can.text((sx0 + sx1) / 2, Y + 1.4, syslab, ha="center", va="center", fontsize=FL, color=MUTED)
        can.plot([sx0 + 0.5, sx1 - 0.5], [Y + 3.4, Y + 3.4], color=INK, lw=0.6)
    G["XB1"] = x - gap
    assert Bd["drawn_nocall"] == [("proteinless", "precision", "dante_ltr"), ("proteinless", "f1", "dante_ltr")], \
        Bd["drawn_nocall"]


def colorbar(can):
    """the one score scale of B and D (viridis, fixed 0-1), right of B, above D"""
    w, h = 2.8, 26.0
    x = G["XB1"] + 13.0
    y0 = G["ROWB0"] + 9 * G["RH"] - h - 2.0
    n = 120
    for i in range(n):
        v = 1 - (i + 0.5) / n                      # top = 1
        can.add_patch(Rectangle((x, y0 + i * h / n), w, h / n + 0.02, fc=CMAP(v), ec="none"))
    can.add_patch(Rectangle((x, y0), w, h, fc="none", ec=INK, lw=0.6))
    for v in (0, 0.5, 1):
        yy = y0 + (1 - v) * h
        can.plot([x + w, x + w + 1.0], [yy, yy], color=INK, lw=0.7)
        can.text(x + w + 1.6, yy, f"{v:g}", ha="left", va="center", fontsize=FT)
    can.text(x + w / 2, y0 - 3.0, "score (B, D)", ha="center", va="baseline", fontsize=FL)


def key_C(fig, can, dry=False):
    """C's key above its plot: the nine tools in the shared order (line samples, LTR_peel heavier), then the 95% band
    and the recall-halved dot, in a grid of aligned columns (as many as fit between the letter C and the letter D);
    returns the number of rows"""
    x0, x1 = 4.6, G["XDL"] - 1.8
    items = [(t, NAME[t]) for t in TOOLS] + [("ci", "95% CI"), ("dot", "recall halved")]
    gl, g1, g2 = 3.6, 0.9, 2.4                    # glyph length, glyph-text gap, gap between columns
    ew = [gl + g1 + text_mm(fig, lab, FL, bold=kind == "ltr_peel") for kind, lab in items]
    for nc in range(len(items), 0, -1):           # most columns whose widest entries fit
        cw = [max(ew[i::nc]) for i in range(nc)]
        if sum(cw) + g2 * (nc - 1) <= x1 - x0:
            break
    nr = -(-len(items) // nc)
    if not dry:
        for i, (kind, lab) in enumerate(items):
            r, c = divmod(i, nc)
            x = x0 + sum(cw[:c]) + g2 * c
            y = G["Y3"] + 1.0 + (r + 0.5) * G["KEYP"]
            if kind in TOOLS:
                can.plot([x, x + gl], [y, y], color=COL[kind], lw=LW_Q if kind == "ltr_peel" else LW_O,
                         solid_capstyle="butt")
            elif kind == "ci":
                can.add_patch(Rectangle((x, y - 0.8), gl, 1.6, fc="#8C8C8C", alpha=BAND_A * 1.6, ec="none"))
            else:
                can.plot([x + gl / 2], [y], "o", ms=3.6, mfc=INK, mec="white", mew=0.6)
            can.text(x + gl + g1, y, lab, ha="left", va="center", fontsize=FL,
                     fontweight="bold" if kind == "ltr_peel" else "normal")
    return nr


def panel_C(fig, can, Cd):
    """fig_benchmark C at README scale: pooled recall of coding elements vs median K2P, one line + 95% band per tool,
    one dot at the pooled recall-halved point, documented identity floors dashed; the plot spans D's rows"""
    x0, x1 = G["X0"], G["XC1"]
    y0, h = G["ROWD0"], 9 * G["RH"]
    ax = poafig.mm_axes(fig, x0, y0, x1 - x0, h)
    for fl, lab, ha, dx in ((85, "85% identity", "left", 0.006), (70, "70%", "right", -0.006)):
        xf = k2p_of_identity(fl)
        ax.axvline(xf, color=INK, lw=0.7, ls=(0, (2.6, 1.8)), zorder=1)
        ax.text(xf + dx, 1.155, lab, ha=ha, va="top", fontsize=FL)
    for t in TOOLS[::-1]:
        v = Cd["S"][t]
        ax.fill_between(v["x"], v["lo"], v["hi"], color=COL[t], alpha=BAND_A, lw=0, zorder=2)
        ax.plot(v["x"], v["r"], color=COL[t], lw=LW_Q if t == "ltr_peel" else LW_O, zorder=3 + (t == "ltr_peel"),
                solid_capstyle="round")
        if v["d50"] not in ("NA", "incomplete") and not v["d50"].startswith(">"):
            ax.plot(float(v["d50"]), v["r0"] / 2, "o", ms=3.6, mfc=COL[t], mec="white", mew=0.6, zorder=7)
    ax.set_xlim(0, 0.40)
    ax.set_ylim(-0.03, 1.16)                 # headroom above 1 for the floor labels; lines at 0 stay visible
    ax.spines["left"].set_bounds(-0.03, 1.0)
    ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4])
    ax.set_xticklabels(["0", "0.1", "0.2", "0.3", "0.4"])
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_yticklabels(["0", "0.25", "0.5", "0.75", "1"])
    clean(ax)
    ax.set_xlabel("median LTR divergence (K2P)", fontsize=FL, labelpad=2.0)
    ax.set_ylabel("coding elements\nfound (recall)", fontsize=FL, labelpad=2.5, linespacing=1.15)
    return ax


def panel_D(fig, can, Dd):
    """fig_benchmark E at README scale, on the same cell size and fonts as B: masking precision / recall / F1 for three
    genomes (blocks of three columns, a missing result spans its genome's block)"""
    R0, RH, w, gap = G["ROWD0"], G["RH"], G["CW"], G["GGAP"]
    names(can, G["XDN"], [R0 + (k + 0.5) * RH for k in range(len(TOOLS))])
    drawn = []
    for j, ds in enumerate(E_DS):
        gx = G["XD0"] + j * (3 * w + gap)
        gw = 3 * w
        can.text(gx + gw / 2, R0 - 7.0, DSNAME[ds], ha="center", va="center", fontsize=FL,
                 fontstyle="italic" if ds in DSITAL else "normal")
        can.plot([gx + 0.5, gx + gw - 0.5], [R0 - 4.6, R0 - 4.6], color=INK, lw=0.6)
        for i, (_, key, mlab) in enumerate(EMET):
            x = gx + i * w
            can.text(x + w / 2, R0 - 2.0, mlab, ha="center", va="center", fontsize=FT)
            col = {t: Dd["mk"][(ds, t)][key] for t in TOOLS if (ds, t) in Dd["mk"]}
            best = best_set(col)
            for k, t in enumerate(TOOLS):
                if t in col:
                    score_cell(can, x, R0 + k * RH, w, RH, col[t], t in best)
        for k, t in enumerate(TOOLS):
            if (ds, t) not in Dd["mk"]:
                kind = Dd["reasons"][(ds, t)][0]
                missing_cell(can, gx, R0 + k * RH, gw, RH, kind)
                drawn.append((ds, t, kind))
    got = {k: {(ds, t) for ds, t, kk in drawn if kk == k} for k in ("failed", "noteval")}
    assert got == E_MISSING, got
    # super-header over the three genome blocks, styled as B's class-system headers: says these genomes are real
    xa, xb = G["XD0"], G["XD0"] + 3 * (3 * w + gap) - gap
    can.text((xa + xb) / 2, R0 - 11.6, "real genomes", ha="center", va="center", fontsize=FL, color=MUTED)
    can.plot([xa + 0.5, xb - 0.5], [R0 - 9.6, R0 - 9.6], color=INK, lw=0.6)


# ============================================================================ legend (claims asserted)
def sim_arrays():
    """the legend says simulated genomes lack tandem arrays: count features of each S2 genome's PrinTE BED whose name
    or class mentions a satellite, tandem array, rDNA or simple repeat (0-based BED; names only are read)"""
    out = []
    for s in ("S2_1", "S2_2", "S2_3"):
        n = k = 0
        with open(B / "sims" / s / "gen1000000_final.bed") as f:
            for line in f:
                n += 1
                nm = line.split("\t")[3].lower() if line.count("\t") >= 3 else ""
                k += any(w in nm for w in ("satel", "tandem", "rdna", "simple"))
        out.append(dict(genome=s, features=n, satellite_tandem_rdna_simple_matches=k))
    return {"legend_sim_features": pd.DataFrame(out)}


def legend_text(ex, Bd, Cd, Cp, Dd):
    # A claims
    found_o = [NAME[t] for t in ex["found_by"] if t != "ltr_peel"]
    miss_t = [NAME[t] for t in TOOLS if t not in ex["tfound_by"]]
    assert len(found_o) == 3 and len(miss_t) == 1 and len(ex["nfound_by"]) == len(TOOLS)
    num("legend", "A: tools other than LTR_peel that rebuild the host / tool missing the tenant",
        f"{', '.join(found_o)} / {miss_t[0]}", "fig_benchmark_numbers.tsv A 'tools matching host / tenant / "
        "non-nested element' (score_elements.match RO 0.9)")
    # B claims
    V = Bd["V"]
    val = lambda c, m, ts: [V[(c, m, t)] for t in ts if V[(c, m, t)] is not None]     # noqa: E731
    for _, c, _ in BCLS:
        for m in ("recall", "f1"):
            assert V[(c, m, "ltr_peel")] > max(val(c, m, TOOLS[1:])), (c, m)
        assert V[(c, "precision", "ltr_peel")] < max(val(c, "precision", TOOLS[1:])), c
        assert V[(c, "precision", "hite")] >= max(val(c, "precision", TOOLS)), c
    nocall = sorted({(c, t) for (c, m, t), v in V.items() if v is None})
    assert nocall == [("proteinless", "dante_ltr")], nocall
    sd = Bd["sd_max"]
    sd_o = max(v for t, v in sd.items() if t != "ltrpred")
    assert sd_o <= 0.023 and sd["ltrpred"] == max(sd.values()) and sd["ltrpred"] > 0.05, sd
    num("legend", "B: largest across-genome SD, every tool but LTRpred / LTRpred (drawn classes, P/R/F1)",
        f"{sd_o:.4f} / {sd['ltrpred']:.4f}", "fig_benchmark_numbers.tsv legend row of the same name "
        "(experiments/figpub/sim_class_prf_joint_summary.tsv *_sd)")
    num("legend", "B: LTRpred's S2_3 run incomplete (no chr6 calls, chr5 partly; gt ltrdigest HMMER error, recorded ok)",
        "one run of 60", "fig_benchmark_numbers.tsv legend row of the same name (experiments/figpub/memo.md section e)")
    # C claims (coding, pooled recall-halved points) and the protein-less exception
    CS = Cd["S"]
    d50c = {t: float(CS[t]["d50"]) for t in TOOLS if CS[t]["d50"][0].isdigit()}
    seven = [t for t in TOOLS if t not in ("ltr_peel", "dante_ltr")]
    lo7, hi7 = min(d50c[t] for t in seven), max(d50c[t] for t in seven)
    q50 = d50c["ltr_peel"]
    assert q50 > 0.3 and hi7 < 0.2 and CS["dante_ltr"]["d50"].startswith(">"), (q50, hi7)
    for t in seven:
        if FLOOR[t] is not None:
            assert d50c[t] < k2p_of_identity(FLOOR[t]), t
    num("legend", "C: coding pooled d50 range of the seven other tools (K2P)", f"{lo7:.4f}-{hi7:.4f}",
        "fig_benchmark_numbers.tsv legend row of the same name (computed, C dots)")
    num("legend", "C: LTR_peel coding pooled d50 (K2P)", f"{q50:.4f}", "fig_benchmark_numbers.tsv C 'LTR_peel coding "
        "pooled d50 (K2P at which pooled recall halves; dot)'")
    PS = Cp["S"]
    eq, eh, dh = float(PS["ltr_peel"]["d50"]), float(PS["hirepeat"]["d50"]), d50c["hirepeat"]
    assert eq < q50 and eh > eq and eh > dh
    num("legend", "C/D: pooled d50 LTR_peel coding / protein-less; HiRepEat coding / protein-less",
        f"{q50:.4f} / {eq:.4f}; {dh:.4f} / {eh:.4f}", "fig_benchmark_numbers.tsv legend row of the same name "
        "(computed from experiments/nonauto/div_elements.tsv; protein-less = published panel D)")
    # D claims
    mk = Dd["mk"]
    rec = [v["recall"] for v in mk.values()]
    pq = [mk[(ds, "ltr_peel")]["precision"] for ds in E_DS]
    strict = ("ltr_retriever", "hite", "dante_ltr")
    ps = [v["precision"] for (d, t), v in mk.items() if t in strict]
    assert min(rec) >= 0.75 and max(pq) < min(ps), (min(rec), pq, ps)
    for ds in E_DS:                                 # HiTE has the best printed F1 in every genome; rice tie
        f1 = {t: float(f"{v['f1']:.2f}") for (d, t), v in mk.items() if d == ds}
        assert f1["hite"] == max(f1.values()), ds
        assert mk[(ds, "ltr_peel")]["f1"] < max(v["f1"] for (d, t), v in mk.items() if d == ds), ds
        tie = sorted(t for t, v in f1.items() if v == f1["hite"])
        assert tie == (["hite", "ltr_retriever"] if ds == "rice" else ["hite"]), (ds, tie)
    st = {k: v[1] for k, v in Dd["reasons"].items()}
    assert st[("finch", "ltr_retriever")] == "resources.tsv status=timeout"
    assert st[("finch", "dante_ltr")] == "resources.tsv status=failed"
    assert st[("finch", "hirepeat")].startswith("masking.tsv ci: incomplete: masking projected")
    num("legend", "D: recall range over the drawn cells (3 genomes)", f"{min(rec):.4f}-{max(rec):.4f}",
        "fig_benchmark_numbers.tsv legend 'E: recall range / precision range over the drawn cells (3 genomes)' "
        "(report/tables/masking.tsv sens)")
    num("legend", "D: LTR_peel precision range (A. thaliana, rice, zebra finch)", f"{min(pq):.4f}-{max(pq):.4f}",
        "report/tables/masking.tsv prec (= fig_benchmark_numbers.tsv E rows)")
    num("legend", "D: LTR_retriever, HiTE, DANTE_LTR precision range (finished genomes)",
        f"{min(ps):.4f}-{max(ps):.4f}", "report/tables/masking.tsv prec (= fig_benchmark_numbers.tsv E rows)")
    num("legend", "D: best F1 per genome (printed, 2 decimals): A. thaliana / rice / zebra finch",
        " / ".join(" and ".join(NAME[t] for t in TOOLS if (ds, t) in mk and
                                f"{mk[(ds, t)]['f1']:.2f}" == f"{mk[(ds, 'hite')]['f1']:.2f}")
                   + f" {mk[(ds, 'hite')]['f1']:.2f}" for ds in E_DS),
        "report/tables/masking.tsv F1 (asserted: HiTE best in each genome, tied with LTR_retriever in rice)")
    num("legend", "D: missing cells, zebra finch (reason as recorded)",
        "; ".join(f"{NAME[t]} {st[('finch', t)]}" for t in ("ltr_retriever", "dante_ltr", "hirepeat")),
        "report/tables/resources.tsv status (36 h cap: bench memo.md); report/tables/masking.tsv ci")
    arr = {r.genome: (int(r.features), int(r.satellite_tandem_rdna_simple_matches))
           for r in read_sd("legend_sim_features").itertuples()}
    assert all(k == 0 for _, k in arr.values()), arr
    num("legend", "B: simulated genomes hold no satellite, tandem-array, rDNA or simple-repeat feature (S2_1, S2_2, "
        "S2_3: features, matches)", "; ".join(f"{s} {n:,}, {k}" for s, (n, k) in arr.items()),
        "sims/<S2_x>/gen1000000_final.bed (this script, case-insensitive satel|tandem|rDNA|simple)")

    txt = (
        f"**LTR_peel finds nested and old LTR retrotransposons that most callers miss, while strict callers are more "
        f"precise and build cleaner repeat libraries.** Panels A–C use genomes simulated with PrinTE, which shares an "
        f"author with LTR_peel, and LTR_peel was tuned on PrinTE simulations and on rice.\n\n"
        f"**A**, A host, its tenant and a non-nested element in a simulated genome, the axis broken between them. A "
        f"host is an element split in two by a later insertion, its tenant. Each tool's calls follow to scale with "
        f"the features it reports (TSDs enlarged). LTR_peel, {found_o[0]}, {found_o[1]} and {found_o[2]} rebuild the "
        f"host, all but {miss_t[0]} the tenant, and every tool the non-nested element.\n\n"
        f"**B**, Precision (prec.), recall (rec.) and F1 by class, averaged over three deeply nested simulated genomes "
        f"(n, true elements). Coding elements carry at least two kinds of LTR retrotransposon protein domain from one "
        f"superfamily, protein-less elements none. Bold marks each column's best and a dash a class the caller never "
        f"calls. LTR_peel leads on recall and F1 in every class, HiTE on precision. The genomes agree (SD at most "
        f"{sd_o:.2f}) except for LTRpred, one of whose runs lost calls to a pipeline error. Simulated genomes lack "
        f"tandem arrays and tandem duplicated genes, and their truth is complete, unlike curated libraries, so "
        f"simulated precision is optimistic.\n\n"
        f"**C**, Recall of coding elements as one simulated burst's LTRs diverge, three genomes pooled (bands, 95% "
        f"intervals over families). Seven tools lose half (dots) by a divergence of {lo7:.2f}–{hi7:.2f}, those with "
        f"a documented identity floor (dashed) well before it. LTR_peel keeps half to {q50:.2f}, and DANTE_LTR, which "
        f"reads proteins, never halves. Protein-less elements of the same series are harder for LTR_peel, half lost "
        f"by {eq:.2f} against {eh:.2f} for HiRepEat. The aging adds substitutions only, so real old copies are likely "
        f"harder for every tool.\n\n"
        f"**D**, Masking of three real genomes by each tool's library, against curated annotation. Precision is the "
        f"share of masked sequence that curators annotate as LTR retrotransposon, recall the share of their LTR "
        f"retrotransposon annotation that is masked. Recall is high for every tool ({min(rec):.2f}–{max(rec):.2f}), so "
        f"tools differ in precision, {min(pq):.2f}–{max(pq):.2f} for LTR_peel against {min(ps):.2f}–{max(ps):.2f} for "
        f"LTR_retriever, HiTE and DANTE_LTR. HiTE has the best F1 in each genome, tied with LTR_retriever in rice. "
        f"Dark gray marks tools that failed on zebra finch (no result within 36 h, or a crash without calls), and "
        f"hatching a library whose masking would have run far past 36 h.")
    for bad in (":", ";", "—"):                    # house style: no colons, semicolons or dash asides
        assert bad not in txt, bad
    return txt


# ============================================================================ checks
def cross_check(path):
    """every recorded A-D number equals the published figure's number for the same quantity (fig_benchmark_numbers.tsv,
    panel E there = D here); returns the count checked"""
    old = pd.read_csv(OLD, sep="\t", dtype=str, keep_default_na=False)
    oldv = {(r.panel, r.quantity): r.value for r in old.itertuples()}
    new = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    n = 0
    for r in new.itertuples():
        if r.panel in OLDPANEL:
            k = (OLDPANEL[r.panel], r.quantity)
            if k not in oldv:
                raise SystemExit(f"{STEM}: {r.panel} '{r.quantity}' not in fig_benchmark_numbers.tsv")
            if oldv[k] != r.value:
                raise SystemExit(f"{STEM}: {r.panel} '{r.quantity}' = {r.value}, published {oldv[k]}")
            n += 1
        elif r.panel == "legend" and (r.panel, r.quantity) in oldv and oldv[(r.panel, r.quantity)] != r.value:
            raise SystemExit(f"{STEM}: legend '{r.quantity}' = {r.value}, published {oldv[(r.panel, r.quantity)]}")
    # protein-less d50 used in the legend = published panel D dots
    for t in ("ltr_peel", "hirepeat"):
        q = f"{NAME[t]} protein-less pooled d50 (K2P at which pooled recall halves; dot)"
        assert oldv[("D", q)] == PL_D50[t], (t, oldv[("D", q)], PL_D50[t])
    return n


def honesty_plot(Bd, Dd):
    """plain bars of the numbers behind B and D (TMP), to compare with the scorecards"""
    fig, axs = plt.subplots(2, 5, figsize=(180 * MM, 90 * MM))
    xs = np.arange(len(TOOLS))
    for j, (_, c, lab) in enumerate(BCLS):
        ax = axs[0, j]
        for i, (m, _) in enumerate(CMET):
            ax.bar(xs + (i - 1) * 0.27, [Bd["V"][(c, m, t)] or 0 for t in TOOLS], 0.27, label=m,
                   color=["#999999", "#555555", "#111111"][i])
        ax.set_title(f"B {lab}", fontsize=6)
    for j, ds in enumerate(E_DS):
        ax = axs[1, j]
        for i, (_, key, _) in enumerate(EMET):
            ax.bar(xs + (i - 1) * 0.27, [Dd["mk"].get((ds, t), {}).get(key, 0) for t in TOOLS], 0.27,
                   color=["#999999", "#555555", "#111111"][i], label=key)
        ax.set_title(f"D {ds} (P, R, F1)", fontsize=6)
    for ax in axs.flat:
        ax.set_xticks(xs)
        ax.set_xticklabels([NAME[t] for t in TOOLS], rotation=90, fontsize=4.5)
        ax.set_ylim(0, 1)
        ax.tick_params(labelsize=5)
    axs[1, 3].axis("off")
    axs[1, 4].axis("off")
    axs[0, 0].legend(fontsize=5)
    fig.tight_layout()
    fig.savefig(TMP / "honesty_plain_bars.png", dpi=200)
    plt.close(fig)


PL_D50 = {}


# ============================================================================ main
def build_source_data():
    """full mode: load and check the benchmark data, write the drawn tables to SD"""
    bench_root()
    TMP.mkdir(parents=True, exist_ok=True)
    T, scal = load_example()
    log("A: example loaded and checked")
    T.update(load_B())
    curve, d50, n_fam, n_el = load_curve("coding")
    _, d50p, _, _ = load_curve("proteinless", tools=("ltr_peel", "hirepeat"))
    T["C_curves"] = pd.DataFrame(curve)
    T["C_d50"] = pd.DataFrame(d50 + d50p)
    scal += [dict(panel="C", name="coding_element_families", value=n_fam),
             dict(panel="C", name="coding_elements_per_timepoint", value=n_el)]
    T.update(load_D())
    T.update(sim_arrays())
    T["scalars"] = pd.DataFrame([dict(x, value=str(x["value"])) for x in scal])   # ints stay ints, floats repr
    log("B-D: tables loaded and checked")
    write_source_data(T)


def main(argv=None):
    global VERBOSE
    ap = argparse.ArgumentParser(description="Build benchmark_summary (LTR_peel README benchmark figure): PDF (figure "
                                             "+ legend), 300-dpi PNG, numbers TSV, 880-px check copy in tmp/.")
    ap.add_argument("-v", "--verbose", action="store_true", help="progress and checks")
    ap.add_argument("--from-source-data", action="store_true",
                    help="redraw from source_data/benchmark_summary/ only (no BENCH_ROOT needed)")
    ap.add_argument("--honesty", action="store_true", help="also write plain bar charts of B and D to tmp/")
    ap.add_argument("--draft", default=None, help="write only a PNG draft to this name inside tmp/")
    a = ap.parse_args(argv)
    VERBOSE = a.verbose
    TMP.mkdir(parents=True, exist_ok=True)
    log("start", True)
    if not a.from_source_data:
        build_source_data()
    ex = example_from_tables()
    Bd = B_from_tables()
    Cd, Cp = C_from_tables()
    PL_D50.update({t: Cp["S"][t]["d50"] for t in ("ltr_peel", "hirepeat")})
    Dd = D_from_tables()
    log(f"tables read from {SD}")
    if a.honesty:
        honesty_plot(Bd, Dd)
        log(f"honesty plot: {TMP / 'honesty_plain_bars.png'}", True)

    probe = plt.figure(figsize=(W * MM, 100 * MM), dpi=300)
    layout(probe)
    plt.close(probe)
    fig = plt.figure(figsize=(W * MM, G["HFIG"] * MM), dpi=300)      # text measured at the output resolution
    can = poafig.key_axes(fig, 0, 0, W, G["HFIG"])
    # band 1: A
    names(can, G["XN"], [row_y(k) for k in range(len(TOOLS))])
    panel_A(fig, can, ex)
    # band 2: B and the shared score scale
    panel_B(fig, can, Bd)
    colorbar(can)
    # band 3: C (key above its plot) and D
    key_C(fig, can)
    panel_C(fig, can, Cd)
    panel_D(fig, can, Dd)
    for x, y, L in ((0, 0.2, "A"), (0, G["Y2"] - 0.2, "B"), (0, G["Y3"] - 0.2, "C"), (G["XDL"], G["Y3"] - 0.2, "D")):
        letter(fig, x, y, L)
    if a.draft:
        p = TMP / Path(a.draft).name
        fig.savefig(p, dpi=300)
        log(f"draft {p} ({W:.0f} x {G['HFIG']:.1f} mm)", True)
        return
    txt = legend_text(ex, Bd, Cd, Cp, Dd)
    out = poafig.write_figure(STEM, fig, txt, None, HERE, str(SD if a.from_source_data else B),
                              str(HERE / f"{STEM}.py"))
    PNG.parent.mkdir(parents=True, exist_ok=True)
    os.replace(HERE / f"{STEM}.png", PNG)           # the PNG is the README's figure
    tsv = HERE / f"{STEM}_numbers.tsv"
    pd.DataFrame(NUM).to_csv(tsv, sep="\t", index=False)
    n = cross_check(tsv)
    df = pd.read_csv(tsv, sep="\t", dtype=str, keep_default_na=False)   # the check passed: say so per row
    df["published"] = [f"= fig_benchmark_numbers.tsv panel {OLDPANEL[p_]}, same quantity and value (checked)"
                       if p_ in OLDPANEL else "" for p_ in df.panel]
    df.to_csv(tsv, sep="\t", index=False)
    from PIL import Image                          # README check: the PNG as GitHub shows it (880 px wide)
    im = Image.open(PNG)
    im.resize((880, round(im.height * 880 / im.width)), Image.LANCZOS).save(TMP / f"{STEM}_880px.png")
    log(f"done {PNG}, {out['pdf']} fig {out['fig_mm']} legend {out['short_mm']} mm, {out['short_words']} words "
        f"{out['short_par_words']}; {len(NUM)} numbers, {n} equal to fig_benchmark_numbers.tsv", True)


if __name__ == "__main__":
    main()
