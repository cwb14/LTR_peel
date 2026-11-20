#!/usr/bin/env python3
"""Fig. S2 figS2_masking: curated-library masking on five real genomes, details.

Panels
  A  rice, one typical 200-kb window: genes, curated LTR-RT and other curated TE hits, then per tool its
     intact calls (LTRs drawn; filled = library-supported, outlined = not) and the bases its library masks (gray =
     curated LTR-RT, tool color = outside curated LTR-RT, open = curated LTR-RT left unmasked).
  B  per genome and tool: curated LTR-RT bases masked (left, within the frame of all curated LTR-RT) and
     bases masked outside curated LTR-RTs (right), Mb.
  C  lib-test masking metrics per tool x genome (precision, recall = sensitivity, F1, specificity, accuracy;
     eval/lib/<ds>/<tool>.intact.tsv; F1 = report/tables/masking.tsv F1); two decimals as Fig. 1; viridis 0-1, column
     best bold (>= 2 tools); missing cells "tool failed" / "not evaluated".
Fig. 1E shows precision, recall and F1 for A. thaliana, rice and zebra finch
(maize and gorilla F1 only here), Fig. 1F curated-family recovery.
Inputs: report/tables/{masking,resources}.tsv, eval/lib/<ds>.ci.tsv, eval/lib/<ds>/<tool>.intact.tsv; raw:
  mask/rice/*/rice.fa.out, calls/rice/<tool>.tsv, eval/adj/rice/<tool>.calls.tsv, data/rice/genes/all.gff3
  (window ranking via figS2_prep.py, cached in tmp/).
Outputs: figS2_masking.{pdf,png}, figS2_masking_numbers.tsv (every number on the figure and legend, with its source),
  tmp/f1_crosscheck.tsv (eval/lib F1, precision, sensitivity vs masking.tsv, every scored cell),
  source_data/figS2_masking/*.tsv (the tables the figure is drawn from; full mode writes them, reads them back and
  draws only from them; --from-source-data redraws from them without BENCH_ROOT, skipping the raw-data checks).
Usage: figS2_masking.py [--reprep] [--from-source-data] [-v]
"""
import argparse
import glob
import os
import sys
from pathlib import Path

sys.dont_write_bytecode = True                       # no __pycache__ next to the read-only imported modules

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # poafig.py sits beside this script
import figS2_prep as P  # noqa: E402
import poafig  # noqa: E402
import matplotlib as mpl  # noqa: E402
import matplotlib.font_manager as fm  # noqa: E402

# the shared matplotlib font cache on this host lists no system fonts: register Liberation in-process (cache untouched)
for _f in sorted(glob.glob("/usr/share/fonts/truetype/liberation2/LiberationSans-*.ttf")
                 + glob.glob("/usr/share/fonts/truetype/liberation2/LiberationMono-*.ttf")):
    fm.fontManager.addfont(_f)
poafig.setup_mpl()
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
from matplotlib.colors import to_rgb  # noqa: E402
from matplotlib.patches import Polygon, Rectangle  # noqa: E402

STEM = "figS2_masking"
SD = HERE / "source_data" / STEM
B = T = None                                         # $BENCH_ROOT and its report/tables, full mode only
MM = poafig.MM
INK = poafig.INK
DS = P.DATASETS
TOOLS = P.TOOLS
LABEL = {"ltr_peel": "LTR_peel", "ltr_retriever": "LTR_retriever", "ltrharvest": "LTRharvest",
         "ltr_finder": "LTR_Finder", "look4ltrs": "Look4LTRs", "hirepeat": "HiRepEat", "hite": "HiTE",
         "ltrpred": "LTRpred", "dante_ltr": "DANTE_LTR"}
COLOR = {"ltr_peel": "#D55E00", "ltr_retriever": "#0072B2", "ltrharvest": "#56B4E9", "ltr_finder": "#009E73",
         "look4ltrs": "#CC79A7", "hirepeat": "#E69F00", "hite": "#000000", "ltrpred": "#A6761D",
         "dante_ltr": "#6A3D9A"}
DSNAME = {"arab": r"$\it{A.\ thaliana}$", "rice": "rice", "finch": "zebra finch", "maize": "maize",
          "gorilla": "gorilla"}
CUR_LTR = "#3C3C3C"       # curated LTR-RT hits (dark gray)
CUR_OTH = "#CFCFCF"       # other curated TE hits (context)
GENE = "#5A5A5A"
MASK_TP = "#A3A3A3"       # masked and curated LTR-RT
MISSC = "#F2F2F2"         # missing cells, not evaluated (hatched)
TOOLFAILC = "#5A5A5A"     # missing cells, tool failed (dark gray, white text)
HATCH = "//////"
HATCHC = "#B8B8B8"
MUTED = "#6B6B6B"
# score map (higher = better), figure-set rule: viridis on a fixed 0-1 scale, dark purple (0) to yellow (1); perceptually
# uniform, color-vision-deficiency safe, with more lightness and hue change across 0.6-1 than cividis
CMAP = mpl.colormaps["viridis"]
WIN_CHROM_NAME = {"c0011": "Chr11"}
# lib-test column, header; order precision, recall, F1, specificity, accuracy
METRICS = [("prec", "prec."), ("sens", "recall"), ("F1", "F1"), ("spec", "spec."), ("accu", "acc.")]
MNAME = {"prec": "precision", "sens": "recall (sensitivity)", "F1": "F1", "spec": "specificity", "accu": "accuracy"}


def tint(c, f):
    """mix color c with white; f = share of white"""
    r, g, b = to_rgb(c)
    return (r + (1 - r) * f, g + (1 - g) * f, b + (1 - b) * f)


def log(msg, verbose=True):
    if verbose:
        print(f"[{STEM}] {msg}", file=sys.stderr, flush=True)


# ================================================================== numbers
NUM = []


def num(panel, ds, tool, quantity, value, source):
    NUM.append(dict(panel=panel, dataset=ds, tool=tool, quantity=quantity, value=value, source=source))
    return value


TOOL_FAIL, NOT_EVAL = "tool failed", "not evaluated"      # the only two missing-cell labels (figure-set rule)
TOOL_FAIL_STATUS = {"timeout": "did not finish within 36 h", "mem_guard": "exceeded 200 GB of memory",
                    "failed": "crashed with no calls"}
# expected cells of each class (the expected list); the script asserts the data give exactly these
EXPECT_TOOL_FAIL = {("finch", "ltr_retriever"), ("maize", "ltr_retriever"), ("gorilla", "ltr_retriever"),
                    ("gorilla", "ltrharvest"), ("gorilla", "ltrpred"), ("gorilla", "look4ltrs"),
                    ("finch", "dante_ltr"), ("gorilla", "dante_ltr")}
EXPECT_NOT_EVAL = {("finch", "hirepeat"), ("gorilla", "ltr_peel"), ("gorilla", "ltr_finder"), ("gorilla", "hirepeat"),
                   ("maize", "ltr_peel"), ("maize", "ltrharvest"), ("maize", "ltr_finder"), ("maize", "look4ltrs"),
                   ("maize", "hirepeat"), ("maize", "hite"), ("maize", "dante_ltr")}


def missing_reasons():
    """(ds, tool) -> (label, source) for every tool x genome with a missing result. Two classes only:
    'tool failed' = resources.tsv status timeout / mem_guard / failed (the tool's own run); 'not evaluated' = the tool
    finished but masking.tsv marks its intact-library masking (or library build) incomplete ('incomplete: ...' ci,
    scripts/summarize.py INCOMPLETE). Exits unless every candidate falls in exactly one class."""
    res = pd.read_csv(T / "resources.tsv", sep="\t")
    m = pd.read_csv(T / "masking.tsv", sep="\t")
    out = {}
    for ds in DS:
        for t in TOOLS:
            r = res[(res.ds == ds) & (res.tool == t)]
            st = r.status.iloc[0] if len(r) else "?"
            row = m[(m.ds == ds) & (m.lib == f"{t}.intact")]
            fail = st in TOOL_FAIL_STATUS
            noeval = bool(len(row)) and row.F1.iloc[0] == "." and str(row.ci.iloc[0]).startswith("incomplete")
            if fail and noeval:
                sys.exit(f"error: {ds} {t} is both 'tool failed' and 'not evaluated'")
            if fail:
                out[ds, t] = (TOOL_FAIL, f"resources.tsv status={st} ({TOOL_FAIL_STATUS[st]})")
            elif noeval:
                out[ds, t] = (NOT_EVAL, f"masking.tsv ci: {row.ci.iloc[0]}")
            elif st != "ok":
                sys.exit(f"error: {ds} {t}: unclassified status {st}")
    got_f = {k for k, v in out.items() if v[0] == TOOL_FAIL}
    got_n = {k for k, v in out.items() if v[0] == NOT_EVAL}
    if got_f != EXPECT_TOOL_FAIL or got_n != EXPECT_NOT_EVAL:
        sys.exit(f"error: missing-cell classes differ from the expected lists: tool failed {sorted(got_f ^ EXPECT_TOOL_FAIL)},"
                 f" not evaluated {sorted(got_n ^ EXPECT_NOT_EVAL)}")
    return out


def load_tables():
    m = pd.read_csv(T / "masking.tsv", sep="\t")
    m = m[(m.role == "primary") & (m.F1 != ".")].copy()
    m["tool"] = m.lib.str.replace(".intact", "", regex=False)
    m = m[m.tool.isin(TOOLS)]
    for c in ("F1", "F1_lo", "F1_hi", "prec", "sens", "FP_mb", "FN_mb"):
        m[c] = m[c].astype(float)
    ci = pd.concat([pd.read_csv(B / f"eval/lib/{ds}.ci.tsv", sep="\t") for ds in DS])
    ci["tool"] = ci.lib.str.replace(".intact", "", regex=False)
    m = m.merge(ci[["ds", "lib", "TP", "FP", "FN"]], on=["ds", "lib"], how="left")
    bad = m[(abs(m.FP / 1e6 - m.FP_mb) > 0.006) | (abs(m.FN / 1e6 - m.FN_mb) > 0.006)]
    if len(bad):
        sys.exit(f"error: eval/lib/*.ci.tsv and masking.tsv disagree:\n{bad}")
    return m


def load_libtest(M, verbose):
    """lib-test metrics of every primary intact-call library (eval/lib/<ds>/<tool>.intact.tsv; files repeat their
    header line, data rows are those whose label is not 'label', one per file). Exits unless the scored set equals
    masking.tsv's scored primary rows and F1, precision and sensitivity (rounded to 4 decimals, as masking.tsv) agree
    for every cell."""
    rows = []
    for ds in DS:
        for t in TOOLS:
            p = B / f"eval/lib/{ds}/{t}.intact.tsv"
            if not p.exists() or p.stat().st_size == 0:
                continue
            d = pd.read_csv(p, sep="\t", dtype=str)
            d = d[d.label != "label"]
            if len(d) != 1:
                sys.exit(f"error: {p}: {len(d)} data rows (expected 1)")
            r = d.iloc[0]
            rows.append(dict(ds=ds, tool=t, **{k: float(r[k]) for k in ("sens", "spec", "accu", "prec", "F1")},
                             **{k: int(r[k]) for k in ("TP", "TN", "FP", "FN")}))
    L = pd.DataFrame(rows)
    got, want = set(zip(L.ds, L.tool)), set(zip(M.ds, M.tool))
    if got != want:
        sys.exit(f"error: eval/lib intact libraries {sorted(got ^ want)} differ from masking.tsv's scored rows")
    x = L.merge(M[["ds", "tool", "F1", "prec", "sens"]], on=["ds", "tool"], suffixes=("", "_masking"))
    for c in ("F1", "prec", "sens"):
        x[f"{c}_equal"] = (x[c].round(4) - x[f"{c}_masking"]).abs() < 1e-9
    x.to_csv(HERE / "tmp/f1_crosscheck.tsv", sep="\t", index=False)
    if not x[["F1_equal", "prec_equal", "sens_equal"]].all().all():
        sys.exit(f"error: eval/lib metrics differ from masking.tsv (tmp/f1_crosscheck.tsv):\n{x}")
    # definitions check: lib-test's specificity and accuracy from its own counts
    tot = L.TP + L.TN + L.FP + L.FN
    if ((L.spec - L.TN / (L.TN + L.FP)).abs().max() > 1e-9) or ((L.accu - (L.TP + L.TN) / tot).abs().max() > 1e-9):
        sys.exit("error: lib-test spec/accu do not equal TN/(TN+FP) and (TP+TN)/total")
    L["total"] = tot
    log(f"F1 cross-check: {len(x)} of {len(x)} cells equal masking.tsv (F1, precision, sensitivity)", verbose)
    return L, len(x)


# ================================================================== panel A data
def window_choice(verbose):
    """first window in the ranked list (figS2_prep.windows) passing the checks: LTR_peel's library masks >= 99% of
    the window's curated LTR-RT bases (no visible LTR_peel miss), and every drawn call's library support
    recomputed from the curated .out (>= 80% of the span) agrees with eval/adj"""
    w = P.windows()
    ref = P.libtest_ci.ltr_hits(str(B / "mask/rice/curated/rice.fa.out"))
    lq = P.libtest_ci.ltr_hits(str(B / "mask/rice/ltr_peel.intact/rice.fa.out"))
    log_rows = []
    for i, r in w.head(30).iterrows():
        c, s, e = r.chrom, int(r.start), int(r.end)
        rm = clip(ref.get(c, ()), s, e)
        tm = clip(lq.get(c, ()), s, e)
        cur = span(rm)
        miss = cur - span(P.intersect(rm, tm))
        ok_sens = miss <= 0.01 * cur
        log_rows.append((i + 1, c, s, e, r.n_calls, r.n_unsup, r.n_genes, round(r.lq_prec, 3), round(r.dist, 3),
                         cur, miss, ok_sens))
        if ok_sens:
            chosen = (c, s, e, i + 1)
            break
    cand = pd.DataFrame(log_rows, columns=["rank", "chrom", "start", "end", "n_calls", "n_unsup", "n_genes",
                                           "lq_prec", "dist", "curated_ltr_bp", "lq_unmasked_bp", "pass_sens"])
    cand.to_csv(HERE / "tmp/rice_window_candidates.tsv", sep="\t", index=False)
    log(f"window: {chosen} (candidates in tmp/rice_window_candidates.tsv)", verbose)
    return chosen, w


def clip(iv, s, e):
    return P.merge([(max(a, s), min(b, e)) for a, b in iv if b > s and a < e])


def span(iv):
    return int((iv[:, 1] - iv[:, 0]).sum()) if len(iv) else 0


def subtract(a, b):
    """a minus b (merged interval arrays)"""
    out = []
    j = 0
    for s, e in a:
        cur = s
        while j < len(b) and b[j, 1] <= cur:
            j += 1
        k = j
        while k < len(b) and b[k, 0] < e:
            if b[k, 0] > cur:
                out.append((cur, b[k, 0]))
            cur = max(cur, b[k, 1])
            k += 1
        if cur < e:
            out.append((cur, e))
    return np.array(out, dtype=np.int64).reshape(-1, 2)


def window_data(c, s, e):
    ren = P.ren_map("rice")
    back = {v: k for k, v in ren.items()}
    D = {}
    g = P.genes("rice")
    D["genes"] = g[(g.chrom == c) & (g.end > s) & (g.start < e)].sort_values("start")
    hits = []
    for t in P.evidence.rm_hits(str(B / "mask/rice/curated/rice.fa.out")):
        if t[4] == c and int(t[6]) > s and int(t[5]) - 1 < e:
            ltr = t[10] != "Unspecified" and P.libtest_ci.LTR_RE.search(t[10]) is not None
            hits.append((int(t[5]) - 1, int(t[6]), t[10], t[9], ltr))
    D["cur"] = pd.DataFrame(hits, columns=["start", "end", "cls", "fam", "ltr"])
    ref = P.libtest_ci.ltr_hits(str(B / "mask/rice/curated/rice.fa.out"))
    refw = clip(ref.get(c, ()), s, e)
    refall = P.merge(ref.get(c, ()))
    D["cur_ltr_bp"] = span(refw)
    calls, masks, checks = {}, {}, []
    for t in TOOLS:
        k = pd.read_csv(B / f"calls/rice/{t}.tsv", sep="\t", dtype={"lltr_end": str, "rltr_start": str})
        k = k[(k.chrom == back[c]) & (k.end > s) & (k.start < e)].copy()
        a = pd.read_csv(B / f"eval/adj/rice/{t}.calls.tsv", sep="\t",
                        usecols=["chrom", "start", "end", "supported", "tgca", "tsd", "k2p", "coding", "gene"])
        a = a[a.chrom == c].rename(columns={"tsd": "tsd_exact"})
        k = k.merge(a.drop(columns="chrom"), on=["start", "end"], how="left")
        if k.supported.isna().any():
            sys.exit(f"error: rice {t}: call without an eval/adj row in the window")
        for _, r in k.iterrows():           # recompute library support from the curated .out (proxies.py rule)
            cov = span(P.intersect(clip(refall, r.start, r.end), P.merge([(r.start, r.end)]))) / (r.end - r.start)
            checks.append((t, r.id, r.start, r.end, int(r.supported), round(cov, 3), int(cov >= 0.8),
                           int(r.tgca), int(r.tsd_exact), r.k2p, int(r.coding), int(r.gene)))
        calls[t] = k
        tst = P.libtest_ci.ltr_hits(str(B / f"mask/rice/{t}.intact/rice.fa.out"))
        mw = clip(tst.get(c, ()), s, e)
        masks[t] = dict(tp=P.intersect(refw, mw), fp=subtract(mw, refw), fn=subtract(refw, mw), all=mw)
    chk = pd.DataFrame(checks, columns=["tool", "id", "start", "end", "supported_adj", "curated_ltr_cov",
                                        "supported_recomputed", "tgca", "tsd", "k2p", "coding", "gene_cds"])
    chk.to_csv(HERE / "tmp/rice_window_call_checks.tsv", sep="\t", index=False)
    if (chk.supported_adj != chk.supported_recomputed).any():
        sys.exit("error: drawn calls' library support differs from the curated .out (tmp/rice_window_call_checks.tsv)")
    D["calls"], D["masks"], D["checks"] = calls, masks, chk
    return D


# ================================================================== drawing helpers
def mm_ax(fig, x, y, w, h):
    return poafig.mm_axes(fig, x, y, w, h)


def track_ax(fig, x, y, w, h, xlim):
    ax = mm_ax(fig, x, y, w, h)
    ax.set_xlim(*xlim)
    ax.set_ylim(h, 0)
    ax.axis("off")
    return ax


def rel_lum(c):
    r, g, b = (x / 12.92 if x <= 0.04045 else ((x + 0.055) / 1.055) ** 2.4 for x in to_rgb(c))
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def text_color(v):
    """white on dark cells, ink on light ones: whichever has the higher WCAG contrast with the cell's relative
    luminance"""
    L = rel_lum(CMAP(v))
    return "white" if (1.05) / (L + 0.05) > (L + 0.05) / (rel_lum(INK) + 0.05) else INK


def missing_box(ax, x, y, w, h, text, fs=5.5):
    """missing result, two styles only (figure-set rule): 'tool failed' = dark gray fill, white text; 'not
    evaluated' = light hatch, gray text"""
    if text == TOOL_FAIL:
        ax.add_patch(Rectangle((x, y), w, h, fc=TOOLFAILC, ec="none", lw=0))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs, color="white")
    elif text == NOT_EVAL:
        ax.add_patch(Rectangle((x, y), w, h, fc=MISSC, ec="none", lw=0, hatch=HATCH, hatchcolor=HATCHC))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs, color="#4D4D4D",
                bbox=dict(boxstyle="square,pad=0.06", fc=MISSC, ec="none"))
    else:
        sys.exit(f"error: unknown missing-cell label {text!r}")


def label_column(fig, x0, y, H, rows, RH, extra=None):
    """tool names (LTR_peel bold) with a square of the tool color, right-aligned to x0; extra: [(label, y)]"""
    lab = mm_ax(fig, 0, y, x0 - 0.8, H)
    lab.set_xlim(0, x0 - 0.8)
    lab.set_ylim(H, 0)
    lab.axis("off")
    for i, t in enumerate(rows):
        yc = i * RH + RH / 2
        lab.text(x0 - 3.4, yc, LABEL[t], ha="right", va="center", fontsize=6,
                 fontweight="bold" if t == "ltr_peel" else "normal")
        lab.add_patch(Rectangle((x0 - 2.7, yc - 0.65), 1.3, 1.3, fc=COLOR[t], ec="none"))
    for text, yy in (extra or []):
        lab.text(x0 - 1.4, yy, text, ha="right", va="center", fontsize=6, color=MUTED)
    return lab


def header_ax(fig, x, y, w, h):
    ax = mm_ax(fig, x, y, w, h)
    ax.set_xlim(0, w)
    ax.set_ylim(h, 0)
    ax.axis("off")
    return ax


# ================================================================== panels A-B
def panel_A(fig, x0, y0, w, D, win):
    c, s, e, _ = win
    bpmm = (e - s) / w
    ROW, MH, GAP = 2.7, 1.5, 1.7    # call row, mask strip, gap between tools (mm)
    tracks = []
    ax = track_ax(fig, x0, y0, w, 4.5, (s, e))       # bp axis on top
    ax.plot([s, e], [4.2, 4.2], color=INK, lw=0.5, solid_capstyle="butt")
    for tk in range(s, e + 1, 25_000):
        major = tk % 50_000 == 0
        ax.plot([tk, tk], [4.2, 3.4 if major else 3.8], color=INK, lw=0.5)
        if major:
            ax.text(tk, 3.0, f"{tk / 1e6:.2f}", ha="center", va="bottom", fontsize=6)
    lab_top = track_ax(fig, 0, y0, x0 - 0.8, 4.5, (0, x0 - 0.8))
    lab_top.text(x0 - 3.4, 2.4, f"rice {WIN_CHROM_NAME.get(c, c)} (Mb)", ha="right", va="center", fontsize=6)
    y = y0 + 6.0
    gh = 2.2                                          # genes: arrows in their orientation
    ax = track_ax(fig, x0, y, w, gh, (s, e))
    for _, g in D["genes"].iterrows():
        a, b = max(g.start, s), min(g.end, e)
        head = min(1.0 * bpmm, (b - a) * 0.6)
        if g.strand == "+":
            pts = [(a, 0.35), (b - head, 0.35), (b, gh / 2), (b - head, gh - 0.35), (a, gh - 0.35)]
        else:
            pts = [(b, 0.35), (a + head, 0.35), (a, gh / 2), (a + head, gh - 0.35), (b, gh - 0.35)]
        ax.add_patch(Polygon(pts, closed=True, fc=GENE, ec="none"))
    tracks.append(("genes", y, gh))
    y += gh + 0.9
    cur = D["cur"]
    for lab, sel, col, hh in (("curated LTR-RT", cur.ltr, CUR_LTR, 2.2), ("other curated TEs", ~cur.ltr, CUR_OTH,
                                                                           1.7)):
        ax = track_ax(fig, x0, y, w, hh, (s, e))
        for _, h in cur[sel].sort_values("start").iterrows():
            a, b = max(h.start, s), min(h.end, e)
            ax.add_patch(Rectangle((a, 0.1), max(b - a - 0.15 * bpmm, 0.12 * bpmm), hh - 0.2, fc=col, ec="none"))
        tracks.append((lab, y, hh))
        y += hh + 0.6
    y += 1.4
    for t in TOOLS:
        k = D["calls"][t]
        rows = poafig.stack_rows(k.start.to_numpy(), k.end.to_numpy(), pad=int(0.6 * bpmm)) if len(k) else []
        nrow = (max(rows) + 1) if len(k) else 1
        hcalls = nrow * ROW
        ax = track_ax(fig, x0, y, w, hcalls + 0.4 + MH, (s, e))
        col, lt = COLOR[t], tint(COLOR[t], 0.75 if t == "hite" else 0.55)
        for (_, r), rr in zip(k.iterrows(), rows):
            yy, hh = rr * ROW + 0.25, ROW - 0.5
            a, b = max(r.start, s), min(r.end, e)
            ltrs = []
            if str(r.lltr_end) not in (".", "nan") and str(r.rltr_start) not in (".", "nan"):
                ltrs = [(max(r.start, s), min(int(r.lltr_end), e)), (max(int(r.rltr_start), s), min(r.end, e))]
            if int(r.supported) == 1:
                ax.add_patch(Rectangle((a, yy), b - a, hh, fc=lt, ec="none"))
                for la, lb in ltrs:
                    if lb > la:
                        ax.add_patch(Rectangle((la, yy), lb - la, hh, fc=col, ec="none"))
            else:
                ax.add_patch(Rectangle((a, yy), b - a, hh, fc="white", ec=col, lw=0.6))
                for la, lb in ltrs:
                    if lb > la:
                        ax.add_patch(Rectangle((la, yy), lb - la, hh, fc="white", ec=col, lw=0.6, hatch="||||||",
                                               hatchcolor=col))
        my = hcalls + 0.4
        m = D["masks"][t]
        for a, b in m["tp"]:
            ax.add_patch(Rectangle((a, my), b - a, MH, fc=MASK_TP, ec="none"))
        for a, b in m["fp"]:
            ax.add_patch(Rectangle((a, my), b - a, MH, fc=col, ec="none"))
        for a, b in m["fn"]:
            ax.add_patch(Rectangle((a, my + 0.12), b - a, MH - 0.24, fc="white", ec=CUR_LTR, lw=0.35))
        tracks.append((t, y, hcalls + 0.4 + MH))
        num("A", "rice", t, "window: masked outside curated LTR-RT (kb)", round(span(m["fp"]) / 1e3, 1),
            "computed: mask/rice/<tool>.intact/rice.fa.out minus curated LTR-RT hits, in the window")
        num("A", "rice", t, "window: curated LTR-RT masked (kb)", round(span(m["tp"]) / 1e3, 1), "computed (as above)")
        num("A", "rice", t, "window: calls drawn (n)", len(k), "calls/rice/<tool>.tsv overlapping the window")
        num("A", "rice", t, "window: calls not library-supported (n)", int((k.supported == 0).sum()),
            "eval/adj/rice/<tool>.calls.tsv supported")
        y += hcalls + 0.4 + MH + GAP
    lab = track_ax(fig, 0, y0, x0 - 0.8, y - y0, (0, x0 - 0.8))
    for name, yt, hh in tracks:
        yc = yt - y0 + hh / 2
        if name in TOOLS:
            lab.text(x0 - 3.4, yc, LABEL[name], ha="right", va="center", fontsize=6,
                     fontweight="bold" if name == "ltr_peel" else "normal")
            lab.add_patch(Rectangle((x0 - 2.7, yc - 0.65), 1.3, 1.3, fc=COLOR[name], ec="none"))
        else:
            lab.text(x0 - 1.4, yc, name, ha="right", va="center", fontsize=6)
    num("A", "rice", ".", "window", f"{WIN_CHROM_NAME.get(c, c)}:{s}-{e} (0-based)", "tmp/rice_window_candidates.tsv")
    num("A", "rice", ".", "window: curated LTR-RT (kb)", round(D["cur_ltr_bp"] / 1e3, 1), "mask/rice/curated/rice.fa.out")
    num("A", "rice", ".", "window: non-transposon genes (n)", len(D["genes"]), "data/rice/genes/all.gff3")
    return y - GAP


def colorbar_key(fig, x, y, w, label):
    """'label' above '0 [bar] 1', centered in the box x..x+w (mm) with its top at y; as fig_benchmark's score key"""
    ax = poafig.key_axes(fig, x, y, w, 6.4)
    xm = w / 2
    ax.text(xm, 1.5, label, ha="center", va="center", fontsize=6)
    barw = 11.0
    x0 = xm - barw / 2
    g = np.linspace(0, 1, 128)
    for k in range(len(g) - 1):
        ax.add_patch(Rectangle((x0 + barw * g[k], 3.6), barw / (len(g) - 1) + 0.02, 1.4, fc=CMAP(g[k]), ec="none"))
    ax.add_patch(Rectangle((x0, 3.6), barw, 1.4, fc="none", ec=INK, lw=0.3))
    ax.text(x0 - 0.7, 4.3, "0", ha="right", va="center", fontsize=6)
    ax.text(x0 + barw + 0.7, 4.3, "1", ha="left", va="center", fontsize=6)


def panel_Mb(fig, x0, y0, cols, colw, RH, M, reasons):
    H = len(TOOLS) * RH
    for j, ds in enumerate(DS):
        m = M[M.ds == ds].set_index("tool")
        cur = int(m.TP.iloc[0] + m.FN.iloc[0])
        num("B", ds, ".", "curated LTR-RT (Mb)", round(cur / 1e6, 1), "eval/lib/<ds>.ci.tsv TP+FN")
        fpmax = float(m.FP.max())
        sc = (colw - 4.4) / (cur + fpmax)               # mm per bp, the same both sides of the line
        xc = cur * sc
        ax = mm_ax(fig, cols[j], y0, colw, H)
        ax.set_xlim(0, colw)
        ax.set_ylim(H, 0)
        ax.axis("off")
        for i, t in enumerate(TOOLS):
            yt, hh = i * RH + 0.3, RH - 0.6
            if t not in m.index:
                missing_box(ax, 0, yt, colw, hh, reasons[(ds, t)][0])
                continue
            r = m.loc[t]
            col = COLOR[t]
            ax.add_patch(Rectangle((xc - r.TP * sc, yt), r.TP * sc, hh, fc=tint(col, 0.5), ec="none"))
            ax.add_patch(Rectangle((0, yt), cur * sc, hh, fc="none", ec="#7A7A7A", lw=0.3))
            ax.add_patch(Rectangle((xc, yt), r.FP * sc, hh, fc=col, ec="none"))
            ax.text(xc + r.FP * sc + 0.5, yt + hh / 2, f"{r.FP / 1e6:.0f}" if r.FP >= 1e6 else f"{r.FP / 1e6:.1f}",
                    ha="left", va="center", fontsize=6, fontweight="bold" if t == "ltr_peel" else "normal")
            num("B", ds, t, "masked outside curated LTR-RT (Mb)", round(r.FP / 1e6, 2), "report/tables/masking.tsv FP_mb")
            num("B", ds, t, "curated LTR-RT masked (Mb)", round(r.TP / 1e6, 2), "eval/lib/<ds>.ci.tsv TP")
        for i, t in enumerate(TOOLS):
            if t in m.index:
                ax.plot([xc, xc], [i * RH, (i + 1) * RH], color=INK, lw=0.45, solid_capstyle="butt")
        hdr = mm_ax(fig, cols[j], y0 - 7.0, colw, 6.8)
        hdr.set_xlim(0, colw)
        hdr.set_ylim(6.8, 0)
        hdr.axis("off")
        hdr.text(colw / 2, 1.9, DSNAME[ds], ha="center", va="center", fontsize=6.5)
        hdr.text(xc - 0.5, 5.5, "found", ha="right", va="center", fontsize=5.5, color=MUTED)
        hdr.text(xc + 0.5, 5.5, "outside", ha="left", va="center", fontsize=5.5, color=MUTED)
        hdr.plot([xc, xc], [6.2, 6.8], color=INK, lw=0.45)
    label_column(fig, x0, y0, H, TOOLS, RH)
    lab = mm_ax(fig, 0, y0 - 3.2, x0 - 0.8, 3.0)
    lab.set_xlim(0, x0 - 0.8)
    lab.set_ylim(3.0, 0)
    lab.axis("off")
    lab.text(x0 - 1.4, 1.5, "Mb", ha="right", va="center", fontsize=6, color=MUTED)
    return y0 + H


# ================================================================== panel C: lib-test metric scorecard
def fmt_score(v):
    """printed cell value: two decimals, as Fig. 1's scorecards (five metrics per genome leave 5.6-mm sub-columns, too
    narrow for three decimals at 6 pt)"""
    return f"{v:.2f}"


def fmt3(v):
    """three decimals, for legend numbers that two decimals would round together (recall 0.990-0.997)"""
    return f"{v:.3f}"


def panel_scores(fig, x0, y0, cols, colw, RH, L, reasons):
    """per genome column: five metric sub-columns (precision, recall, F1, specificity, accuracy), one row per tool; a
    library without a result is one cell across the five metrics ('tool failed' / 'not evaluated'). Returns bottom y."""
    H = len(TOOLS) * RH
    sw = colw / len(METRICS)
    pad = 0.13
    for j, ds in enumerate(DS):
        ax = mm_ax(fig, cols[j], y0, colw, H)
        ax.set_xlim(0, colw)
        ax.set_ylim(H, 0)
        ax.axis("off")

        def cell(k, yy, v, bold=False):
            ax.add_patch(Rectangle((k * sw + pad, yy + pad), sw - 2 * pad, RH - 2 * pad, fc=CMAP(v), ec="none"))
            ax.text(k * sw + sw / 2, yy + RH / 2, fmt_score(v), ha="center", va="center", fontsize=6,
                    color=text_color(v), fontweight="bold" if bold else "normal")
        l = L[L.ds == ds].set_index("tool")
        # bold = the column's best printed value among the tools (ties all bold; higher is better for all five), only
        # when >= 2 tools have a value (main-figure rule)
        best = {c: max(float(fmt_score(float(v))) for v in l[c]) for c, _ in METRICS}
        for i, t in enumerate(TOOLS):
            if t not in l.index:
                missing_box(ax, pad, i * RH + pad, colw - 2 * pad, RH - 2 * pad, reasons[(ds, t)][0], fs=6)
                continue
            for k, (c, _) in enumerate(METRICS):
                v = float(l.loc[t, c])
                bold = len(l) >= 2 and float(fmt_score(v)) == best[c]
                cell(k, i * RH, v, bold)
                src = f"eval/lib/{ds}/{t}.intact.tsv {c}"
                if c == "F1":
                    src += " (= report/tables/masking.tsv F1 to 4 decimals, tmp/f1_crosscheck.tsv)"
                num("C", ds, t, f"{MNAME[c]} (printed {fmt_score(v)}{', bold' if bold else ''})", round(v, 4), src)
        hdr = header_ax(fig, cols[j], y0 - 7.0, colw, 6.8)
        hdr.text(colw / 2, 1.9, DSNAME[ds], ha="center", va="center", fontsize=6.5)
        for k, (_, h) in enumerate(METRICS):
            hdr.text(k * sw + sw / 2, 5.3, h, ha="center", va="center", fontsize=6, color=MUTED)
    label_column(fig, x0, y0, H, TOOLS, RH)
    return y0 + H


# ================================================================== keys
def key_A(fig, x, y, w):
    ax = poafig.key_axes(fig, x, y, w, 7.0)
    row1 = [
        (poafig.g_poly([(0, -0.7), (1.6, -0.7), (2.4, 0), (1.6, 0.7), (0, 0.7)], GENE), "gene"),
        (poafig.g_box(CUR_LTR, w=2.4, h=1.4), "curated LTR-RT"),
        (poafig.g_box(CUR_OTH, w=2.4, h=1.2), "other curated TE"),
        (glyph_element(True), "element, library-supported"),
        (glyph_element(False), "element, unsupported"),
    ]
    row2 = [
        (poafig.g_box(MASK_TP, w=2.4, h=1.4), "masked, curated LTR-RT"),
        (glyph_tools(), "masked, not curated LTR-RT"),
        (poafig.g_box("white", w=2.4, h=1.1, ec=CUR_LTR, lw=0.35), "curated LTR-RT, unmasked"),
    ]
    poafig.hkey(ax, row1, 0, 1.5, w, row_h=3.0, gap=3.2, fontsize=6)
    poafig.hkey(ax, row2, 0, 4.5, w, row_h=3.0, gap=3.2, fontsize=6)


def glyph_element(sup):
    col = COLOR["ltr_peel"]

    def draw(ax, x, y, measure=False):
        W = 5.0
        if not measure:
            if sup:
                ax.add_patch(Rectangle((x, y - 0.75), W, 1.5, fc=tint(col, 0.55), ec="none"))
                for xx in (x, x + W - 0.9):
                    ax.add_patch(Rectangle((xx, y - 0.75), 0.9, 1.5, fc=col, ec="none"))
            else:
                ax.add_patch(Rectangle((x, y - 0.75), W, 1.5, fc="white", ec=col, lw=0.6))
                for xx in (x, x + W - 0.9):
                    ax.add_patch(Rectangle((xx, y - 0.75), 0.9, 1.5, fc="white", ec=col, lw=0.6, hatch="||||||",
                                           hatchcolor=col))
        return W
    return draw


def glyph_tools():
    cols = [COLOR[t] for t in ("ltr_peel", "ltr_retriever", "ltr_finder", "hite")]

    def draw(ax, x, y, measure=False):
        W = 2.4
        if not measure:
            for k, c in enumerate(cols):
                ax.add_patch(Rectangle((x + k * W / len(cols), y - 0.7), W / len(cols), 1.4, fc=c, ec="none"))
        return W
    return draw


# ================================================================== legend
WORD = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven", 8: "eight", 9: "nine"}


def legend_text(V):
    return f"""**Fig. S2. On real genomes, the library built from LTR_peel's elements masks nearly all curated LTR retrotransposon sequence but also much sequence that curators did not call LTR retrotransposon, so its masking precision, specificity and accuracy trail the strictest callers.** Each caller's elements served as a repeat library to mask its own genome, scored base by base against the curated annotation. Fig. 1E shows precision, recall and F1 for *A. thaliana*, rice and zebra finch, and Fig. 1F the recovery of curated families. Dark gray cells mark callers that failed (no result within 36 h or 200 GB of memory, or a crash without calls). Hatched cells were not evaluated because masking with the caller's library, or building it, would have run past 36 h or was stopped there.

**A**, A typical 200-kb window of rice chromosome 11. Elements carry their long terminal repeats (darker ends) and are outlined where curators do not support them. Beneath each caller, what its library masks here, gray inside curated LTR retrotransposons and colored outside them. {V['A']}

**B**, Curated LTR retrotransposon sequence masked by each library (left, framed by all curated sequence) and sequence masked outside it (right, Mb). {V['B']}

**C**, Precision, recall (sensitivity), F1, specificity and accuracy of each library's masking. Yellower is higher, and the best value in each column is bold where at least two callers were scored. {V['C']}"""


def legend_values(M, L, D, PCH):
    lq = M[M.tool == "ltr_peel"].set_index("ds")
    fp = {t: span(D["masks"][t]["fp"]) / 1e3 for t in TOOLS}
    k = D["calls"]["ltr_peel"]
    assert (k.supported == 1).all()
    V = {}
    V["A"] = (f"All {WORD[len(k)]} LTR_peel elements are supported, yet its library masks {fp['ltr_peel']:.0f} kb outside "
              f"curated LTR retrotransposons here, HiTE's {fp['hite']:.0f} kb.")
    sens = [lq.loc[d, "sens"] for d in ("arab", "rice", "finch")]
    V["B"] = (f"LTR_peel masks {min(sens) * 100:.1f}–{max(sens) * 100:.1f}% of the curated sequence, like the other "
              f"permissive callers, but {lq.loc['arab', 'FP_mb']:.0f}, {lq.loc['rice', 'FP_mb']:.0f} and "
              f"{lq.loc['finch', 'FP_mb']:.0f} Mb besides in *A. thaliana*, rice and zebra finch.")
    l = L.set_index(["ds", "tool"])
    pm = num("legend", "maize", ".", "curated LTR-RT share of the genome (%)", round(PCH["maize"] * 100),
             "computed: lib-test (TP+FN)/(TP+TN+FP+FN), eval/lib/maize/ltrpred.intact.tsv")
    rec = [fmt3(float(l.loc[(d, "ltr_peel"), "sens"])) for d in ("arab", "rice", "finch")]   # 3 decimals: the cells
    # print 1.00, 0.99, 0.99
    V["C"] = (f"Outside maize ({pm}% curated LTR retrotransposon), most of each genome is not LTR retrotransposon, "
              f"which compresses specificity and accuracy toward 1. LTR_peel's recall is {min(rec)}–{max(rec)}, but its "
              f"specificity, {l.loc[('arab', 'ltr_peel'), 'spec']:.2f}, {l.loc[('rice', 'ltr_peel'), 'spec']:.2f} and "
              f"{l.loc[('finch', 'ltr_peel'), 'spec']:.2f} in *A. thaliana*, rice and zebra finch, is below HiTE's "
              f"{l.loc[('arab', 'hite'), 'spec']:.2f}, {l.loc[('rice', 'hite'), 'spec']:.2f} and "
              f"{l.loc[('finch', 'hite'), 'spec']:.2f}.")
    for ds in ("arab", "rice", "finch"):                 # check the claims in that sentence
        assert l.loc[(ds, "ltr_peel"), "spec"] < l.loc[(ds, "hite"), "spec"], ds
        num("legend", ds, "ltr_peel", "recall (as printed in the legend, 3 decimals)",
            fmt3(float(l.loc[(ds, "ltr_peel"), "sens"])), f"eval/lib/{ds}/ltr_peel.intact.tsv sens")
    for ds, t in (("arab", "ltr_peel"), ("rice", "ltr_peel"), ("finch", "ltr_peel"), ("arab", "hite"),
                  ("rice", "hite"), ("finch", "hite")):
        num("legend", ds, t, "specificity (as printed in the legend)", fmt_score(float(l.loc[(ds, t), "spec"])),
            f"eval/lib/{ds}/{t}.intact.tsv spec")
    return V


# ================================================================== source data
def write_table(df, name):
    """source_data/<stem>/<name>.tsv, gzipped (mtime 0) when over 1 MB"""
    raw = df.to_csv(sep="\t", index=False)
    if len(raw.encode()) > 1_000_000:
        df.to_csv(SD / f"{name}.tsv.gz", sep="\t", index=False, compression={"method": "gzip", "mtime": 0})
    else:
        (SD / f"{name}.tsv").write_text(raw)


def read_table(name, **kw):
    p = SD / f"{name}.tsv"
    return pd.read_csv(p if p.exists() else SD / f"{name}.tsv.gz", sep="\t", **kw)


def build_source_data(reprep, verbose):
    """full mode: load the full dataset, run the raw-data checks, write the source-data tables"""
    global B, T
    B = P.init_bench()
    T = B / "report/tables"
    M = load_tables()
    L, nchk = load_libtest(M, verbose)
    reasons = missing_reasons()
    if reprep:
        P.windows(force=True)
    win, _ = window_choice(verbose)
    D = window_data(*win[:3])
    log(f"window data: {D['cur_ltr_bp']} bp curated LTR-RT; call checks ok", verbose)
    SD.mkdir(parents=True, exist_ok=True)
    c, s, e, rank = win
    write_table(pd.DataFrame([("window_chrom", c), ("window_start", s), ("window_end", e), ("window_rank", rank),
                              ("window_curated_ltr_bp", D["cur_ltr_bp"]), ("f1_crosscheck_cells", nchk)],
                             columns=["name", "value"]), "scalars")
    write_table(D["genes"][["start", "end", "strand"]], "A_genes")
    write_table(D["cur"][["start", "end"]].assign(ltr_rt=D["cur"].ltr.astype(int)), "A_curated_hits")
    write_table(pd.concat([D["calls"][t][["start", "end", "lltr_end", "rltr_start", "supported"]].assign(tool=t)
                           for t in TOOLS])[["tool", "start", "end", "lltr_end", "rltr_start", "supported"]],
                "A_calls")
    write_table(pd.DataFrame([(t, kind, a, b) for t in TOOLS for kind in ("tp", "fp", "fn")
                              for a, b in D["masks"][t][kind]], columns=["tool", "kind", "start", "end"]), "A_masks")
    write_table(M[["ds", "tool", "TP", "FN", "FP", "sens", "FP_mb"]], "B_masking")
    write_table(L[["ds", "tool", "prec", "sens", "F1", "spec", "accu", "TP", "TN", "FP", "FN"]], "C_libtest")
    write_table(pd.DataFrame([(ds, t, lab, src) for (ds, t), (lab, src) in reasons.items()],
                             columns=["ds", "tool", "label", "source"]), "BC_missing")


def read_source_data():
    """the tables the figure is drawn from: M (masking, B), L (lib-test, C), reasons (missing cells), win, D (A),
    nchk"""
    sc = dict(read_table("scalars", dtype=str).values)
    win = (sc["window_chrom"], int(sc["window_start"]), int(sc["window_end"]), int(sc["window_rank"]))
    D = {"genes": read_table("A_genes"), "cur_ltr_bp": int(sc["window_curated_ltr_bp"])}
    cur = read_table("A_curated_hits")
    D["cur"] = cur.assign(ltr=cur.ltr_rt.astype(bool))
    calls = read_table("A_calls", dtype={"lltr_end": str, "rltr_start": str})
    D["calls"] = {t: calls[calls.tool == t].reset_index(drop=True) for t in TOOLS}
    mk = read_table("A_masks")
    D["masks"] = {t: {kind: mk[(mk.tool == t) & (mk.kind == kind)][["start", "end"]].to_numpy(np.int64).reshape(-1, 2)
                      for kind in ("tp", "fp", "fn")} for t in TOOLS}
    M = read_table("B_masking")
    L = read_table("C_libtest")
    L["total"] = L.TP + L.TN + L.FP + L.FN
    miss = read_table("BC_missing")
    reasons = {(r.ds, r.tool): (r.label, r.source) for r in miss.itertuples()}
    return M, L, reasons, win, D, int(sc["f1_crosscheck_cells"])


# ================================================================== main
def main(argv=None):
    ap = argparse.ArgumentParser(description="Build Fig. S2 (figS2_masking): figure, legend page and numbers.")
    ap.add_argument("--reprep", action="store_true", help="recompute the cached window table")
    ap.add_argument("--from-source-data", action="store_true",
                    help="redraw from source_data/figS2_masking/ only (no BENCH_ROOT; raw-data checks skipped)")
    ap.add_argument("-v", "--verbose", action="store_true", help="progress and checks")
    a = ap.parse_args(argv)
    v = a.verbose
    if not a.from_source_data:
        build_source_data(a.reprep, v)
    M, L, reasons, win, D, nchk = read_source_data()
    for ds in DS:                                          # every missing result must have a reason
        for t in TOOLS:
            if t not in set(M[M.ds == ds].tool) and (ds, t) not in reasons:
                sys.exit(f"error: no reason for the missing masking result {ds} {t}")
    PCH = {}                                               # curated LTR-RT share of the scored genome (lib-test counts)
    for ds in DS:
        l = L[L.ds == ds]
        pl = (l.TP + l.FN) / l.total
        if pl.max() - pl.min() > 1e-4:
            sys.exit(f"error: {ds}: curated share differs between libraries {pl.tolist()}")
        PCH[ds] = float(pl.median())
    num("C", ".", ".", "F1 cross-check: cells where eval/lib F1, precision, sensitivity equal masking.tsv",
        f"{nchk} of {nchk}", "tmp/f1_crosscheck.tsv")

    W, Hfig = 180.0, 152.5
    x0, xr = 24.0, 176.0
    fig = plt.figure(figsize=(W * MM, Hfig * MM))
    # ---- A
    key_A(fig, 2.0, 1.0, xr - 2.0)
    yA0 = 9.5
    poafig.panel_letter(fig, 1.0, yA0, "A")
    yA = panel_A(fig, x0, yA0, xr - x0, D, win)
    # ---- shared genome columns for B and C
    gap = 3.2
    colw = (xr - x0 - 4 * gap) / len(DS)
    cols = [x0 + j * (colw + gap) for j in range(len(DS))]
    # ---- B
    yB = yA + 12.0
    poafig.panel_letter(fig, 1.0, yB - 7.4, "B")
    yBend = panel_Mb(fig, x0, yB, cols, colw, 2.5, M, reasons)
    # ---- C (scorecard) with its color key in the label column
    yC = yBend + 12.0
    poafig.panel_letter(fig, 1.0, yC - 7.4, "C")
    colorbar_key(fig, 6.0, yC - 7.4, x0 - 7.0, "score")
    yCend = panel_scores(fig, x0, yC, cols, colw, 2.6, L, reasons)
    log(f"bottom at {yCend:.1f} mm of {Hfig}", v)

    for (ds, t), (lab, src) in sorted(reasons.items()):
        num("B-C", ds, t, "missing result (label)", lab, src)

    # ---- legend
    V = legend_values(M, L, D, PCH)
    for ds in ("arab", "rice", "finch"):
        num("legend", ds, "ltr_peel", "sensitivity (curated sequence masked)",
            float(M[(M.ds == ds) & (M.tool == "ltr_peel")].sens.iloc[0]), "report/tables/masking.tsv sens")
        num("legend", ds, "ltr_peel", "masked outside curated LTR-RT (Mb)",
            float(M[(M.ds == ds) & (M.tool == "ltr_peel")].FP_mb.iloc[0]), "report/tables/masking.tsv FP_mb")
    short = legend_text(V)
    src = B if B is not None else SD
    out = poafig.write_figure(STEM, fig, short, None, HERE, str(src), str(HERE / f"{STEM}.py"))
    pd.DataFrame(NUM).to_csv(HERE / f"{STEM}_numbers.tsv", sep="\t", index=False)
    log(f"wrote {out['pdf']} ({out['fig_mm']}, legend {out['short_mm']} mm, {out['short_words']} words, per "
        f"paragraph {out['short_par_words']})", True)


if __name__ == "__main__":
    main()
