#!/usr/bin/env python3
"""Fig. S5 figS5_resources — run time, CPU time and peak memory of nine LTR retrotransposon callers.

One row per tool (shared order), one column per measure (A wall time, B CPU time, C peak memory); each cell draws that
tool's 13 runs against genome size (log-log): real genomes as filled circles joined in size order, each simulated
genome (3 seeds x 2 regimes) as an open circle, failed runs (36 h cap, 200 GB memory guard, 0 calls and an error
exit) as one dark gray 'tool failed' marker at the value reached (reasons in the legend only); the other eight tools'
real-genome lines in gray behind. Caps drawn as dashed lines.

Inputs (read-only): report/tables/resources.tsv (status, threads, load_at_start; drawn values checked against it),
runs/<ds>/<tool>/resources.json (the table's source; full precision for the log axis), data/<ds>/genome.fa.fai.
Outputs: figS5_resources.{pdf,png} (poafig.write_figure: figure + one legend page), figS5_resources_numbers.tsv.
check_claims() stops the run, before anything is written, if a number worded in the legend no longer holds.
Usage: figS5_resources.py [-v]
"""
import argparse
import glob
import json
import os
import sys
from pathlib import Path

sys.dont_write_bytecode = True                       # no __pycache__ next to the read-only imported modules

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # poafig.py sits beside this script
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
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

STEM = "figS5_resources"
B = Path(os.environ["BENCH_ROOT"])
HERE = Path(__file__).resolve().parent
TAB = B / "report/tables/resources.tsv"
MM = poafig.MM
INK = poafig.INK
MUTED = "#6B6B6B"
CONTEXT = "#C8C8C8"       # other tools' lines
CAPC = "#8C8C8C"          # cap / guard lines
CONN = "#A0A0A0"          # genome-name connectors
GUIDE = "#EFEFEF"         # faint vertical guides at each genome size
FAILC = "#5A5A5A"         # tool failed (missing-result style, tool side)
FAIL_STYLE = dict(marker="X", ms=4.2, mfc=FAILC, mec="white", mew=0.3)

TOOLS = ["ltrquest", "ltr_retriever", "ltrharvest", "ltr_finder", "look4ltrs", "hirepeat", "hite", "ltrpred",
         "dante_ltr"]
LABEL = {"ltrquest": "LTRquest", "ltr_retriever": "LTR_retriever", "ltrharvest": "LTRharvest",
         "ltr_finder": "LTR_Finder", "look4ltrs": "Look4LTRs", "hirepeat": "HiRepEat", "hite": "HiTE",
         "ltrpred": "LTRpred", "dante_ltr": "DANTE_LTR"}
COLOR = {"ltrquest": "#D55E00", "ltr_retriever": "#0072B2", "ltrharvest": "#56B4E9", "ltr_finder": "#009E73",
         "look4ltrs": "#CC79A7", "hirepeat": "#E69F00", "hite": "#000000", "ltrpred": "#A6761D",
         "dante_ltr": "#6A3D9A"}
REAL = ["arab", "asue", "rice", "finch", "pann", "maize", "gorilla"]
SIMS = {"S1": ["S1_1", "S1_2", "S1_3"], "S2": ["S2_1", "S2_2", "S2_3"]}
SIM_DS = SIMS["S1"] + SIMS["S2"]
# genome sizes (bp); asserted against the .fai sums
GENOME_BP = {"rice": 374471240, "arab": 119668634, "maize": 2182075994, "finch": 1141165377, "gorilla": 3545850631,
           "asue": 262563758, "pann": 1805882605}
NAME = {"arab": r"$\it{A.\ thaliana}$", "S1": "simulated S1", "asue": r"$\it{A.\ suecica}$", "S2": "simulated S2",
        "rice": "rice", "finch": "zebra finch", "pann": r"$\it{P.\ annua}$", "maize": "maize", "gorilla": "gorilla"}
CAP_H = 36.0              # wall-time cap (h)
GUARD_GB = 200.0          # PSS memory guard (GB)
THREADS = 16
CPU_CEIL = CAP_H * THREADS

# measures: key, column header, y-limits, ticks (value, label), cap line (value, style, label)
MEASURES = [
    ("wall_h", "wall time", (1.1e-3, 90.0),
     [(1 / 60, "1 min"), (1 / 6, "10 min"), (1.0, "1 h"), (10.0, "10 h")],
     (CAP_H, (0, (3, 2)), "36 h cap")),
    ("cpu_h", "CPU time", (8e-3, 1600.0),
     [(1 / 60, "1 min"), (1.0, "1 h"), (10.0, "10 h"), (100.0, "100 h")],
     (CPU_CEIL, (0, (1, 1.5)), "16 threads × 36 h")),
    ("mem_gb", "peak memory", (0.45, 420.0),
     [(1.0, "1 GB"), (10.0, "10 GB"), (100.0, "100 GB")],
     (GUARD_GB, (0, (3, 2)), "200 GB guard")),
]

# geometry (mm)
FIG_W = 180.0
NAME_W = 19.5             # tool-name column (name + color square)
TICK_W = 9.0              # y tick labels of each column
GAP_W = 1.2
RIGHT_PAD = 2.0
CELL_W = (FIG_W - NAME_W - 3 * TICK_W - 2 * GAP_W - RIGHT_PAD) / 3
LABEL_MARGIN = 1.5        # right part of each cell past the data (room for the gorilla markers)
X_LO, X_HI = 0.095, 3.95  # Gb, data part of the cell
NAME_GAP = 2.5            # min spacing of the rotated genome names (mm)
Y_KEY = 3.0
Y_HEAD = 11.0
Y_NAMES = 30.0            # bottom of the rotated genome names
Y0 = 34.5                 # top of the first row of cells
ROW_H = 18.0
ROW_GAP = 2.6
FIG_H = Y0 + 9 * ROW_H + 8 * ROW_GAP + 10.0


def log(msg, verbose):
    if verbose:
        print(msg, file=sys.stderr, flush=True)


# ---------------------------------------------------------------------------------------------------------- data
def genome_size_gb(ds):
    bp = sum(int(line.split("\t")[1]) for line in open(B / f"data/{ds}/genome.fa.fai"))
    if ds in GENOME_BP and bp != GENOME_BP[ds]:
        sys.exit(f"error: {ds} genome size {bp} != expected size {GENOME_BP[ds]}")
    return bp / 1e9


def load():
    """117 runs: values from runs/<ds>/<tool>/resources.json with summarize.py's formulas, checked against the
    rounded report table; status, threads and load from the table."""
    tab = pd.read_csv(TAB, sep="\t", dtype=str)
    tab = tab[tab.tool.isin(TOOLS) & tab.ds.isin(REAL + SIM_DS)]
    if len(tab) != len(TOOLS) * len(REAL + SIM_DS):
        sys.exit(f"error: resources.tsv has {len(tab)} of {len(TOOLS) * len(REAL + SIM_DS)} expected rows")
    size = {d: genome_size_gb(d) for d in REAL + SIM_DS}
    rows = []
    for r in tab.itertuples():
        p = B / f"runs/{r.ds}/{r.tool}/resources.json"
        j = json.load(open(p))
        wall = j["wall_s"] / 3600
        cpu = max(j["cpu_s"], j.get("cpu_proc_s") or 0) / 3600
        mem = j["peak_pss_gb"]
        for v, col, fmt in ((wall, "wall_h", "{:.2f}"), (cpu, "cpu_h", "{:.2f}"), (mem, "peak_pss_gb", "{:.1f}")):
            if fmt.format(v) != getattr(r, col):
                sys.exit(f"error: {r.ds}/{r.tool} {col}: json {v} does not round to the table's {getattr(r, col)}")
        if j["status"] != r.status:
            sys.exit(f"error: {r.ds}/{r.tool} status json {j['status']} != table {r.status}")
        rows.append(dict(ds=r.ds, tool=r.tool, status=r.status, wall_h=wall, cpu_h=cpu, mem_gb=mem,
                         rss_gb=float(r.peak_rss_gb), threads=int(r.threads), load=int(r.load_at_start),
                         size_gb=size[r.ds], real=r.ds in REAL,
                         group="real" if r.ds in REAL else r.ds.split("_")[0]))
    df = pd.DataFrame(rows)
    bad = set(df.status) - {"ok", "timeout", "mem_guard", "failed"}
    if bad:
        sys.exit(f"error: unexpected status values {bad}")
    return df, size


def facts(df):
    """every number worded in the legend, computed from the table (also written to the numbers file)"""
    ok = df[df.status == "ok"]
    gor = df[df.ds == "gorilla"]
    gor_ok = gor[gor.status == "ok"]
    lq = df[df.tool == "ltrquest"].set_index("ds")
    fold = {}
    for d in REAL:
        w = ok[ok.ds == d].wall_h
        fold[d] = w.max() / w.min()
    lp = ok[ok.tool == "ltrpred"]
    s2_over_rice = {}
    for t in TOOLS:
        s2 = df[(df.tool == t) & df.ds.isin(SIMS["S2"])].wall_h
        rice = df[(df.tool == t) & (df.ds == "rice")].wall_h.iloc[0]
        s2_over_rice[t] = s2.min() / rice
    return dict(
        not_ok={(r.ds, r.tool, r.status) for r in df[df.status != "ok"].itertuples()},
        s2_over_rice=s2_over_rice,
        s2_slower_all_seeds=sorted(t for t, v in s2_over_rice.items() if v > 1),
        n_runs=len(df), n_stopped=int(df.status.isin(["timeout", "mem_guard"]).sum()),
        n_timeout=int((df.status == "timeout").sum()), n_guard=int((df.status == "mem_guard").sum()),
        n_failed=int((df.status == "failed").sum()),
        gor_stopped=int(gor.status.isin(["timeout", "mem_guard"]).sum()),
        gor_ok_wall_min=gor_ok.wall_h.min(), gor_ok_wall_max=gor_ok.wall_h.max(),
        fold_min=min(fold.values()), fold_max=max(fold.values()), fold=fold,
        lq_status_all_ok=bool((df[df.tool == "ltrquest"].status == "ok").all()),
        lq_gor_wall=lq.loc["gorilla", "wall_h"], lq_gor_cpu=lq.loc["gorilla", "cpu_h"],
        max_cpu_run=df.loc[df.cpu_h.idxmax(), ["ds", "tool"]].tolist(), max_cpu=df.cpu_h.max(),
        lq_maize_mem=lq.loc["maize", "mem_gb"], max_ok_mem_run=ok.loc[ok.mem_gb.idxmax(), ["ds", "tool"]].tolist(),
        max_ok_mem=ok.mem_gb.max(),
        lf_mem_max=df[df.tool == "ltr_finder"].mem_gb.max(),
        l4_gor_mem=df[(df.tool == "look4ltrs") & (df.ds == "gorilla")].mem_gb.iloc[0],
        l4_gor_wall=df[(df.tool == "look4ltrs") & (df.ds == "gorilla")].wall_h.iloc[0],
        ltrpred_cores_max=(lp.cpu_h / lp.wall_h).max(), ltrpred_cores_min=(lp.cpu_h / lp.wall_h).min(),
        harvest_threads=sorted(df[df.tool == "ltrharvest"].threads.unique().tolist()),
        other_threads=sorted(df[df.tool != "ltrharvest"].threads.unique().tolist()),
        load_min=int(df.load.min()), load_max=int(df.load.max()),
    )


def check_claims(f):
    """stop before writing anything if a legend statement no longer holds"""
    claims = [
        (f["n_runs"] == 117, "117 runs"),
        (f["n_timeout"] == 5 and f["n_guard"] == 1 and f["n_failed"] == 2, "5 timeouts, 1 guard, 2 failed"),
        (f["not_ok"] == {("finch", "ltr_retriever", "timeout"), ("maize", "ltr_retriever", "timeout"),
                         ("gorilla", "ltr_retriever", "timeout"), ("gorilla", "ltrharvest", "timeout"),
                         ("gorilla", "ltrpred", "timeout"), ("gorilla", "look4ltrs", "mem_guard"),
                         ("finch", "dante_ltr", "failed"), ("gorilla", "dante_ltr", "failed")},
         "the failed runs named in the legend"),
        (f["gor_stopped"] == 4, "gorilla: 4 of 9 tools stopped at a limit"),
        (f"{f['gor_ok_wall_min']:.1f}" == "0.5" and f"{f['gor_ok_wall_max']:.1f}" == "33.4",
         "gorilla finished runs 0.5-33.4 h"),
        (f["lq_status_all_ok"], "LTRquest finished every genome"),
        (f"{f['lq_gor_wall']:.1f}" == "33.4" and f"{f['lq_gor_cpu']:.0f}" == "507", "LTRquest gorilla 33.4 h, 507 CPU h"),
        (f["max_cpu_run"] == ["gorilla", "ltrquest"], "LTRquest gorilla = most CPU of any run"),
        (f["max_ok_mem_run"] == ["maize", "ltrquest"] and f"{f['max_ok_mem']:.1f}" == "92.6",
         "LTRquest maize 92.6 GB = highest finished-run memory"),
        (f"{f['lf_mem_max']:.1f}" == "1.9", "LTR_Finder at most 1.9 GB"),
        (f["l4_gor_mem"] > GUARD_GB and f"{f['l4_gor_wall']:.1f}" == "13.2", "Look4LTRs gorilla over 200 GB at 13.2 h"),
        (f["ltrpred_cores_max"] < 1.05, "LTRpred about one core"),
        (f["harvest_threads"] == [1] and f["other_threads"] == [THREADS], "16 threads; LTRharvest 1"),
        (f["load_min"] == 14 and f["load_max"] == 250, "load at start 14-250"),
        (round(f["fold_min"]) == 26 and round(f["fold_max"]) == 176, "26- to 176-fold spread per real genome"),
        (set(f["s2_slower_all_seeds"]) >= {"hite", "ltrpred", "dante_ltr"}
         and all(f["s2_over_rice"][t] > 1.5 for t in ("hite", "ltrpred", "dante_ltr")),
         "HiTE, LTRpred, DANTE_LTR slower on every S2 seed than on rice (> 1.5x)"),
    ]
    failed = [c for ok, c in claims if not ok]
    if failed:
        sys.exit("error: legend claims no longer hold: " + "; ".join(failed))


# ---------------------------------------------------------------------------------------------------------- drawing
def col_left(j):
    return NAME_W + TICK_W + j * (TICK_W + CELL_W + GAP_W)


def x_mm(j, gb):
    """genome size (Gb) -> mm from the figure's left edge, for column j"""
    data_w = CELL_W - LABEL_MARGIN
    return col_left(j) + (np.log10(gb) - np.log10(X_LO)) / (np.log10(X_HI) - np.log10(X_LO)) * data_w


def x_lim_right():
    data_w = CELL_W - LABEL_MARGIN
    return 10 ** (np.log10(X_LO) + (np.log10(X_HI) - np.log10(X_LO)) * CELL_W / data_w)


def fig_line(fig, xs, ys, **kw):
    W, H = FIG_W, FIG_H
    fig.add_artist(Line2D([x / W for x in xs], [1 - y / H for y in ys], transform=fig.transFigure, **kw))


def fig_text(fig, x, y, s, **kw):
    fig.text(x / FIG_W, 1 - y / FIG_H, s, **kw)


def draw_cell(ax, df, tool, key, meas, row, col, guides):
    _, _, ylim, ticks, cap = meas
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(X_LO, x_lim_right())
    ax.set_ylim(*ylim)
    # faint guides at each genome's size, so positions read down the column
    for gb in guides:
        ax.axvline(gb, color=GUIDE, lw=0.3, zorder=0)
    # cap / guard / ceiling, over the data part only (the reasons sit right of it)
    ax.hlines(cap[0], X_LO, X_HI, color=CAPC, lw=0.5, ls=cap[1], zorder=1)
    if row == 0:
        ax.text(X_LO * 1.06, cap[0] / 1.18, cap[2], fontsize=6, color=MUTED, va="top", ha="left", zorder=4)
    # other tools: real-genome lines through finished runs
    for t2 in TOOLS:
        if t2 == tool:
            continue
        s = df[(df.tool == t2) & df.real & (df.status == "ok")].sort_values("size_gb")
        ax.plot(s.size_gb, s[key], "-", color=CONTEXT, lw=0.35, zorder=2, solid_capstyle="round")
    c = COLOR[tool]
    lw = 1.25 if tool == "ltrquest" else 0.75
    me = df[df.tool == tool]
    real = me[me.real].sort_values("size_gb")
    ok = real[real.status == "ok"]
    ax.plot(ok.size_gb, ok[key], "-", color=c, lw=lw, zorder=5, solid_joinstyle="round")
    ax.plot(ok.size_gb, ok[key], "o", ms=2.9, mfc=c, mec="white", mew=0.35, zorder=6)
    sim = me[~me.real]
    ax.plot(sim.size_gb, sim[key], "o", ms=2.7, mfc="white", mec=c, mew=0.55, zorder=6)
    # tool failed (timeout 36 h, memory guard 200 GB, 0 calls + error exit): one style, at the value reached;
    # the reasons are in the legend only 
    for r in real[real.status != "ok"].itertuples():
        y = getattr(r, key)
        if key == "wall_h" and r.status == "timeout":
            y = CAP_H                       # drawn on the cap (measured 36.00-36.01 h)
        ax.plot([r.size_gb], [y], **FAIL_STYLE, zorder=7)
    # axes
    ax.set_yticks([v for v, _ in ticks])
    ax.set_yticklabels([lab for _, lab in ticks])
    ax.yaxis.set_minor_locator(mpl.ticker.NullLocator())
    ax.xaxis.set_minor_locator(mpl.ticker.NullLocator())
    xt = [(0.1, "100 Mb"), (0.3, "300 Mb"), (1.0, "1 Gb"), (3.0, "3 Gb")]
    ax.set_xticks([v for v, _ in xt])
    ax.set_xticklabels([lab for _, lab in xt] if row == len(TOOLS) - 1 else [])
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.spines["bottom"].set_bounds(X_LO, X_HI)
    ax.tick_params(axis="both", length=1.8, pad=1.2)


def draw_key(fig):
    items = [
        (dict(marker="o", ls="-", color=INK, lw=0.75, ms=2.9, mfc=INK, mec="white", mew=0.35), "real genome"),
        (dict(marker="o", ls="none", ms=2.7, mfc="white", mec=INK, mew=0.55), "simulated genome"),
        (dict(ls="none", **FAIL_STYLE), "tool failed"),
        (dict(ls="-", color=CONTEXT, lw=0.6), "other tools"),
        (dict(ls=(0, (3, 2)), color=CAPC, lw=0.6), "limit"),
    ]
    W, H = FIG_W, FIG_H
    ax = fig.add_axes([0, 1 - (Y_KEY + 4) / H, 1, 4 / H])
    ax.set_axis_off()
    ax.set_xlim(0, W)
    ax.set_ylim(0, 4)
    x = NAME_W + TICK_W
    renderer = fig.canvas.get_renderer()
    for style, lab in items:
        if style.get("ls") not in ("none",) and "marker" not in style:
            ax.plot([x, x + 5], [2, 2], **style)
        elif style.get("ls") == "-":
            ax.plot([x, x + 5], [2, 2], **style)
            style = {**style, "ls": "none"}
            ax.plot([x + 2.5], [2], **style)
        else:
            ax.plot([x + 2.5], [2], **style)
        t = ax.text(x + 6.2, 2, lab, fontsize=6.5, va="center", ha="left")
        bb = t.get_window_extent(renderer)
        x = x + 6.2 + bb.width / fig.dpi / MM + 5.0


def draw_genome_names(fig, j, size):
    """rotated genome names above the column at their true x, pushed apart to NAME_GAP where sizes crowd, each
    joined to its true x by a thin connector"""
    order = sorted(["arab", "asue", "rice", "finch", "pann", "maize", "gorilla", "S1", "S2"],
                   key=lambda k: size[k])
    xt = np.array([x_mm(j, size[k]) for k in order])
    xl = xt.copy()
    for _ in range(50):                      # relax: push apart, then re-center each crowded group on its targets
        for i in range(1, len(xl)):
            xl[i] = max(xl[i], xl[i - 1] + NAME_GAP)
        groups, start = [], 0
        for i in range(1, len(xl) + 1):
            if i == len(xl) or xl[i] - xl[i - 1] > NAME_GAP + 1e-6:
                groups.append((start, i))
                start = i
        new = xl.copy()
        for a, b in groups:
            shift = xt[a:b].mean() - xl[a:b].mean()
            new[a:b] = xl[a:b] + shift
        if np.allclose(new, xl):
            break
        xl = new
    for k, a, b in zip(order, xl, xt):
        fig_text(fig, a, Y_NAMES, NAME[k], fontsize=6, rotation=90, ha="center", va="bottom")
        fig_line(fig, [a, a, b, b], [Y_NAMES + 0.4, Y_NAMES + 1.4, Y0 - 1.4, Y0 - 0.2], color=CONN, lw=0.3)


def draw(df, size_all, verbose):
    size = dict(size_all)
    size["S1"] = float(np.mean([size_all[d] for d in SIMS["S1"]]))
    size["S2"] = float(np.mean([size_all[d] for d in SIMS["S2"]]))
    fig = plt.figure(figsize=(FIG_W * MM, FIG_H * MM))
    draw_key(fig)
    for j, meas in enumerate(MEASURES):
        cl = col_left(j)
        poafig.panel_letter(fig, cl - TICK_W, Y_HEAD - 1.6, "ABC"[j])
        fig_text(fig, cl - TICK_W + 4.0, Y_HEAD, meas[1], fontsize=6.5, va="center", ha="left")
        draw_genome_names(fig, j, size)
        for i, tool in enumerate(TOOLS):
            y = Y0 + i * (ROW_H + ROW_GAP)
            ax = poafig.mm_axes(fig, cl, y, CELL_W, ROW_H)
            draw_cell(ax, df, tool, meas[0], meas, i, j, [size[k] for k in REAL + ["S1", "S2"]])
            if i == len(TOOLS) - 1:
                ax.set_xlabel("genome size", fontsize=6.5, labelpad=1.5)
                ax.xaxis.set_label_coords(0.5 * (CELL_W - LABEL_MARGIN) / CELL_W, -0.25)
    for i, tool in enumerate(TOOLS):
        y = Y0 + i * (ROW_H + ROW_GAP) + ROW_H / 2
        # tool label = dark name + tool-color square 
        fig_text(fig, NAME_W - 3.6, y, LABEL[tool], fontsize=6.5, color=INK, ha="right", va="center",
                 fontweight="bold" if tool == "ltrquest" else "normal")
        sq = 2.0
        fig.add_artist(Rectangle(((NAME_W - 0.4 - sq) / FIG_W, 1 - (y + sq / 2) / FIG_H), sq / FIG_W, sq / FIG_H,
                                 transform=fig.transFigure, fc=COLOR[tool], ec="none"))
    log(f"figure {FIG_W:.0f} x {FIG_H:.1f} mm; cell {CELL_W:.1f} x {ROW_H:.1f} mm", verbose)
    return fig


# ---------------------------------------------------------------------------------------------------------- outputs
def numbers(df, f, size):
    src_run = "report/tables/resources.tsv (drawn at full precision from runs/{ds}/{tool}/resources.json)"
    rows = []
    for r in df.sort_values(["tool", "size_gb"], key=lambda s: s.map(TOOLS.index) if s.name == "tool" else s).itertuples():
        drawn = {"ok": "filled circle" if r.real else "open circle", "timeout": "tool failed (36 h cap)",
                 "mem_guard": "tool failed (200 GB guard)", "failed": "tool failed (0 calls, error exit)"}[r.status]
        src = src_run.format(ds=r.ds, tool=r.tool)
        rows += [("wall time (h)", r.ds, LABEL[r.tool], r.status, f"{r.wall_h:.4f}", "A " + drawn, src),
                 ("CPU time (h)", r.ds, LABEL[r.tool], r.status, f"{r.cpu_h:.4f}", "B " + drawn, src),
                 ("peak memory PSS (GB)", r.ds, LABEL[r.tool], r.status, f"{r.mem_gb:.3f}", "C " + drawn, src),
                 ("peak memory RSS (GB)", r.ds, LABEL[r.tool], r.status, f"{r.rss_gb:.1f}", "not drawn", src),
                 ("threads", r.ds, LABEL[r.tool], r.status, str(r.threads), "legend", "report/tables/resources.tsv"),
                 ("load at start (of 256 cores)", r.ds, LABEL[r.tool], r.status, str(r.load), "legend (range)",
                  "report/tables/resources.tsv")]
    for d in REAL + SIM_DS:
        rows.append(("genome size (Gb)", d, ".", ".", f"{size[d]:.4f}", "x axis", f"data/{d}/genome.fa.fai"))
    T = "computed from report/tables/resources.tsv"
    leg = [
        ("runs drawn", f["n_runs"], T), ("runs stopped at the 36 h cap", f["n_timeout"], T),
        ("runs stopped by the 200 GB guard", f["n_guard"], T), ("runs with 0 calls and an error exit", f["n_failed"], T),
        ("gorilla: tools stopped at a limit (of 9)", f["gor_stopped"], T),
        ("gorilla: shortest finished run (h)", f"{f['gor_ok_wall_min']:.2f}", T),
        ("gorilla: longest finished run (h)", f"{f['gor_ok_wall_max']:.2f}", T),
        ("fold spread of finished wall times, smallest over real genomes", f"{f['fold_min']:.1f}", T),
        ("fold spread of finished wall times, largest over real genomes", f"{f['fold_max']:.1f}", T),
        *[(f"fold spread of finished wall times, {d}", f"{f['fold'][d]:.1f}", T) for d in REAL],
        *[(f"{LABEL[t]}: fastest S2 seed wall / rice wall",
           f"{f['s2_over_rice'][t]:.2f}", T) for t in TOOLS],
        ("LTRquest gorilla wall (h)", f"{f['lq_gor_wall']:.2f}", T),
        ("LTRquest gorilla CPU (h), most of any run", f"{f['lq_gor_cpu']:.1f}", T),
        ("LTRquest maize peak PSS (GB), most of any finished run", f"{f['lq_maize_mem']:.1f}", T),
        ("LTR_Finder largest peak PSS (GB)", f"{f['lf_mem_max']:.1f}", T),
        ("Look4LTRs gorilla peak PSS at the stop (GB)", f"{f['l4_gor_mem']:.1f}", T),
        ("Look4LTRs gorilla wall at the stop (h)", f"{f['l4_gor_wall']:.2f}", T),
        ("LTRpred CPU/wall, finished runs, min", f"{f['ltrpred_cores_min']:.3f}", T),
        ("LTRpred CPU/wall, finished runs, max", f"{f['ltrpred_cores_max']:.3f}", T),
        ("load at job start, min (of 256 cores)", f["load_min"], T),
        ("load at job start, max (of 256 cores)", f["load_max"], T),
        ("threads per tool", THREADS, "bench/memo.md Fixed decisions; resources.tsv threads"),
        ("LTRharvest threads (gt has no thread option)", 1, "resources.tsv threads; bench/memo.md DNFs"),
        ("wall-time cap (h)", CAP_H, "bench/memo.md Fixed decisions"),
        ("memory guard, PSS (GB)", GUARD_GB, "bench/memo.md Fixed decisions"),
        ("CPU ceiling line, 16 threads x 36 h (h)", CPU_CEIL, "derived: threads x cap"),
        ("finch LTR_retriever: calls complete at (h)", 2.6, "bench/memo.md DNFs"),
        ("maize LTR_retriever: calls complete at (h)", 22.1, "bench/memo.md DNFs"),
    ]
    for item, val, src in leg:
        rows.append((item, ".", ".", ".", str(val), "legend", src))
    out = pd.DataFrame(rows, columns=["item", "ds", "tool", "status", "value", "drawn_as", "source"])
    out.to_csv(HERE / f"{STEM}_numbers.tsv", sep="\t", index=False)
    return len(out)


def legend(f):
    return (
        "**Fig. S5. Run time and memory of nine LTR retrotransposon callers.** "
        "Cost rises with genome size, but on each real genome the fastest and slowest finished callers differ "
        f"{f['fold_min']:.0f}- to {f['fold_max']:.0f}-fold. On the 3.5 Gb gorilla genome {f['gor_stopped']} of 9 "
        "callers stopped at the 36 h or 200 GB limit, and one found no elements.\n\n"
        "**A**, Wall-clock time of every run against genome size, one row per caller, the other callers in gray. "
        "Real genomes (five with curated repeat libraries, two allopolyploids) are joined in order of size; open "
        "circles, three simulated genomes per regime. HiTE, LTRpred and DANTE_LTR took longer on the LTR-rich "
        "simulated S2 genomes than on the larger rice genome. LTRquest finished every genome, taking "
        f"{f['lq_gor_wall']:.1f} h on gorilla.\n\n"
        "**B**, CPU time of the same runs, a lower bound (the larger of two measures, each of which misses some "
        f"short-lived or killed processes). LTRquest used {f['lq_gor_cpu']:.0f} CPU hours on gorilla, the most of "
        "any run. LTRharvest runs on one thread, and LTRpred used about one core although given 16.\n\n"
        "**C**, Peak memory, counting memory shared between processes once. LTR_Finder never needed more than "
        f"{f['lf_mem_max']:.1f} GB. LTRquest reached {f['lq_maize_mem']:.1f} GB on maize, the most of any finished "
        "run.\n\n"
        "Gray crosses, tool failed, drawn at the value reached when the run ended. Five runs stopped at the 36 h "
        "limit: LTR_retriever on zebra finch, maize and gorilla, and LTRharvest and LTRpred on gorilla. The three "
        "gorilla runs were all still in the single-threaded LTRharvest search, which LTR_retriever and LTRpred also "
        "run; on zebra finch and maize LTR_retriever had finished its calls (at 2.6 h and 22.1 h) and spent the rest "
        "in its default genome annotation. Look4LTRs passed the 200 GB memory guard on gorilla. DANTE_LTR found no "
        "complete element in the two animal genomes and exited with an error. For runs stopped at a limit the true "
        "cost is higher. Every caller ran with 16 threads on a shared server "
        f"(load at job start {f['load_min']} to {f['load_max']} of 256 cores), so wall times are noisy. PrinTE, "
        "which made the simulated genomes, shares an author with LTRquest, and LTRquest was tuned on PrinTE "
        "simulations and on rice."
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("-v", "--verbose", action="store_true", help="report sizes and checks")
    a = ap.parse_args()
    df, size = load()
    f = facts(df)
    check_claims(f)
    log(f"{len(df)} runs; statuses {df.status.value_counts().to_dict()}", a.verbose)
    fig = draw(df, size, a.verbose)
    n = numbers(df, f, size)
    res = poafig.write_figure(STEM, fig, legend(f), None, HERE, str(B), str(HERE / f"{STEM}.py"))
    print(f"wrote {res['pdf']}, {res['png'][0]}, {STEM}_numbers.tsv ({n} rows); figure {res['fig_mm'][0]} mm, "
          f"legend page {res['short_mm']} mm, {res['short_words']} words {res['short_par_words']}")


if __name__ == "__main__":
    main()
