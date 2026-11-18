"""Plotting helpers for the benchmark figures: matplotlib style, mm-based axes, keys and the figure + legend PDF writer."""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import re
import sys
from pathlib import Path

import numpy as np

INK = "#1A1A1A"
MUTED = "#6B6B6B"
MM = 1 / 25.4


def log(msg, verbose=True):
    if verbose:
        print(msg, file=sys.stderr, flush=True)


def setup_mpl():
    import matplotlib as mpl
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Liberation Sans", "Arial"],
        "mathtext.fontset": "custom",
        "mathtext.rm": "Liberation Sans",
        "mathtext.it": "Liberation Sans:italic",
        "mathtext.bf": "Liberation Sans:bold",
        "mathtext.cal": "Liberation Sans:italic",
        "mathtext.sf": "Liberation Sans",
        "mathtext.tt": "Liberation Mono",
        "axes.facecolor": "none",
        "font.size": 6.5,
        "axes.labelsize": 6.5,
        "axes.linewidth": 0.5,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "xtick.major.width": 0.5,
        "ytick.major.width": 0.5,
        "xtick.major.size": 2,
        "ytick.major.size": 2,
        "xtick.color": INK,
        "ytick.color": INK,
        "axes.edgecolor": INK,
        "axes.labelcolor": INK,
        "text.color": INK,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "hatch.linewidth": 0.4,
        "legend.frameon": False,
        "figure.dpi": 100,
        "savefig.dpi": 300,
    })

def stack_rows(starts, ends, pad):
    """Greedy IGV-'expanded' row assignment; elements closer than pad (bp) go to separate rows."""
    order = np.argsort(starts, kind="stable")
    row_end = []
    rows = np.zeros(len(starts), dtype=int)
    for i in order:
        for r, re_ in enumerate(row_end):
            if starts[i] >= re_ + pad:
                rows[i] = r
                row_end[r] = ends[i]
                break
        else:
            rows[i] = len(row_end)
            row_end.append(ends[i])
    return rows

def mm_axes(fig, x, y, w, h, **kw):
    """Axes placed in mm from the figure's top-left corner."""
    W, H = fig.get_size_inches() / MM
    return fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H], **kw)


def panel_letter(fig, x, y, letter):
    W, H = fig.get_size_inches() / MM
    fig.text(x / W, 1 - y / H, letter, fontsize=8, fontweight="bold", va="top", ha="left")


def hkey(ax, items, x0, y, width, row_h=2.6, gap=2.6, fontsize=6):
    """Horizontal key with measured label widths, wrapping to a new row past width (axes in mm,
    y increasing downward). items: (draw(ax, x, y) -> glyph width mm, label)."""
    r = ax.figure.canvas.get_renderer()
    inv = ax.transData.inverted()
    x = x0
    for draw, lab in items:
        t = ax.text(0, y, lab, va="center", ha="left", fontsize=fontsize)
        bb = t.get_window_extent(renderer=r)
        (a, _), (b, _) = inv.transform([(bb.x0, bb.y0), (bb.x1, bb.y0)])
        tw = abs(b - a)
        gw = draw(None, 0, 0, measure=True)
        if x > x0 and x + gw + 0.7 + tw > x0 + width:
            x = x0
            y += row_h
        draw(ax, x, y)
        t.set_position((x + gw + 0.7, y))
        x += gw + 0.7 + tw + gap
    return y


def g_box(fc, w=1.4, h=1.4, ec="none", lw=0.0, hatch=None, alpha=1.0):
    from matplotlib.patches import Rectangle

    def draw(ax, x, y, measure=False):
        if not measure:
            ax.add_patch(Rectangle((x, y - h / 2), w, h, fc=fc, ec=ec, lw=lw, hatch=hatch, alpha=alpha))
        return w
    return draw


def g_line(col, w=2.2, lw=0.9, ls="-"):
    def draw(ax, x, y, measure=False):
        if not measure:
            ax.plot([x, x + w], [y, y], color=col, lw=lw, ls=ls, solid_capstyle="butt")
        return w
    return draw


def g_poly(pts, fc, alpha=1.0):
    """pts: (dx, dy) in mm relative to (x, y)."""
    from matplotlib.patches import Polygon
    w = max(p[0] for p in pts)

    def draw(ax, x, y, measure=False):
        if not measure:
            ax.add_patch(Polygon([(x + a, y + b) for a, b in pts], fc=fc, ec="none", alpha=alpha))
        return w
    return draw

def key_axes(fig, x, y, w, h):
    ax = mm_axes(fig, x, y, w, h)
    ax.set_xlim(0, w)
    ax.set_ylim(h, 0)
    ax.axis("off")
    return ax


def fig_line(fig, pts_mm, **kw):
    from matplotlib.lines import Line2D
    W, H = fig.get_size_inches() / MM
    xs = [p[0] / W for p in pts_mm]
    ys = [1 - p[1] / H for p in pts_mm]
    fig.add_artist(Line2D(xs, ys, transform=fig.transFigure, **kw))



LEGEND_FONT = "Liberation Sans"   # same family as the figures
LEGEND_PT = 10.0                  # legend font size; never shrunk (cut words instead)
LEGEND_PITCH = 1.2                # baseline-to-baseline distance, x font size
LEGEND_PARA_GAP = 0.5             # extra space between paragraphs, x line pitch
LEGEND_MARGIN_MM = 8.0            # page margin on every side
LEGEND_SUP = (0.7, 0.36)          # superscript size and baseline rise, x font size
LEGEND_ASC, LEGEND_DESC = 0.905, 0.212   # Liberation Sans ascent / descent (em), for the page height
LEGEND_MAX_MM = 297.0             # warn when a legend page is taller than A4
LEGEND_ITAL_CORR = 0.05           # gap after italic before upright text (em); < pdftotext's word-break gap


def _legend_prop(italic, bold, sup):
    from matplotlib.font_manager import FontProperties
    size = LEGEND_PT * (LEGEND_SUP[0] if sup else 1.0)
    return FontProperties(family=LEGEND_FONT, style="italic" if italic else "normal",
                          weight="bold" if bold else "normal", size=size)


_FONT_CACHE = {}


def _legend_font(style):
    """FT2Font exactly as the PDF backend loads it (same file lookup, size set at 72 dpi -> points)."""
    from matplotlib import font_manager as fm
    if style not in _FONT_CACHE:
        prop = _legend_prop(*style)
        try:
            files = fm.fontManager._find_fonts_by_props(prop)      # what backend_pdf uses
        except AttributeError:                                   # older matplotlib
            files = fm.findfont(prop, fallback_to_default=False)
        font = fm.get_font(files)
        if font.family_name != LEGEND_FONT:
            raise SystemExit(f"write_figure: font '{LEGEND_FONT}' not found (got '{font.family_name}')")
        _FONT_CACHE[style] = (font, prop)
    font, prop = _FONT_CACHE[style]
    font.clear()
    font.set_size(prop.get_size_in_points(), 72)
    return font, prop


def _advance_pt(text, style):
    """Advance width (pt) of one same-style text run, measured as backend_pdf places the glyphs."""
    from matplotlib import _text_helpers
    try:
        from matplotlib.ft2font import Kerning, LoadFlags
        kern, noh = Kerning.UNFITTED, LoadFlags.NO_HINTING
    except ImportError:                                          # matplotlib < 3.10
        from matplotlib.ft2font import KERNING_UNFITTED as kern, LOAD_NO_HINTING as noh
    font, _ = _legend_font(style)
    last = None
    for last in _text_helpers.layout(text, font, kern_mode=kern):
        pass
    if last is None:
        return 0.0
    g = last.ft_object.load_glyph(last.glyph_idx, flags=noh)
    return last.x + g.linearHoriAdvance / 65536


def _parse_markup(par):
    """One paragraph of markup -> list of (char, (italic, bold, sup)).
    *italic*, **bold**, ^superscript^, backslash escapes the next character. Line breaks are spaces."""
    heading = par.startswith("## ")
    s = re.sub(r"\s+", " ", par[3:] if heading else par).strip()
    out, it, bd, sup, i = [], False, heading, False, 0
    while i < len(s):
        c = s[i]
        if c == "\\" and i + 1 < len(s):
            out.append((s[i + 1], (it, bd, sup)))
            i += 2
            continue
        if s.startswith("**", i):
            bd = not bd
            i += 2
            continue
        if c == "*":
            it = not it
        elif c == "^":
            sup = not sup
        else:
            out.append((c, (it, bd, sup)))
        i += 1
    if it or sup or bd != heading:
        raise ValueError(f"write_figure: unclosed *, ** or ^ markup in: {par[:70]}...")
    # no-break spaces: italic 'P.' + species; relations ('≥ 5', 'P = 0.14'); number + unit ('5 kb'); 'Fig. S4'
    text = "".join(c for c, _ in out)
    for m in re.finditer(" ", text):
        j = m.start()
        prev = text[:j].rsplit(" ", 1)[-1]
        nxt = text[j + 1:].split(" ", 1)[0]
        if (out[j][1][0] and re.fullmatch(r"\(?[A-Z]\.", prev)) or prev in _NB_AFTER or nxt in _NB_BEFORE \
                or (prev[-1:].isdigit() and nxt.rstrip(",.;:)") in _UNITS) or prev in ("Fig.", "Figs"):
            out[j] = (" ", out[j][1])
    return out


_NB_AFTER = {"≥", "≤", "=", "<", ">", "~", "≈", "±", "×", "n", "P"}   # never break after these
_NB_BEFORE = {"=", "≈", "±"}                                     # never break before these
_UNITS = {"bp", "kb", "Mb", "Gb", "kyr", "Myr", "yr", "mm", "cm", "%", "×", "pt", "kb)", "Mb)"}


def _check_glyphs(chars):
    missing = {}
    for c, st in chars:
        if c in "  ":
            continue
        font, _ = _legend_font(st)
        if font.get_char_index(ord(c)) == 0:
            missing.setdefault(c, set()).add(font.style_name)
    if missing:
        m = "; ".join(f"'{c}' (U+{ord(c):04X}, {', '.join(sorted(v))})" for c, v in missing.items())
        raise ValueError(f"write_figure: {LEGEND_FONT} lacks {m}. For exponents write ^…^ markup, "
                         f"e.g. 3×10^−8^ (U+2212 minus).")


def _runs(chars):
    """Merge consecutive same-style characters into (text, style) runs (no-break spaces drawn as spaces)."""
    runs = []
    for c, st in chars:
        c = " " if c == " " else c
        if runs and runs[-1][1] == st:
            runs[-1][0] += c
        else:
            runs.append([c, st])
    return [(t, st) for t, st in runs]


def _width(chars):
    runs = _runs(chars)
    corr = sum(1 for (_, a), (t, b) in zip(runs, runs[1:]) if a[0] and not b[0] and not t.startswith(" "))
    return sum(_advance_pt(t, st) for t, st in runs) + corr * LEGEND_ITAL_CORR * LEGEND_PT


def _tokens(chars, avail):
    """Split into words at spaces: (space char or None, [chars]); a word wider than the line is broken after
    '/' (paths) or, failing that, between characters."""
    toks, cur, sp = [], [], None
    for c, st in chars + [(" ", (False, False, False))]:
        if c == " ":
            if cur:
                toks.append((sp, cur))
            cur, sp = [], (c, st)
        else:
            cur.append((c, st))
    out = []
    for sp, w in toks:
        if _width(w) <= avail:
            out.append((sp, w))
            continue
        piece = []
        for ch in w:
            if piece and _width(piece + [ch]) > avail:
                cut = max((i for i, c in enumerate(piece) if c[0] == "/"), default=len(piece) - 1)
                out.append((sp, piece[:cut + 1]))
                sp, piece = None, piece[cut + 1:]
            piece.append(ch)
        out.append((sp, piece))
    return out


def _wrap(chars, avail):
    """Greedy ragged-right line breaking with exact per-run measurement. Returns lists of chars per line."""
    lines, line = [], []
    for sp, w in _tokens(chars, avail):
        cand = line + ([sp] if (line and sp) else []) + w
        if line and _width(cand) > avail:
            lines.append(line)
            line = list(w)
        else:
            line = cand
    if line:
        lines.append(line)
    return lines


def _legend_page(text, width_in):
    """Lay out one legend page (markup text) at the given width; height fitted to the text.
    Returns (Figure, height_mm, words per paragraph)."""
    from matplotlib.figure import Figure
    pts = width_in * 72
    margin = LEGEND_MARGIN_MM / 25.4 * 72
    avail = pts - 2 * margin
    pitch = LEGEND_PT * LEGEND_PITCH
    pars = [p.strip() for p in re.split(r"\n\s*\n", text.strip()) if p.strip()]
    placed, words, y = [], [], margin + LEGEND_ASC * LEGEND_PT
    for k, par in enumerate(pars):
        chars = _parse_markup(par)
        _check_glyphs(chars)
        words.append(len("".join(c for c, _ in chars).replace(" ", " ").split()))
        if k:
            y += pitch * (1 + LEGEND_PARA_GAP)
        for j, line in enumerate(_wrap(chars, avail)):
            if j:
                y += pitch
            x, prev = margin, None
            for t, st in _runs(line):
                if prev is not None and prev[0] and not st[0] and not t.startswith(" "):
                    x += LEGEND_ITAL_CORR * LEGEND_PT          # italic correction before upright punctuation
                placed.append((x, y, t, st))
                x += _advance_pt(t, st)
                prev = st
    h_pt = y + LEGEND_DESC * LEGEND_PT + margin
    fig = Figure(figsize=(width_in, h_pt / 72))
    for x, yy, t, st in placed:
        _, prop = _legend_font(st)
        rise = LEGEND_SUP[1] * LEGEND_PT if st[2] else 0.0
        fig.text(x / pts, 1 - (yy - rise) / h_pt, t, fontproperties=prop, ha="left", va="baseline",
                 color=INK, parse_math=False)
    return fig, h_pt / 72 * 25.4, words


def write_figure(stem, figs, short, long, outdir, analysis_dir, script, date=None):
    """Write <stem>.pdf (figure page(s), then the legend page; vector, fonttype 42) and the figure-only
    PNG(s) at 300 dpi (<stem>.png, or <stem>_p<k>.png for several figure pages). No SVG.

    figs: a matplotlib Figure or a list of them (the figure pages). short / long: legend text in a small markup:
    *italic*, **bold**, ^superscript^ (e.g. 3×10^−8^), a blank line starts a new paragraph, a paragraph starting
    with '## ' is a bold heading line, backslash escapes a markup character. Pass long=None: the PDF has one legend page, `short`, ending in the line
    'Analysis: <host>:<analysis_dir> · Script: <script> · <date>' (date defaults to today, ISO). Older scripts
    that pass a long legend still get it as its own page, which then carries that last line.
    Legend pages are as wide as the first figure page, 10-pt Liberation Sans, ragged right, height fitted to
    the text. Only writes files; deletes nothing. Returns a dict: paths, page sizes (mm), word counts
    (long_* are None without a long legend)."""
    import datetime
    import socket
    import matplotlib as mpl
    from matplotlib.backends.backend_pdf import PdfPages
    figs = list(figs) if isinstance(figs, (list, tuple)) else [figs]
    if not figs:
        raise ValueError("write_figure: no figure pages")
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    date = date or datetime.date.today().isoformat()
    esc = lambda s: re.sub(r"([\\*^])", r"\\\1", str(s))    # noqa: E731
    prov = (f"Analysis: {esc(socket.gethostname())}:{esc(analysis_dir)} · "
            f"Script: {esc(script)} · {esc(date)}")
    has_long = bool(long and long.strip())
    texts = [short, long.rstrip() + "\n\n" + prov] if has_long else [short.rstrip() + "\n\n" + prov]
    width_in = figs[0].get_size_inches()[0]
    if any(abs(f.get_size_inches()[0] - width_in) > 1e-6 for f in figs):
        log(f"write_figure: {stem}: figure pages differ in width; legend pages use the first ({width_in * 25.4:.0f} mm)")
    with mpl.rc_context({"pdf.fonttype": 42, "ps.fonttype": 42}):
        pages = [_legend_page(t, width_in) for t in texts]
        pngs = []
        for k, f in enumerate(figs, 1):
            p = outdir / (f"{stem}.png" if len(figs) == 1 else f"{stem}_p{k}.png")
            f.savefig(p, dpi=300)
            pngs.append(str(p))
        pdf = outdir / f"{stem}.pdf"
        with PdfPages(pdf, metadata={"Title": stem}) as pp:
            for f in figs:
                pp.savefig(f, dpi=300)
            for f, _, _ in pages:
                pp.savefig(f)
    for name, (_, h, _) in zip(("short", "long"), pages):
        if h > LEGEND_MAX_MM:
            log(f"write_figure: {stem}: {name} legend page is {h:.0f} mm tall (> {LEGEND_MAX_MM:.0f} mm); cut words")
    short_par = pages[0][2] if has_long else pages[0][2][:-1]    # without the provenance line
    out = dict(pdf=str(pdf), png=pngs, n_pages=len(figs) + len(pages),
               fig_mm=[tuple(round(v * 25.4, 1) for v in f.get_size_inches()) for f in figs],
               short_mm=round(pages[0][1], 1), short_words=sum(short_par), short_par_words=short_par,
               long_mm=None, long_words=None, long_par_words=None)
    if has_long:
        out.update(long_mm=round(pages[1][1], 1), long_words=sum(pages[1][2][:-1]),
                   long_par_words=pages[1][2][:-1])
    return out
