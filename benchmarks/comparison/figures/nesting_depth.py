#!/usr/bin/env python3
"""nesting_depth — LTR_peel README figure: what the depth0/1/2 FASTA records hold and what deleting the masked letters leaves.

One real nest of three LTR retrotransposons on A. thaliana chromosome 4 (whole-genome LTR_peel run of the benchmark),
drawn in the visual language of the README workflow figure (tiles, letter headers, element glyph, pills, identity colors):
  A  each element's FASTA record as written, drawn at its genome position, inserts stacked beneath their hosts (IGV expanded);
     every masked run (N, R) is drawn where it lies in the record, at true proportions, in the hue of the element it stands for
  B  the same records after deleting every non-ACGT letter (what `sed '/^>/!s/[^ATCGatcg]//g'` does), same bp scale
  C  the letters of the depth2 record across its insert, before and after deletion
Protein domains come from the run's GFF3 (protein_match under each element); TSDs from the depth tables' tsd column.
Inputs (read-only, under $BENCH_ROOT): runs/arab/ltr_peel/arab_LTRs_depth{0,1,2}_clean_ltr.{tsv,fa},
runs/arab/ltr_peel/arab_LTRs_all_depth_LTR_cleaned.gff3, data/arab/genome.fa (+ .fai); record coordinates are checked
against the genome first. The drawn data are written to source_data/nesting_depth/ and the figure is drawn from those
tables as read back; --from-source-data redraws from them alone (no $BENCH_ROOT; the checks against the raw files are skipped).
Outputs: docs/images/nesting_depth.png (300 dpi), nesting_depth.pdf (figure + legend page) and nesting_depth_numbers.tsv
beside this script.
Usage: nesting_depth.py [-v] [--from-source-data]
"""
import argparse
import glob
import os
import re
import shutil
import sys
from pathlib import Path

sys.dont_write_bytecode = True                       # no __pycache__ next to the imported poafig.py

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import pandas as pd  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # poafig.py sits beside this script
import poafig  # noqa: E402
import matplotlib as mpl  # noqa: E402,F401
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
from matplotlib.patches import Circle, FancyBboxPatch, Polygon, Rectangle  # noqa: E402

STEM = "nesting_depth"
HERE = Path(__file__).resolve().parent
TMP = HERE / "tmp" / STEM                              # write_figure output before the PNG and PDF are moved
SD = HERE / "source_data" / STEM                       # drawn data, the only input of the drawing code
PNG_OUT = HERE.parents[2] / "docs" / "images" / f"{STEM}.png"
PFX = "arab_LTRs"
RUN_REL = Path("runs/arab/ltr_peel")                   # paths under $BENCH_ROOT, as cited in the numbers table
GFF_REL = RUN_REL / f"{PFX}_all_depth_LTR_cleaned.gff3"
GENOME_REL = Path("data/arab/genome.fa")
NAMES_REL = Path("data/arab/names.tsv")
RUN = GFF = GENOME = None                              # set from $BENCH_ROOT in full mode
CHROM, CHROM_NAME = "c0004", "chromosome 4"          # data/arab/names.tsv: c0004 = NC_003075.7
# the nest: host, its insert, the insert's insert
NEST = [("H", 2, "c0004:4407852-4436793"), ("C", 1, "c0004:4416804-4434093"), ("G", 0, "c0004:4425022-4431463")]
SED_STRIP = re.compile(r"[^ATCGatcg]")                 # what sed '/^>/!s/[^ATCGatcg]//g' deletes from sequence lines
MASK_DEPTH = {"N": 0, "R": 1, "D": 2, "Y": 3}         # mask letter = the masked element's own depth (docs/outputs.md 4.2)

# ---- the README workflow figure's visual language
INK, MUTED, GRAY = "#1A1A1A", "#5A5A5A", "#8C8C8C"
TILE_C, TILE_R = "#F3F3F3", 2.2                         # tile fill, corner rounding (mm)
PT = 0.3528                                            # mm per point
# one identity color per element (dark = LTR, light = internal region, mid = pill and the block that masks it)
EL = {2: dict(ltr="#8A3D00", int="#EE9A5E", mid="#D55E00"),
      1: dict(ltr="#8E4A72", int="#EBC2D9", mid="#CC79A7"),
      0: dict(ltr="#00684C", int="#9AD9C6", mid="#009E73")}
RIBBON = "#F5DFD0"                                     # panel C ribbons: the depth2 host's own bases, pale vermillion

# ---- layout (mm, from the top-left corner)
W_MM = 180.0
TX0, TX1 = 0.8, 179.2                  # tiles' outer edges
GUT = 4.0                              # gutter between tiles (workflow: 4 mm)
T1Y0, T1Y1 = 0.8, 57.6                 # top row of tiles (A, B)
T2Y0, T2Y1 = 60.8, 93.2                # bottom tile (C)
H_MM = T2Y1 + 0.8
HDY = 5.2                              # header line below a tile's top (workflow: 6.0 - 0.8)
PILL_X = TX0 + 1.8                     # pills align with the header circle's left edge
A_X0 = 27.5                            # panel A data start (= genome window start), also panel C sequence start
A_LO, A_HI = 4_407_000, 4_437_400      # genome window drawn (0-based edges)
S = 2.5                                # mm per kb, panels A and B (same bp scale)
ROW_Y = {2: 17.6, 1: 29.6, 0: 41.6}    # row centers by depth
EH = 4.6                               # element height (LTR blocks, masks); internal region 0.6 of it (workflow glyph)
TSD_S = 1.4                            # TSD square side
AXIS_Y = 48.2
C_Y1, C_Y2 = 78.6, 88.6                # panel C text rows (vertical centers)
SEQ_PT = 9.0                           # sequence text
MONO_ADV = 0.60009 * SEQ_PT / 72 * 25.4    # Liberation Mono advance (1229/2048 em) in mm
DOTS = 2                               # mono cells taken by an elision mark (three gray dots)
FS_NAME, FS_PILL, FS_NUM, FS_LAB, FS_TICK = 10.5, 7.5, 7.5, 7.0, 6.5


def log(msg, verbose=True):
    if verbose:
        print(msg, file=sys.stderr, flush=True)


# ---------------------------------------------------------------- data
def read_fai(path):
    fai = {}
    for line in open(str(path) + ".fai"):
        n, ln, off, lb, lw = line.split("\t")[:5]
        fai[n] = (int(ln), int(off), int(lb), int(lw))
    return fai


def fetch(chrom, start0, end0):
    """Genome slice, 0-based half-open, read through the .fai offsets (no full load)."""
    ln, off, lb, lw = read_fai(GENOME)[chrom]
    start0, end0 = max(0, start0), min(ln, end0)
    a = off + (start0 // lb) * lw + start0 % lb
    b = off + (end0 // lb) * lw + end0 % lb
    with open(GENOME, "rb") as fh:
        fh.seek(a)
        return fh.read(b - a).decode().replace("\n", "").replace("\r", "").upper()


def read_record(fa, key):
    """One FASTA record whose header starts with key + '#'; returns (header, sequence as written)."""
    hdr, buf, hit = None, [], False
    with open(fa) as fh:
        for line in fh:
            if line.startswith(">"):
                if hit:
                    break
                hit = line[1:].split("#")[0] == key
                if hit:
                    hdr = line[1:].rstrip("\n")
            elif hit:
                buf.append(line.strip())
    if hdr is None:
        raise SystemExit(f"{STEM}: record {key} not found in {fa}")
    return hdr, "".join(buf)


def runs(seq):
    """Runs of one letter other than ACGT, as (letter, start0, end0) in record coordinates."""
    return [(m.group()[0], m.start(), m.end()) for m in re.finditer(r"([^ACGT])\1*", seq.upper())]


def read_gff(keys):
    """Each element's GFF3 features (streamed): {key: {'id', 'ltr': [(s, e)], 'dom': [(name, s, e)]}}, 1-based closed."""
    out, id2key = {}, {}
    with open(GFF) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            att = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            if f[2] == "LTR_retrotransposon" and att.get("Name") in keys:
                if f[0] != CHROM or f"{f[0]}:{f[3]}-{f[4]}" != att["Name"]:
                    raise SystemExit(f"{STEM}: GFF3 element {att['Name']} sits at {f[0]}:{f[3]}-{f[4]}")
                id2key[att["ID"]] = att["Name"]
                out[att["Name"]] = dict(id=att["ID"], ltr=[], dom=[], strand=f[6])
            elif att.get("Parent") in id2key:
                el = out[id2key[att["Parent"]]]
                if f[2] == "long_terminal_repeat":
                    el["ltr"].append((int(f[3]), int(f[4])))
                elif f[2] == "protein_match":
                    el["dom"].append((att["Name"], int(f[3]), int(f[4])))
    missing = [k for k in keys if k not in out]
    if missing:
        raise SystemExit(f"{STEM}: not in {GFF}: {missing}")
    return out


def load(verbose):
    els = {}
    gff = read_gff([k for _, _, k in NEST])
    for name, depth, key in NEST:
        tsv = RUN / f"{PFX}_depth{depth}_clean_ltr.tsv"
        fa = RUN / f"{PFX}_depth{depth}_clean_ltr.fa"
        t = pd.read_csv(tsv, sep="\t").rename(columns={"#seq_id": "seq_id"})
        row = t[t.seq_id.str.split("#").str[0] == key]
        if len(row) != 1:
            raise SystemExit(f"{STEM}: {key} not found once in {tsv}")
        r = row.iloc[0]
        hdr, seq = read_record(fa, key)
        chrom, s, e = re.match(r"^([^:]+):(\d+)-(\d+)$", key).groups()
        s, e = int(s), int(e)
        el = dict(name=name, depth=depth, key=key, chrom=chrom, s=s, e=e, header=hdr, seq=seq, tsv=str(tsv), fa=str(fa),
                  seq_len=int(r.seq_len), ltr5_end=int(r.ltr5_end), ltr3_start=int(r.ltr3_start),
                  ltr5_len=int(r.ltr5_len), ltr3_len=int(r.ltr3_len), confidence=r.confidence, motif=r.motif,
                  tsd=str(r.tsd), tsd_offset=str(r.tsd_offset), domains_col=str(r.domains),
                  cls=r.seq_id.split("#")[1], orientation=r.orientation, nest_status=r.nest_status, k2p=float(r.k2p))
        el["runs"] = runs(seq)
        el["stripped"] = SED_STRIP.sub("", seq)
        el["gff"] = gff[key]
        els[name] = el
    check(els, verbose)
    return els


def check(els, verbose):
    """Structural and format checks on the drawn records; exits before drawing if any fails."""
    H, C, G = els["H"], els["C"], els["G"]
    for el in els.values():
        if len(el["seq"]) != el["seq_len"] or el["seq_len"] != el["e"] - el["s"] + 1:
            raise SystemExit(f"{STEM}: {el['key']}: record length != seq_len != end-start+1")
        if el["orientation"] != "+":
            raise SystemExit(f"{STEM}: {el['key']}: drawn as stored forward; orientation is {el['orientation']}")
        g = fetch(el["chrom"], el["s"] - 1, el["e"])          # 1-based closed -> 0-based half-open
        bad = sum(1 for a, b in zip(el["seq"].upper(), g) if a in "ACGT" and a != b)
        if bad or g.count("N"):
            raise SystemExit(f"{STEM}: {el['key']}: {bad} unmasked bases differ from the genome, {g.count('N')} genome N")
        if not (el["seq"][:2].upper() == "TG" and el["seq"][-2:].upper() == "CA"):
            raise SystemExit(f"{STEM}: {el['key']}: ends are not TG..CA")
        el["tsd5"], el["tsd3"] = fetch(el["chrom"], el["s"] - 6, el["s"] - 1), fetch(el["chrom"], el["e"], el["e"] + 5)
        # TSD: drawn only where the depth table lists one; the listed TSD must sit in the genome at both ends,
        # shifted by tsd_offset (o5 bp into the element at the 5' end, -o3 bp into the flank at the 3' end)
        if el["tsd"] not in (".", "", "nan"):
            o5, o3 = (int(v) for v in el["tsd_offset"].split(","))
            n = len(el["tsd"])
            t5 = fetch(el["chrom"], el["s"] - 1 - n + o5, el["s"] - 1 + o5)
            t3 = fetch(el["chrom"], el["e"] - o3, el["e"] - o3 + n)
            if not (t5 == el["tsd"].upper() == t3):
                raise SystemExit(f"{STEM}: {el['key']}: table TSD {el['tsd']} (offset {el['tsd_offset']}) not found "
                                 f"in the genome: {t5}/{t3}")
            el["has_tsd"] = True
        else:
            el["has_tsd"] = False
        # GFF3 agrees with the table: LTR coordinates, strand, and the domain list
        gf = el["gff"]
        ltrs = sorted(gf["ltr"])
        exp_ltr = [(el["s"], el["s"] + el["ltr5_end"] - 1), (el["s"] + el["ltr3_start"] - 1, el["e"])]
        if ltrs != exp_ltr or gf["strand"] != el["orientation"]:
            raise SystemExit(f"{STEM}: {el['key']}: GFF3 LTRs {ltrs} / strand {gf['strand']} differ from the table "
                             f"{exp_ltr} / {el['orientation']}")
        col = sorted((m.group(1), int(m.group(2)), int(m.group(3)))
                     for m in re.finditer(r"([A-Za-z]+)\|[^@;]*@(\d+)-(\d+)", el["domains_col"]))
        if sorted(gf["dom"]) != col:
            raise SystemExit(f"{STEM}: {el['key']}: GFF3 domains {sorted(gf['dom'])} differ from the table {col}")
        # domains in record coordinates (0-based half-open), clipped to unmasked letters; nothing drawn in a mask
        el["dom"] = []
        for nm, ds, de in sorted(gf["dom"], key=lambda d: d[1]):
            a, b = ds - el["s"], de - el["s"] + 1
            if not (el["ltr5_end"] <= a and b <= el["ltr3_start"] - 1):
                raise SystemExit(f"{STEM}: {el['key']}: domain {nm} {ds}-{de} is not inside the internal region")
            pieces = [(a, b)]
            for _, ra, rb in el["runs"]:
                pieces = [q for p0, p1 in pieces for q in ((p0, min(p1, ra)), (max(p0, rb), p1)) if q[1] > q[0]]
            el["dom"] += [(nm, p0, p1, ds, de) for p0, p1 in pieces]
        log(f"  {el['name']} depth{el['depth']} {el['key']} {el['seq_len']} bp, {el['confidence']}, "
            f"flanks {el['tsd5']}/{el['tsd3']}, table TSD {el['tsd']} ({el['tsd_offset']}), "
            f"runs {[(c, b - a) for c, a, b in el['runs']]}, stripped {len(el['stripped'])} bp, "
            f"domains {[(d[0], d[2] - d[1]) for d in el['dom']]}", verbose)
    # the TSD facts the legend states: depth2 and depth0 have one, depth1 has none
    if not (H["has_tsd"] and G["has_tsd"] and not C["has_tsd"]):
        raise SystemExit(f"{STEM}: TSDs differ from the legend: H {H['tsd']}, C {C['tsd']}, G {G['tsd']}")
    # nest links resolve both ways
    links = (f"nest-outer:{C['key']}" in H["nest_status"] and f"nest-outer:{G['key']}" in H["nest_status"]
             and f"nest-outer:{G['key']}" in C["nest_status"] and f"nest-inner:{H['key']}" in C["nest_status"]
             and f"nest-inner:{C['key']}" in G["nest_status"] and f"nest-inner:{H['key']}" in G["nest_status"])
    if not links:
        raise SystemExit(f"{STEM}: nest_status links do not resolve")
    # every masked run is exactly an insert listed in the tables (no assembly gaps, no unlisted masks)
    exp_H = [("R", C["s"] - H["s"], G["s"] - H["s"]), ("N", G["s"] - H["s"], G["e"] - H["s"] + 1),
             ("R", G["e"] - H["s"] + 1, C["e"] - H["s"] + 1)]
    exp_C = [("N", G["s"] - C["s"], G["e"] - C["s"] + 1)]
    if H["runs"] != exp_H or C["runs"] != exp_C or G["runs"]:
        raise SystemExit(f"{STEM}: masked runs differ from the inserts: {H['runs']} {C['runs']} {G['runs']}")
    # each mask letter stands for the element of that depth (so it takes that element's color)
    for el in (H, C):
        for c, a, b in el["runs"]:
            ins = G if (a, b) == (G["s"] - el["s"], G["e"] - el["s"] + 1) else C     # else a piece of C
            if MASK_DEPTH[c] != ins["depth"]:
                raise SystemExit(f"{STEM}: {el['key']}: run of {c} stands for depth{ins['depth']} element")
    for el in els.values():
        if len(el["stripped"]) != el["seq_len"] - sum(b - a for _, a, b in el["runs"]):
            raise SystemExit(f"{STEM}: {el['key']}: stripped length mismatch")
    # the legend calls all three Athila elements
    if any(el["cls"] != "LTR/Gypsy/Athila" for el in els.values()):
        raise SystemExit(f"{STEM}: legend says three Athila elements; classes are {[el['cls'] for el in els.values()]}")
    # each insert lies inside its host's internal region (not in an LTR)
    for host, ins in ((H, C), (C, G)):
        if not (host["s"] + host["ltr5_end"] <= ins["s"] and ins["e"] < host["s"] + host["ltr3_start"] - 1):
            raise SystemExit(f"{STEM}: {ins['key']} is not inside the internal region of {host['key']}")


# ---------------------------------------------------------------- source data
C_HOST, C_RUN_END = 12, 6               # panel C: host bases shown per side, letters shown at each end of a run
SD_READ = dict(sep="\t", keep_default_na=False)       # '' stays '' (elision marks hold no letters)


def fa_rel(depth):
    return str(RUN_REL / f"{PFX}_depth{depth}_clean_ltr.fa")


def tsv_rel(depth):
    return str(RUN_REL / f"{PFX}_depth{depth}_clean_ltr.tsv")


def c_letters(H):
    """Panel C letters of the depth2 record, piece by piece: the host bases on each side of the insert, each masked run
    shortened to its first and last letters, elision marks (no letters) for what is left out. Record positions 1-based."""
    seq, (r1, r2, r3) = H["seq"], H["runs"]
    k, m = C_HOST, C_RUN_END
    rows = []

    def add(line, kind, run, a, b, show=True):         # record edges a..b, 0-based half-open
        rows.append(dict(line=line, kind=kind, run=run, letters=seq[a:b] if show else "", rec_start=a + 1, rec_end=b))

    for line in ("as_written", "after_deleting"):
        add(line, "elision", 0, 0, r1[1] - k, show=False)
        add(line, "host", 0, r1[1] - k, r1[1])
        if line == "as_written":
            for i, (c, a, b) in enumerate((r1, r2, r3), 1):
                add(line, "masked", i, a, a + m)
                add(line, "elision", i, a + m, b - m, show=False)
                add(line, "masked", i, b - m, b)
        add(line, "host", 0, r3[2], r3[2] + k)
        add(line, "elision", 0, r3[2] + k, H["seq_len"], show=False)
    df = pd.DataFrame(rows)
    df.insert(1, "piece", df.groupby("line").cumcount() + 1)
    return df


def write_source_data(els):
    """The drawn data as tidy TSVs in SD: element records, masked runs, protein domains, panel C letters."""
    order = [els[n] for n, _, _ in NEST]
    elements = pd.DataFrame([dict(
        element=el["name"], depth=el["depth"], key=el["key"], header=el["header"], chrom=el["chrom"],
        start=el["s"], end=el["e"], seq_len=el["seq_len"], ltr5_end=el["ltr5_end"], ltr3_start=el["ltr3_start"],
        ltr5_len=el["ltr5_len"], ltr3_len=el["ltr3_len"], tsd=el["tsd"], tsd_offset=el["tsd_offset"] if el["has_tsd"] else ".",
        has_tsd=el["has_tsd"], family_class=el["cls"], len_after_deleting=len(el["stripped"])) for el in order])
    masked = pd.DataFrame([dict(element=el["name"], depth=el["depth"], letter=c, rec_start=a + 1, rec_end=b)
                           for el in order for c, a, b in el["runs"]])
    domains = pd.DataFrame([dict(element=el["name"], depth=el["depth"], domain=nm, rec_start=a + 1, rec_end=b,
                                 genome_start=ds, genome_end=de) for el in order for nm, a, b, ds, de in el["dom"]])
    SD.mkdir(parents=True, exist_ok=True)
    for name, df in (("AB_elements", elements), ("AB_masked_runs", masked), ("AB_domains", domains),
                     ("C_letters", c_letters(els["H"]))):
        df.to_csv(SD / f"{name}.tsv", sep="\t", index=False)


def read_source_data():
    """Element dicts and the panel C letters, from the source-data tables only (record edges back to 0-based)."""
    el_t = pd.read_csv(SD / "AB_elements.tsv", **SD_READ)
    runs_t = pd.read_csv(SD / "AB_masked_runs.tsv", **SD_READ)
    dom_t = pd.read_csv(SD / "AB_domains.tsv", **SD_READ)
    letters = pd.read_csv(SD / "C_letters.tsv", **SD_READ, dtype={"letters": str})
    els = {}
    for r in el_t.itertuples(index=False):
        el = dict(name=r.element, depth=int(r.depth), key=r.key, header=r.header, chrom=r.chrom, s=int(r.start),
                  e=int(r.end), seq_len=int(r.seq_len), ltr5_end=int(r.ltr5_end), ltr3_start=int(r.ltr3_start),
                  ltr5_len=int(r.ltr5_len), ltr3_len=int(r.ltr3_len), tsd=str(r.tsd), tsd_offset=str(r.tsd_offset),
                  has_tsd=str(r.has_tsd) == "True", cls=r.family_class, len_stripped=int(r.len_after_deleting),
                  fa=fa_rel(int(r.depth)), tsv=tsv_rel(int(r.depth)))
        el["runs"] = sorted(((x.letter, int(x.rec_start) - 1, int(x.rec_end))
                             for x in runs_t[runs_t.element == r.element].itertuples(index=False)), key=lambda t: t[1])
        el["dom"] = [(x.domain, int(x.rec_start) - 1, int(x.rec_end), int(x.genome_start), int(x.genome_end))
                     for x in dom_t[dom_t.element == r.element].itertuples(index=False)]
        els[r.element] = el
    return els, letters


# ---------------------------------------------------------------- drawing helpers
def fig_axes():
    fig = plt.figure(figsize=(W_MM * poafig.MM, H_MM * poafig.MM))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W_MM)
    ax.set_ylim(H_MM, 0)
    ax.axis("off")
    return fig, ax


def gx(p0):
    """0-based genome edge -> x (mm), panel A. 1-based base p spans edges p-1 .. p."""
    return A_X0 + (p0 - A_LO) / 1000 * S


def segments(el):
    """The record as drawn segments in record coordinates: (kind, a, b), 0-based half-open."""
    marks = [("ltr", 0, el["ltr5_end"]), ("ltr", el["ltr3_start"] - 1, el["seq_len"])]
    out, cuts = [], sorted({0, el["ltr5_end"], el["ltr3_start"] - 1, el["seq_len"]}
                           | {a for _, a, _ in el["runs"]} | {b for _, _, b in el["runs"]})
    for a, b in zip(cuts[:-1], cuts[1:]):
        kind = "int"
        for _, x0, x1 in marks:
            if x0 <= a and b <= x1:
                kind = "ltr"
        for c, x0, x1 in el["runs"]:
            if x0 <= a and b <= x1:
                kind = c
        out.append((kind, a, b))
    return out


def stripped_pos(el, p):
    """Record edge p (0-based) -> edge in the record after deleting the masked letters."""
    return p - sum(min(max(p - a, 0), b - a) for _, a, b in el["runs"])


def rect(ax, x0, x1, y, h, fc, z=3, ec="none", lw=0.0):
    ax.add_patch(Rectangle((x0, y - h / 2), x1 - x0, h, fc=fc, ec=ec, lw=lw, zorder=z))


def tsd(ax, x, y, s=TSD_S, z=6):
    """Open square (target-site duplication), centred at x (workflow glyph)."""
    ax.add_patch(Rectangle((x - s / 2, y - s / 2), s, s, fc="white", ec=INK, lw=0.8, zorder=z))


def draw_element(ax, x_of, el, segs, doms, y, letters=True):
    """Workflow element glyph from real coordinates: LTRs full height in the element's dark shade, internal region
    0.6 high in its light shade, masked runs full height in the mid hue of the element they stand for (white letter),
    black domain bars, open TSD squares just outside both ends when the table lists a TSD."""
    col = EL[el["depth"]]
    ih = EH * 0.6
    for kind, a, b in segs:
        x0, x1 = x_of(a), x_of(b)
        if kind == "ltr":
            rect(ax, x0, x1, y, EH, col["ltr"], z=5)
        elif kind == "int":
            rect(ax, x0, x1, y, ih, col["int"], z=3)
        else:
            rect(ax, x0, x1, y, EH, EL[MASK_DEPTH[kind]]["mid"], z=4)
            if letters and x1 - x0 > 3.0:
                ax.text((x0 + x1) / 2, y + 0.15, kind, ha="center", va="center", fontsize=8, fontweight="bold",
                        color="white", zorder=6)
    for _, a, b in doms:
        rect(ax, x_of(a), x_of(b), y, ih * 0.62, INK, z=4)
    xl, xr = x_of(segs[0][1]), x_of(segs[-1][2])
    if el["has_tsd"]:
        tsd(ax, xl - TSD_S / 2 - 0.35, y)
        tsd(ax, xr + TSD_S / 2 + 0.35, y)
    return xr + (TSD_S + 0.35 if el["has_tsd"] else 0.0)


def text_w(ax, s, fs, weight="normal"):
    """Width in mm of a text string at font size fs."""
    fig = ax.figure
    t = ax.text(0, 0, s, fontsize=fs, fontweight=weight)
    bb = t.get_window_extent(renderer=fig.canvas.get_renderer())
    t.remove()
    return bb.width / fig.dpi * 25.4


def pill(ax, x, y, s, fc, tc="white", fs=FS_PILL, weight="bold", pad=1.3, h=None, ha="left", z=8):
    """Rounded badge with text (workflow figure); (x, y) = left-centre or centre. Returns right x."""
    w = text_w(ax, s, fs, weight) + 2 * pad
    h = fs * PT * 1.75 if h is None else h
    x0 = x if ha == "left" else x - w / 2
    ax.add_patch(FancyBboxPatch((x0, y - h / 2), w, h, boxstyle=f"round,pad=0,rounding_size={h / 2}", fc=fc,
                                ec="none", lw=0, zorder=z))
    ax.text(x0 + w / 2, y + 0.05, s, ha="center", va="center", fontsize=fs, color=tc, fontweight=weight, zorder=z + 1)
    return x0 + w


def tile(ax, x0, x1, y0, y1):
    """Rounded light-gray tile (workflow figure)."""
    ax.add_patch(FancyBboxPatch((x0, y0), x1 - x0, y1 - y0, boxstyle=f"round,pad=0,rounding_size={TILE_R}",
                                fc=TILE_C, ec="none", zorder=0))


def header(ax, x0, y0, letter, name):
    """Tile header (workflow figure, panel letter in place of the step number)."""
    cx, cy = x0 + 4.4, y0 + HDY
    ax.add_patch(Circle((cx, cy), 2.6, fc=INK, ec="none", zorder=5))
    ax.text(cx, cy + 0.1, letter, ha="center", va="center", fontsize=8, color="white", fontweight="bold", zorder=6)
    ax.text(cx + 4.1, cy, name, ha="left", va="center", fontsize=FS_NAME, fontweight="bold", color=INK)


def row_label(ax, y, depth, sub):
    """Depth pill in the element's mid hue, a muted sub-label beneath."""
    pill(ax, PILL_X, y - 1.9, f"depth{depth}", EL[depth]["mid"])
    ax.text(PILL_X, y + 2.6, sub, ha="left", va="center", fontsize=FS_LAB, color=MUTED)


def cells(pieces):
    """Width of text pieces in monospace cells; an elision mark (None) takes DOTS cells."""
    return sum(DOTS if p is None else len(p) for p, _ in pieces)


def draw_pieces(ax, x, y, pieces):
    """Monospace text pieces side by side; None draws an elision mark (three gray dots at mid-height)."""
    pos = 0
    for p, col in pieces:
        if p is None:
            for j in range(3):
                cx = x + (pos + DOTS / 2) * MONO_ADV + (j - 1) * 0.72
                ax.add_patch(Circle((cx, y + 0.1), 0.23, fc=col, ec="none", zorder=3))
            pos += DOTS
        else:
            ax.text(x + pos * MONO_ADV, y, p, ha="left", va="center", fontsize=SEQ_PT, family="Liberation Mono",
                    color=col, zorder=3)
            pos += len(p)


def down_triangle(ax, x, y_tip, color, w=1.8, h=1.6):
    ax.add_patch(Polygon([(x - w / 2, y_tip - h), (x + w / 2, y_tip - h), (x, y_tip)], closed=True, fc=color, ec="none",
                         zorder=7))


# ---------------------------------------------------------------- figure
def build(els, letters):
    """Draws from the source-data tables as read back (read_source_data); returns the figure and the numbers rows."""
    fig, ax = fig_axes()
    H, C, G = els["H"], els["C"], els["G"]
    by_depth = {el["depth"]: el for el in els.values()}
    num = []          # (item, value, source)

    # tile edges: A ends after its longest length label, B takes the rest of the row
    a_right = max(gx(by_depth[d]["e"]) + (TSD_S + 0.35 if by_depth[d]["has_tsd"] else 0) + 1.4
                  + text_w(ax, f"{by_depth[d]['seq_len']:,} bp", FS_NUM) for d in (2, 1, 0))
    XA1 = a_right + 3.0
    XB0 = XA1 + GUT
    B_X0 = XB0 + 3.6 + TSD_S + 0.35          # panel B data start (room for a left TSD square)
    tile(ax, TX0, XA1, T1Y0, T1Y1)
    tile(ax, XB0, TX1, T1Y0, T1Y1)
    tile(ax, TX0, TX1, T2Y0, T2Y1)
    header(ax, TX0, T1Y0, "A", "FASTA records")
    header(ax, XB0, T1Y0, "B", "N and R deleted")
    header(ax, TX0, T2Y0, "C", "Letters at the join")

    # panel A: records as written, at their genome positions
    sub = {2: "two layers inside", 1: "one layer inside", 0: "nothing inside"}
    ax_x = lambda el: (lambda p: gx(el["s"] - 1 + p))        # noqa: E731  record edge -> x
    for d in (2, 1, 0):
        el, y = by_depth[d], ROW_Y[d]
        row_label(ax, y, d, sub[d])
        xr = draw_element(ax, ax_x(el), el, segments(el), [(nm, a, b) for nm, a, b, _, _ in el["dom"]], y)
        ax.text(xr + 1.4, y, f"{el['seq_len']:,} bp", ha="left", va="center", fontsize=FS_NUM, color=INK)
        num.append((f"A: depth{d} record length as written, >{el['header']}", el["seq_len"], el["fa"]))
        num.append((f"A/B: depth{d} {el['key']} 5' LTR / 3' LTR length drawn (bp; record 1-{el['ltr5_end']}, "
                    f"{el['ltr3_start']}-{el['seq_len']})", f"{el['ltr5_len']}/{el['ltr3_len']}", el["tsv"]))
        for c, a, b in el["runs"]:
            num.append((f"A/C: depth{d} record {el['key']}: run of {c} at record positions {a + 1}-{b} (1-based); "
                        f"drawn in the depth{MASK_DEPTH[c]} element's hue", b - a, el["fa"]))
        for nm, a, b, ds, de in el["dom"]:
            num.append((f"A/B: depth{d} {el['key']} domain {nm} drawn, {CHROM}:{ds}-{de} (GFF3 1-based closed) = "
                        f"record {a + 1}-{b}", b - a, str(GFF_REL)))
        num.append((f"A/B: depth{d} {el['key']} TSD drawn (table tsd, tsd_offset)",
                    f"{el['tsd']} ({el['tsd_offset']})" if el["has_tsd"] else "none (tsd '.')", el["tsv"]))
    # guides: each insert's ends, from the host strip down to the insert's own strip
    for ins in (C, G):
        for pos in (ins["s"] - 1, ins["e"]):
            x = gx(pos)
            ax.plot([x, x], [ROW_Y[2] + EH / 2, ROW_Y[ins["depth"]] - EH / 2], color=GRAY, lw=0.6,
                    ls=(0, (1.2, 1.2)), zorder=1, solid_capstyle="butt")
    # genome axis
    ax.plot([gx(A_LO), gx(A_HI)], [AXIS_Y, AXIS_Y], color=GRAY, lw=0.7, solid_capstyle="butt", zorder=2)
    for t in range(4_410_000, A_HI + 1, 5_000):
        x = gx(t)
        ax.plot([x, x], [AXIS_Y, AXIS_Y + 1.0], color=GRAY, lw=0.7, zorder=2)
        ax.text(x, AXIS_Y + 1.5, f"{t / 1e6:.3f}", ha="center", va="top", fontsize=FS_TICK, color=INK)
        num.append(("A: axis tick (Mb)", f"{t / 1e6:.3f}", "axis"))
    ax.text((gx(A_LO) + gx(A_HI)) / 2, AXIS_Y + 6.4, rf"$\it{{A.\,thaliana}}$ {CHROM_NAME} (Mb)", ha="center",
            va="center", fontsize=FS_LAB, color=INK)
    num.append((f"A: genome window drawn, {CHROM_NAME} = {CHROM} (0-based edges)", f"{A_LO}-{A_HI}",
                f"{GENOME_REL}; {NAMES_REL}"))

    # panel B: after deleting every non-ACGT letter, same bp scale, left-aligned
    bx = lambda p: B_X0 + p / 1000 * S        # noqa: E731
    for d in (2, 1, 0):
        el, y = by_depth[d], ROW_Y[d]
        L = el["len_stripped"]
        l5, l3 = el["ltr5_end"], el["seq_len"] - el["ltr3_start"] + 1
        segs = [("ltr", 0, l5), ("int", l5, L - l3), ("ltr", L - l3, L)]
        doms = [(nm, stripped_pos(el, a), stripped_pos(el, b)) for nm, a, b, _, _ in el["dom"]]
        xr = draw_element(ax, bx, el, segs, doms, y, letters=False)
        if el["runs"]:
            join = el["runs"][0][1]               # kept bases before the first masked run = the join point
            c0 = el["runs"][0][0]
            down_triangle(ax, bx(join), y - EH / 2 - 0.35, EL[MASK_DEPTH[c0]]["mid"])
            num.append((f"B: depth{d} join point (kept bases before the insert); marker in the depth{MASK_DEPTH[c0]} "
                        f"element's hue", join, el["fa"]))
        ax.text(xr + 1.4, y, f"{L:,} bp", ha="left", va="center", fontsize=FS_NUM, color=INK)
        num.append((f"B: depth{d} record length after deleting non-ACGT letters, {el['key']}", L,
                    el["fa"] + " (same as sed '/^>/!s/[^ATCGatcg]//g')"))
    sb0, sb1 = B_X0, B_X0 + 5 * S
    ax.plot([sb0, sb1], [AXIS_Y, AXIS_Y], color=GRAY, lw=0.9, solid_capstyle="butt")
    ax.text(sb1 + 1.5, AXIS_Y, "5 kb", ha="left", va="center", fontsize=FS_LAB, color=INK)
    num.append(("B: scale bar (kb), same bp scale as A", 5, "axis"))
    num.append(("A/B: bp scale (mm per kb)", S, "layout"))

    # panel C: the letters of the depth2 record across its insert
    def piece(r):
        if r.kind == "elision":
            return (None, MUTED)                  # None = elision mark
        return (r.letters, INK if r.kind == "host" else EL[MASK_DEPTH[r.letters[0]]]["mid"])

    wr = letters[letters.line == "as_written"].sort_values("piece").reset_index(drop=True)
    af = letters[letters.line == "after_deleting"].sort_values("piece").reset_index(drop=True)
    hosts = wr[wr.kind == "host"]
    lh, rh = hosts.iloc[0], hosts.iloc[-1]
    left, right, k = lh.letters, rh.letters, len(lh.letters)
    m = len(wr[wr.kind == "masked"].iloc[0].letters)
    pieces = [piece(r) for r in wr.itertuples(index=False)]
    spans = []                                    # (mask letter, full run length, first cell, end cell)
    for _, g in wr[wr.run > 0].groupby("run", sort=True):
        i0, i1 = g.index.min(), g.index.max() + 1
        spans.append((g[g.kind == "masked"].letters.iloc[0][0], int(g.rec_end.max() - g.rec_start.min() + 1),
                      cells(pieces[:i0]), cells(pieces[:i1])))
    pill(ax, PILL_X, C_Y1 - 1.9, "depth2", EL[2]["mid"])
    ax.text(PILL_X, C_Y1 + 2.6, "as written", ha="left", va="center", fontsize=FS_LAB, color=MUTED)
    ax.text(PILL_X, C_Y2, "after deleting", ha="left", va="center", fontsize=FS_LAB, color=MUTED)
    xc = A_X0
    draw_pieces(ax, xc, C_Y1, pieces)
    for c, n, i0, i1 in spans:                    # full run lengths above each run
        hue = EL[MASK_DEPTH[c]]["mid"]
        x0, x1 = xc + i0 * MONO_ADV + 0.35, xc + i1 * MONO_ADV - 0.35
        yb = C_Y1 - 3.2
        ax.plot([x0, x0, x1, x1], [yb + 0.9, yb, yb, yb + 0.9], color=hue, lw=0.8, solid_capstyle="butt")
        pill(ax, (x0 + x1) / 2, yb - 2.6, f"{n:,} {c}", hue, fs=7.0, ha="center")
    after = [piece(r) for r in af.itertuples(index=False)]
    draw_pieces(ax, xc, C_Y2, after)
    jx = xc + (DOTS + k) * MONO_ADV
    down_triangle(ax, jx, C_Y2 - 2.3, EL[MASK_DEPTH[spans[0][0]]]["mid"], w=1.6, h=1.4)
    # ribbons: the left flank stays, the right flank slides left onto the join
    yt, yb2 = C_Y1 + 2.1, C_Y2 - 2.6
    lx0 = xc + DOTS * MONO_ADV
    rx0 = xc + (cells(pieces) - DOTS - k) * MONO_ADV
    ax.add_patch(Polygon([(lx0, yt), (jx, yt), (jx, yb2), (lx0, yb2)], closed=True, fc=RIBBON, ec="none", zorder=1))
    ax.add_patch(Polygon([(rx0, yt), (rx0 + k * MONO_ADV, yt), (jx + k * MONO_ADV, yb2), (jx, yb2)], closed=True,
                         fc=RIBBON, ec="none", zorder=1))
    num.append(("C: host bases shown on each side of the insert (depth2 record)", k, H["fa"]))
    num.append(("C: masked letters shown at each end of each run", m, H["fa"]))
    num.append((f"C: host bases left of the insert, record positions {lh.rec_start}-{lh.rec_end}", left, H["fa"]))
    num.append((f"C: host bases right of the insert, record positions {rh.rec_start}-{rh.rec_end}", right, H["fa"]))
    return fig, num


LEGEND = """**Every element in an LTR_peel FASTA is one record from LTR to LTR in which each element nested inside it is \
masked by letters, so deleting those letters leaves the host's own sequence.** One real nest of three *Athila* elements on \
*Arabidopsis thaliana* chromosome 4.

**A** Each element's record at its place in the genome, with its inserts stacked beneath it. An element goes to the file \
for the number of layers nested inside it (depth0, depth1, depth2), whether or not it sits inside another. Each insert \
is masked by its own depth: N if nothing is inside it, R if one element is. Each element keeps one color, also on the \
block that masks it; black bars are protein domains and open squares target-site duplications, which the depth1 element \
lacks.

**B** The same records after deleting every letter other than A, C, G and T: each keeps both LTRs and its own internal \
sequence; triangles mark where the insert was cut out. An assembly gap (N) inside an element is deleted too.

**C** The letters of the depth2 record across its insert, before and after deletion; each run is shortened to its ends, \
with its full length above."""


def main():
    global RUN, GFF, GENOME
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("-v", "--verbose", action="store_true", help="print the structural checks of the three records")
    ap.add_argument("--from-source-data", action="store_true",
                    help=f"redraw from source_data/{STEM}/ only (no BENCH_ROOT; raw-data checks skipped)")
    args = ap.parse_args()
    if args.from_source_data:
        if not (SD / "C_letters.tsv").is_file():
            raise SystemExit(f"{STEM}: no source data in {SD}")
        analysis = str(SD)
        log(f"{STEM}: drawing from {SD}")
    else:
        if not os.environ.get("BENCH_ROOT"):
            raise SystemExit(f"{STEM}: set BENCH_ROOT to the benchmark root, or use --from-source-data")
        bench = Path(os.environ["BENCH_ROOT"])
        RUN, GFF, GENOME = bench / RUN_REL, bench / GFF_REL, bench / GENOME_REL
        analysis = str(RUN)
        log(f"{STEM}: reading the nest from {RUN}")
        write_source_data(load(args.verbose))
    els, letters = read_source_data()
    fig, num = build(els, letters)
    num += [("legend: elements in the nest (each other's nest-outer / nest-inner)", len(els),
             "; ".join(sorted({el["tsv"] for el in els.values()}))),
            ("legend: Arabidopsis thaliana chromosome (TAIR10.1 NC_003075.7 = c0004)", 4, str(NAMES_REL)),
            ("legend: family class of all three elements", "LTR/Gypsy/Athila", "seq_id column of the depth tables")]
    res = poafig.write_figure(STEM, fig, LEGEND, None, TMP, analysis, str(HERE / f"{STEM}.py"))
    PNG_OUT.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(res["png"][0], PNG_OUT)
    shutil.move(res["pdf"], HERE / f"{STEM}.pdf")
    pd.DataFrame(num, columns=["item", "value", "source"]).to_csv(HERE / f"{STEM}_numbers.tsv", sep="\t", index=False)
    log(f"{STEM}: wrote {PNG_OUT}, {HERE / f'{STEM}.pdf'}, {HERE / f'{STEM}_numbers.tsv'} (figure "
        f"{tuple(map(float, res['fig_mm'][0]))} mm; legend {res['short_mm']} mm, {res['short_words']} words "
        f"{res['short_par_words']})")


if __name__ == "__main__":
    main()
