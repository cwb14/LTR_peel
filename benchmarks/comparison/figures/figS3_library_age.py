#!/usr/bin/env python3
"""Fig. S3 figS3_library_age: curated libraries annotate old LTR retrotransposons less completely (approach 4, post hoc).

Panels:
  A  one old, library-unsupported, protein-confirmed LTR_peel element of A. thaliana beside the youngest clean copy of
     the same curated family, in element coordinates: LTRs, protein domains, curated hits (family, identity).
  B  every protein-confirmed locus (pooled callers) of rice, zebra finch, maize and gorilla as a row, sorted by LTR
     divergence in four bins: curated LTR coverage along the element scaled to one length, beside the same loci moved
     20-200 kb (chance); mean profiles. The A. thaliana version is in the main figure (fig_benchmark H).
Numbers: report/tables/library_age_bias.tsv; drawings: eval/adj/<ds>/<tool>.calls.tsv, mask/<ds>/curated/<ds>.fa.out
(1-based), runs/arab/dante_ltr/dante.gff3 (1-based), eval/k2p/arab/ltr_peel/annot.tsv, calls/arab/ltr_peel.ids.tsv,
data/arab/genome.fa (+ .fai). Frozen rules imported from scripts/adj_summary.py and scripts/proxies.py (unchanged).
Outputs (this directory): figS3_library_age.{pdf,png}, figS3_library_age_numbers.tsv; caches in tmp/.
Source data: source_data/figS3_library_age/ (A_copies, A_hits, A_domains: the two drawn copies; B_profiles: every
drawn locus row; legend_supported_share). The full mode writes them and draws from the copies read back;
--from-source-data draws from them without BENCH_ROOT.
Usage: figS3_library_age.py [--prep-only] [--redo-prep] [--from-source-data] [-v]
"""
import argparse
import csv
import glob
import gzip
import math
import os
import re
import statistics
import sys
from collections import defaultdict
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

B = Path(os.environ["BENCH_ROOT"]) if os.environ.get("BENCH_ROOT") else None
HERE = Path(__file__).resolve().parent
TMP = HERE / "tmp"
STEM = "figS3_library_age"
SD = HERE / "source_data" / STEM
GZ_OVER = 1 << 20                           # source-data tables larger than this (bytes, uncompressed) are gzipped
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # poafig.py sits beside this script
su = proxies = None                         # frozen rules from BENCH_ROOT/scripts, imported in full mode only


def import_frozen():
    global su, proxies
    sys.path.insert(0, str(B / "scripts"))
    import adj_summary as su  # noqa: E402
    import proxies  # noqa: E402

DS_B = ["rice", "finch", "maize", "gorilla"]          # panel B; A. thaliana is fig_benchmark H
TOOLS = ["ltr_peel", "ltr_retriever", "ltrharvest", "ltr_finder", "look4ltrs", "hirepeat", "hite", "ltrpred",
         "dante_ltr"]                                 # the locus pool (every primary caller with calls)
GEN_NAME = {"rice": "rice", "finch": "zebra finch", "maize": "maize", "gorilla": "gorilla"}
BINS = ["<0.01", "0.01-0.03", "0.03-0.1", ">=0.1"]
NB5, NBI, NB3 = 10, 30, 10                  # profile bins: 5' LTR, internal region, 3' LTR
NBIN = NB5 + NBI + NB3
# the example (selection rule and checks in the code below); arab ids are 1-based closed in the id, 0-based half-open in start
OLD_ID, YOUNG_ID = "c0005:11527310-11538273", "c0001:16626039-16637234"
VERBOSE = False


def log(msg):
    if VERBOSE:
        print(f"{STEM}: {msg}", file=sys.stderr, flush=True)


# ============================================================================ prep: loci, coverage profiles (cached)
def load_coding(ds, tool):
    """coding-unit calls of one tool: the columns the frozen locus rule and the drawings need"""
    out = []
    p = B / f"eval/adj/{ds}/{tool}.calls.tsv"
    with open(p) as f:
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
        """covered bases in [0, x) on c, vectorized over x"""
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
    """bin edges along one element: NB5 bins over [s, a), NBI over [a, b), NB3 over [b, e) (genome orientation)"""
    return np.round(np.concatenate([np.linspace(s, a, NB5 + 1)[:-1], np.linspace(a, b, NBI + 1)[:-1],
                                    np.linspace(b, e, NB3 + 1)])).astype(np.int64)


def profile(cov, c, ed):
    v = cov.cum(c, ed)
    w = np.diff(ed).astype(float)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(w > 0, np.diff(v) / np.maximum(w, 1), np.nan)


def prep(ds):
    """coding-unit loci pooled over callers (adj_summary.age_bias rule), with a per-locus coverage profile of a
    representative member and of its moved span; checked against report/tables/library_age_bias.tsv"""
    pool = [t for t in TOOLS if (B / f"eval/adj/{ds}/{t}.calls.tsv").exists()]
    rows = {t: load_coding(ds, t) for t in pool}
    mem = [(r[0], r[1], r[2], (t, k)) for t in pool for k, r in enumerate(rows[t])]
    log(f"{ds}: {len(mem)} coding calls from {len(pool)} callers; building loci")
    lab, n = su.loci(mem)
    by = defaultdict(list)
    for c, s, e, key in mem:
        by[lab[key]].append(key)
    log(f"{ds}: {n} loci; loading curated LTR intervals")
    cov = Cover(B / f"mask/{ds}/curated/{ds}.fa.out")
    k2p, sup, nul, prof, nprof, rep = [], [], [], [], [], []
    for lid in sorted(by):
        keys = by[lid]
        rs = [rows[t][k] for t, k in keys]
        s_ = sum(r[6] for r in rs) >= len(rs) / 2                       # frozen: >= half of members supported
        ks = [float(r[3]) for r in rs if r[3] != "NA"]
        km = statistics.median(ks) if ks else None
        nv = [r[8] == "1" for r in rs if r[7] != "NA" and r[8] in ("0", "1")]
        nl = (sum(nv) >= len(nv) / 2) if nv else None                   # same majority rule on the moved spans
        # representative: the member with K2P closest to the locus median (pool order breaks ties)
        if km is None:
            r = rs[0]
        else:
            r = min(rs, key=lambda x: abs(float(x[3]) - km) if x[3] != "NA" else math.inf)
        c, s, e, a, b = r[0], r[1], r[2], r[4], r[5]
        if not (s < a < b < e):
            a, b = s + (e - s) // 10, e - (e - s) // 10
        ed = edges(s, e, a, b)
        p = profile(cov, c, ed)
        if r[7] != "NA":
            d = int(r[7])
            q = profile(cov, c, ed + d)
        else:
            q = np.full(NBIN, np.nan)
        k2p.append(np.nan if km is None else km)
        sup.append(s_)
        nul.append(np.nan if nl is None else float(nl))
        prof.append(p)
        nprof.append(q)
        rep.append(f"{c}:{s}-{e}")
    k2p = np.array(k2p)
    sup = np.array(sup)
    nul = np.array(nul)
    np.savez_compressed(TMP / f"prep_{ds}.npz", k2p=k2p, sup=sup, null=nul, prof=np.array(prof, np.float32),
                        nprof=np.array(nprof, np.float32), rep=np.array(rep))
    # check against the frozen table
    tab = {r["age_bin"]: r for r in csv.DictReader(open(B / "report/tables/library_age_bias.tsv"), delimiter="\t")
           if r["ds"] == ds}
    for i, (lo, hi) in enumerate([(0, .01), (.01, .03), (.03, .1), (.1, math.inf)]):
        m = (k2p >= lo) & (k2p < hi)
        t = tab[BINS[i]]
        if int(t["n_loci"]) != int(m.sum()) or abs(float(t["supported_share"]) - sup[m].mean()) > 5e-4:
            raise SystemExit(f"{STEM}: {ds} {BINS[i]}: loci {m.sum()} / share {sup[m].mean():.4f} do not reproduce "
                             f"library_age_bias.tsv ({t['n_loci']} / {t['supported_share']})")
    log(f"{ds}: prep done, reproduces library_age_bias.tsv")


def load_prep(ds):
    z = np.load(TMP / f"prep_{ds}.npz")
    return {k: z[k] for k in z.files}


# ============================================================================ example (panel A, B)
def rm_hits(ds, c, s, e):
    """curated RepeatMasker hits overlapping [s, e): 1-based closed in the file -> 0-based half-open"""
    out = []
    with open(B / f"mask/{ds}/curated/{ds}.fa.out") as f:
        for line in f:
            t = line.split()
            if len(t) > 10 and t[0].isdigit() and t[4] == c:
                a, b = int(t[5]) - 1, int(t[6])
                if b > s and a < e:
                    out.append(dict(s=a, e=b, strand=t[8], name=t[9], cls=t[10], div=float(t[1]),
                                    dl=float(t[2]), ins=float(t[3])))
    return out


def domains(ds, c, s, e):
    """REXdb protein domains (Class_I|LTR) inside [s, e): GFF3 1-based closed -> 0-based half-open"""
    out = []
    with open(B / f"runs/{ds}/dante_ltr/dante.gff3") as f:
        for line in f:
            if line.startswith("#"):
                continue
            t = line.rstrip("\n").split("\t")
            if len(t) < 9 or t[0] != c or t[2] != "protein_domain":
                continue
            a, b = int(t[3]) - 1, int(t[4])
            if a < s or b > e:
                continue
            fc = re.search(r"Final_Classification=([^;]+)", t[8]).group(1)
            if fc.startswith("Class_I|LTR"):
                out.append(dict(s=a, e=b, strand=t[6], name=re.search(r"Name=([^;]+)", t[8]).group(1), cls=fc))
    return out


def example(cid):
    """one arab LTR_peel call: record, Kmer2LTR LTRs and CIGAR, hits, domains, end sequences; structural checks"""
    rec = next(r for r in csv.DictReader(open(B / "eval/adj/arab/ltr_peel.calls.tsv"), delimiter="\t")
               if r["id"] == cid)
    call = next(r for r in csv.DictReader(open(B / "calls/arab/ltr_peel.ids.tsv"), delimiter="\t") if r["id"] == cid)
    k = next(r for r in csv.DictReader(open(B / "eval/k2p/arab/ltr_peel/annot.tsv"), delimiter="\t")
             if r["seq_id"] == cid)
    c, s, e = rec["chrom"], int(rec["start"]), int(rec["end"])
    idx = proxies.fai(str(B / "data/arab/genome.fa"))
    with open(B / "data/arab/genome.fa", "rb") as fh:
        left = proxies.fetch(fh, idx, c, s - 14, s + 10)
        right = proxies.fetch(fh, idx, c, e - 10, e + 14)
    doms = domains("arab", c, s, e)
    hits = [h for h in rm_hits("arab", c, s, e)]
    strand = call["strand"] if call["strand"] in "+-" else doms[0]["strand"]
    ex = dict(id=cid, chrom=c, s=s, e=e, strand=strand, k2p=float(rec["k2p"]), supported=rec["supported"] == "1",
              coding=rec["coding"] == "1", tsd=call["tsd"], left=left, right=right, doms=doms, hits=hits,
              l5=(int(k["ltr5_start"]) - 1, int(k["ltr5_end"])), l3=(int(k["ltr3_start"]) - 1, int(k["ltr3_end"])),
              cigar=k["cigar"], n_x=int(k["n_ts"]) + int(k["n_tv"]), aln_len=int(k["aln_len"]))
    # coverage of the span by curated LTR hits (the support rule)
    cov = Cover(B / "mask/arab/curated/arab.fa.out")
    ex["cov"] = float((cov.cum(c, [e])[0] - cov.cum(c, [s])[0]) / (e - s))
    # structural checks: exact TSD both sides, TG..CA, Kmer2LTR LTRs at the call ends, one strand and one
    # superfamily, canonical domain order
    tsd = call["tsd"]
    chk = dict(tsd=len(tsd) in (4, 5, 6) and left[14 - len(tsd):14] == tsd and right[10:10 + len(tsd)] == tsd,
               tgca=left[14:16] == "TG" and right[8:10] == "CA",
               ltr_bounds=ex["l5"][0] == 0 and ex["l3"][1] == e - s,
               one_strand=len({d["strand"] for d in doms}) == 1,
               one_superfamily=len({d["cls"].split("|")[2] for d in doms}) == 1)
    names = [d["name"] for d in sorted(doms, key=lambda d: d["s"])]
    if strand == "-":
        names = names[::-1]
    ref = ["GAG", "PROT", "RT", "RH", "aRH", "INT", "CHD"]
    chk["domain_order"] = all(n in ref for n in names) and [ref.index(n) for n in names] == \
        sorted(ref.index(n) for n in names)
    chk["core_domains"] = {"GAG", "PROT", "RT", "RH", "INT"} <= set(names)
    if not all(chk.values()):
        raise SystemExit(f"{STEM}: example {cid} fails a structural check: {chk}")
    ex["checks"] = chk
    ex["names"] = names
    return ex


def family_name(n):
    return re.sub(r"(_I|_LTR|_INT|I|LTR)$", "", n)


# ============================================================================ source data
BIN_COLS = ([f"ltr5_{i + 1:02d}" for i in range(NB5)] + [f"int_{i + 1:02d}" for i in range(NBI)]
            + [f"ltr3_{i + 1:02d}" for i in range(NB3)])
COPY_COLS = ["copy", "id", "chrom", "start", "end", "strand", "k2p", "tsd", "ltr5_start", "ltr5_end", "ltr3_start",
             "ltr3_end", "coverage"]
LEGEND_BINS = ("<0.01", ">=0.1")


def write_table(df, name):
    """source_data/<stem>/<name>.tsv, or .tsv.gz (gzip mtime 0, so reruns are byte-identical) above GZ_OVER"""
    SD.mkdir(parents=True, exist_ok=True)
    text = df.to_csv(sep="\t", index=False, na_rep="NA", lineterminator="\n").encode()
    plain, gz = SD / f"{name}.tsv", SD / f"{name}.tsv.gz"
    if len(text) > GZ_OVER:
        with open(gz, "wb") as fh, gzip.GzipFile(filename="", mode="wb", fileobj=fh, mtime=0) as z:
            z.write(text)
        plain.unlink(missing_ok=True)
    else:
        plain.write_bytes(text)
        gz.unlink(missing_ok=True)
    log(f"wrote {name} ({len(df):,} rows)")


def read_table(name, **kw):
    p = SD / f"{name}.tsv"
    if not p.exists():
        p = SD / f"{name}.tsv.gz"
    if not p.exists():
        raise SystemExit(f"{STEM}: source-data table {name}.tsv(.gz) not found in {SD}")
    return pd.read_csv(p, sep="\t", **kw)


def write_source_data():
    """the drawn and printed values: A, the two example copies (record, curated LTR hits, protein domains; genome
    coordinates 0-based half-open, LTR bounds relative to the copy start); B, every drawn locus row (K2P measured)
    in drawing order with its binned coverage profile, observed and moved; the legend's supported shares"""
    old, young = example(OLD_ID), example(YOUNG_ID)
    copies, hits, doms = [], [], []
    for role, ex in (("supported", young), ("unsupported", old)):
        copies.append(dict(zip(COPY_COLS, (role, ex["id"], ex["chrom"], ex["s"], ex["e"], ex["strand"], ex["k2p"],
                                           ex["tsd"], *ex["l5"], *ex["l3"], ex["cov"]))))
        hits += [dict(copy=role, start=h["s"], end=h["e"], name=h["name"], perc_div=h["div"]) for h in ex["hits"]
                 if proxies.LTR.search(h["cls"])]
        doms += [dict(copy=role, start=d["s"], end=d["e"], name=d["name"]) for d in ex["doms"]]
    write_table(pd.DataFrame(copies, columns=COPY_COLS), "A_copies")
    write_table(pd.DataFrame(hits, columns=["copy", "start", "end", "name", "perc_div"]), "A_hits")
    write_table(pd.DataFrame(doms, columns=["copy", "start", "end", "name"]), "A_domains")
    parts = []
    for ds in DS_B:
        d = load_prep(ds)
        ok = ~np.isnan(d["k2p"])
        order = np.lexsort((d["k2p"][ok],))
        for track, M in (("loci", d["prof"]), ("moved", d["nprof"])):
            t = pd.DataFrame(M[ok][order], columns=BIN_COLS)
            t.insert(0, "ds", ds)
            t.insert(1, "row", np.arange(1, int(ok.sum()) + 1))
            t.insert(2, "k2p", d["k2p"][ok][order])
            t.insert(3, "track", track)
            parts.append(t)
    write_table(pd.concat(parts, ignore_index=True), "B_profiles")
    age = [r for r in csv.DictReader(open(B / "report/tables/library_age_bias.tsv"), delimiter="\t")
           if r["ds"] in DS_B and r["age_bin"] in LEGEND_BINS]
    write_table(pd.DataFrame([dict(ds=r["ds"], age_bin=r["age_bin"], supported_share=float(r["supported_share"]))
                              for r in age]), "legend_supported_share")


def read_source_data():
    """the two example copies (as example() returns them, fields drawn only), the per-genome locus profiles (as
    load_prep() returns them, drawn rows only) and the legend's supported shares"""
    cp = read_table("A_copies", keep_default_na=False, dtype={"tsd": str})
    hits = read_table("A_hits")
    doms = read_table("A_domains")
    ex = {}
    for r in cp.to_dict("records"):
        role = r["copy"]
        ex[role] = dict(id=r["id"], chrom=r["chrom"], s=int(r["start"]), e=int(r["end"]), strand=r["strand"],
                        k2p=float(r["k2p"]), tsd=r["tsd"], l5=(int(r["ltr5_start"]), int(r["ltr5_end"])),
                        l3=(int(r["ltr3_start"]), int(r["ltr3_end"])), cov=float(r["coverage"]),
                        hits=[dict(s=int(h["start"]), e=int(h["end"]), name=h["name"], div=float(h["perc_div"]))
                              for h in hits[hits["copy"] == role].to_dict("records")],
                        doms=[dict(s=int(d["start"]), e=int(d["end"]), name=d["name"])
                              for d in doms[doms["copy"] == role].to_dict("records")])
    bp = read_table("B_profiles", dtype={c: np.float32 for c in BIN_COLS})
    prof = {}
    for ds in DS_B:
        lo, mv = bp[(bp["ds"] == ds) & (bp["track"] == "loci")], bp[(bp["ds"] == ds) & (bp["track"] == "moved")]
        if len(lo) != len(mv) or not np.array_equal(lo["k2p"].to_numpy(), mv["k2p"].to_numpy()):
            raise SystemExit(f"{STEM}: B_profiles: {ds} loci and moved rows do not pair up")
        prof[ds] = dict(k2p=lo["k2p"].to_numpy(np.float64), prof=lo[BIN_COLS].to_numpy(np.float32),
                        nprof=mv[BIN_COLS].to_numpy(np.float32))
    sh = read_table("legend_supported_share")
    share = {(r["ds"], r["age_bin"]): float(r["supported_share"]) for r in sh.to_dict("records")}
    return dict(old=ex["unsupported"], young=ex["supported"], prof=prof, share=share)


def check_read_back(sd):
    """full mode: the profiles read back equal the caches bit for bit (drawn rows, drawing order)"""
    for ds in DS_B:
        d = load_prep(ds)
        ok = ~np.isnan(d["k2p"])
        order = np.lexsort((d["k2p"][ok],))
        r = sd["prof"][ds]
        if not (np.array_equal(r["k2p"], d["k2p"][ok][order])
                and np.array_equal(r["prof"], d["prof"][ok][order], equal_nan=True)
                and np.array_equal(r["nprof"], d["nprof"][ok][order], equal_nan=True)):
            raise SystemExit(f"{STEM}: B_profiles read back for {ds} differ from tmp/prep_{ds}.npz")


# ============================================================================ drawing
def setup_fonts():
    """Liberation Sans registered in-process; poafig._advance_pt patched for matplotlib >= 3.11 (poafig.py itself is not modified)"""
    import poafig
    import matplotlib.font_manager as fm
    for f in sorted(glob.glob("/usr/share/fonts/truetype/liberation2/LiberationSans-*.ttf")
                    + glob.glob("/usr/share/fonts/truetype/liberation2/LiberationMono-*.ttf")):
        fm.fontManager.addfont(f)
    poafig.setup_mpl()

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

    import inspect
    from matplotlib import _text_helpers
    if "kern_mode" not in inspect.signature(_text_helpers.layout).parameters:
        poafig._advance_pt = _advance_pt_mpl311
    return poafig


INK, MUTED, HITC, NULLC = "#1A1A1A", "#6B6B6B", "#4D4D4D", "#8C8C8C"
LTR5C, LTR3C, INTC = "#D55E00", "#EE9F6B", "#F6D9C7"
DOMC = {"GAG": "#D9CCEA", "PROT": "#C2ADDB", "RT": "#A68BC9", "RH": "#8E6FBA", "aRH": "#7A5AA9", "INT": "#5E3F8E",
        "CHD": "#4A2F73"}
LOCIC = "#6A3D9A"           # protein-confirmed loci (B profiles): the domain hue; genomes are named, not colored
FS, FSS = 6.0, 5.0           # text sizes (pt): labels, small labels


def lighten(hexc, f):
    import matplotlib.colors as mc
    r, g, b = mc.to_rgb(hexc)
    return (r + (1 - r) * f, g + (1 - g) * f, b + (1 - b) * f)


def fmt_k(x):
    return f"{x:.3f}"


def revcomp(s):
    return s[::-1].translate(str.maketrans("ACGTN", "TGCAN"))


def to_elem(ex, a, b):
    """genome [a, b) -> element coordinates 5'->3' (0 = element start)"""
    if ex["strand"] == "-":
        return ex["e"] - b, ex["e"] - a
    return a - ex["s"], b - ex["s"]


def draw_copy(ax, ex, y, x0, mm_per_bp, label, numbers):
    """one element in element coordinates: curated hits (top), element (middle), domains (bottom); y = top (mm)"""
    from matplotlib.patches import Rectangle, Polygon
    L = ex["e"] - ex["s"]
    X = lambda p: x0 + p * mm_per_bp                     # noqa: E731
    yh, ye, yd = y + 3.2, y + 7.4, y + 10.6             # hit box top, element centre, domain box top
    # row labels
    ax.text(1.0, ye, label[0], ha="left", va="bottom", fontsize=FS, fontweight="bold")
    ax.text(1.0, ye + 0.5, label[1], ha="left", va="top", fontsize=FS)
    ax.text(x0 - 1.6, yh + 1.0, "curated hits", ha="right", va="center", fontsize=FSS, color=MUTED)
    ax.text(x0 - 1.6, yd + 1.0, "protein domains", ha="right", va="center", fontsize=FSS, color=MUTED)
    # flanks and internal region
    ax.plot([X(-260), X(L + 260)], [ye, ye], color=NULLC, lw=0.6, solid_capstyle="butt", zorder=1)
    l5 = to_elem(ex, ex["s"] + ex["l5"][0], ex["s"] + ex["l5"][1]) if ex["strand"] == "+" else \
        to_elem(ex, ex["s"] + ex["l3"][0], ex["s"] + ex["l3"][1])
    l3 = to_elem(ex, ex["s"] + ex["l3"][0], ex["s"] + ex["l3"][1]) if ex["strand"] == "+" else \
        to_elem(ex, ex["s"] + ex["l5"][0], ex["s"] + ex["l5"][1])
    ax.add_patch(Rectangle((X(l5[1]), ye - 0.8), X(l3[0]) - X(l5[1]), 1.6, fc=INTC, ec="none", zorder=2))
    for (a, b), col in ((l5, LTR5C), (l3, LTR3C)):
        w = X(b) - X(a)
        tip = min(1.2, w * 0.4)
        ax.add_patch(Polygon([(X(a), ye - 1.3), (X(b) - tip, ye - 1.3), (X(b), ye), (X(b) - tip, ye + 1.3),
                              (X(a), ye + 1.3)], closed=True, fc=col, ec="none", zorder=3))
    tl = len(ex["tsd"]) * mm_per_bp
    for xx in (X(0) - max(tl, 0.5), X(L)):
        ax.add_patch(Rectangle((xx, ye - 0.7), max(tl, 0.5), 1.4, fc="white", ec=INK, lw=0.45, zorder=4))
    # curated hits: every hit, labeled with its family part and identity to the consensus (100 - divergence)
    for h in sorted(ex["hits"], key=lambda h: h["s"]):
        a, b = to_elem(ex, max(h["s"], ex["s"] - 200), min(h["e"], ex["e"] + 200))
        ax.add_patch(Rectangle((X(a), yh), X(b) - X(a), 2.0, fc=HITC, ec="white", lw=0.3, zorder=3))
        fam = family_name(h["name"]).replace("_ATH", "")
        part = "LTR" if re.search(r"LTR$", h["name"]) else "int"
        ident = round(100 - h["div"])
        numbers.append((f"A {label[0]} curated hit {h['name']} identity (%)", ident,
                        f"mask/arab/curated/arab.fa.out (100 - perc div. {h['div']}); {h['s'] + 1}-{h['e']}"))
        w = X(b) - X(a)
        txt = next((t for t in (f"{fam} {part}", fam, part) if 0.95 * len(t) + 1.0 < w), "")
        if txt:
            ax.text((X(a) + X(b)) / 2, yh + 1.0, txt, ha="center", va="center", fontsize=FSS, color="white", zorder=4)
        ax.text((X(a) + X(b)) / 2, yh - 0.3, f"{ident}%", ha="center", va="bottom", fontsize=FSS, color=INK)
    # domains in element order
    for d in ex["doms"]:
        a, b = to_elem(ex, d["s"], d["e"])
        ax.add_patch(Rectangle((X(a), yd), X(b) - X(a), 2.0, fc=DOMC.get(d["name"], "#7A5AA9"), ec="none", zorder=3))
        ax.text((X(a) + X(b)) / 2, yd + 2.3, d["name"], ha="center", va="top", fontsize=FSS, color=INK)


def build(sd, analysis_dir):
    poafig = setup_fonts()
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    from matplotlib.colors import LinearSegmentedColormap
    MM = poafig.MM
    numbers = []

    def num(name, value, source):
        numbers.append((name, value, source))

    W, H = 180.0, 117.0
    fig = plt.figure(figsize=(W * MM, H * MM))

    # ------------------------------------------------------------------ A: the example pair
    old, young = sd["old"], sd["young"]
    yA = 3.0
    axA = poafig.mm_axes(fig, 0, yA, W, 34)
    axA.set_xlim(0, W)
    axA.set_ylim(34, 0)
    axA.axis("off")
    x0 = 31.0
    Lmax = max(young["e"] - young["s"], old["e"] - old["s"])
    mm_per_bp = (W - 4 - x0) / (Lmax + 300)
    draw_copy(axA, young, 0.5, x0, mm_per_bp, ("supported copy", f"K2P {fmt_k(young['k2p'])}"), numbers)
    draw_copy(axA, old, 17.0, x0, mm_per_bp, ("unsupported copy", f"K2P {fmt_k(old['k2p'])}"), numbers)
    sx = W - 4 - 1000 * mm_per_bp                     # scale bar: 1 kb
    axA.plot([sx, W - 4], [32.8, 32.8], color=INK, lw=0.8, solid_capstyle="butt")
    axA.text((sx + W - 4) / 2, 32.4, "1 kb", ha="center", va="bottom", fontsize=FSS)
    poafig.panel_letter(fig, 1.0, yA - 1.0, "A")
    for nm, ex in (("supported", young), ("unsupported", old)):
        num(f"A {nm} copy K2P (LTR-LTR, substitutions per site)", ex["k2p"],
            "eval/adj/arab/ltr_peel.calls.tsv k2p (Kmer2LTR)")
        num(f"A {nm} copy curated LTR coverage of the span", round(ex["cov"], 3),
            "computed from mask/arab/curated/arab.fa.out (proxies.ref_ltr rule)")
        num(f"A {nm} copy length (bp)", ex["e"] - ex["s"], "eval/adj/arab/ltr_peel.calls.tsv len")
        num(f"A {nm} copy TSD (open boxes)", ex["tsd"], "calls/arab/ltr_peel.ids.tsv tsd; checked against "
            "data/arab/genome.fa")

    # ------------------------------------------------------------------ B: every locus as a row
    # four genome columns over the width of A's element track (left edges line up at x0), key row beneath
    yB = 41.0
    poafig.panel_letter(fig, 1.0, yB - 0.5, "B")
    XR, CG = W - 4, 3.0
    CW = (XR - x0 - 3 * CG) / 4
    colx = [x0 + i * (CW + CG) for i in range(4)]
    BW = (CW - 1.4) / 2                            # observed / moved block width
    cmap = LinearSegmentedColormap.from_list("cov", ["#FFFFFF", "#262626"])
    yP, hP = yB + 9.2, 8.0                          # profile axes
    yH, hH = yP + hP + 1.4, 48.0                    # heat axes
    for i, ds in enumerate(DS_B):
        first = i == 0
        d = sd["prof"][ds]
        k = d["k2p"]
        ok = ~np.isnan(k)
        bins = np.digitize(k[ok], [0.01, 0.03, 0.1])
        order = np.lexsort((k[ok],))
        P, Q = d["prof"][ok][order], d["nprof"][ok][order]
        bb = bins[order]
        n_tot = int(ok.sum())
        fig.text((colx[i] + CW / 2) / W, 1 - yB / H, f"{GEN_NAME[ds]}, n = {n_tot:,}", ha="center", va="top",
                 fontsize=FS, color=INK)
        num(f"B {ds} protein-confirmed loci drawn (K2P measured)", n_tot,
            "report/tables/library_age_bias.tsv n_loci (sum of the four bins); rows rebuilt from eval/adj")
        for j, (M, title) in enumerate(((P, "loci"), (Q, "moved"))):
            bx = colx[i] + j * (BW + 1.4)
            axg = poafig.mm_axes(fig, bx, yB + 6.0, BW, 2.4)          # element glyph
            axg.set_xlim(0, NBIN)
            axg.set_ylim(-1, 1)
            axg.axis("off")
            axg.add_patch(Rectangle((NB5, -0.35), NBI, 0.7, fc=INTC, ec="none"))
            axg.add_patch(Rectangle((0, -0.9), NB5, 1.8, fc=LTR5C, ec="none"))
            axg.add_patch(Rectangle((NB5 + NBI, -0.9), NB3, 1.8, fc=LTR3C, ec="none"))
            fig.text((bx + BW / 2) / W, 1 - (yB + 5.7) / H, title, ha="center", va="bottom", fontsize=FSS,
                     color=INK if j == 0 else MUTED)
            axp = poafig.mm_axes(fig, bx, yP, BW, hP)                 # mean profiles
            xs = np.arange(NBIN) + 0.5
            if j == 0:
                for b_, shade in ((0, 0.0), (3, 0.55)):
                    prof_ = np.nanmean(M[bb == b_], axis=0)
                    axp.plot(xs, prof_, color=lighten(LOCIC, shade), lw=0.9)
                    num(f"B {ds} mean curated coverage along loci, K2P class {BINS[b_]} (mean over bins)",
                        round(float(np.nanmean(prof_)), 3), "computed from eval/adj + mask/<ds>/curated (tmp cache)")
            else:
                prof_ = np.nanmean(M, axis=0)
                axp.plot(xs, prof_, color=NULLC, lw=0.9)
                num(f"B {ds} mean curated coverage along moved loci (mean over bins)",
                    round(float(np.nanmean(prof_)), 3), "computed from eval/adj null_shift + mask/<ds>/curated")
            axp.set_xlim(0, NBIN)
            axp.set_ylim(0, 1.05)
            for s_ in ("top", "right"):
                axp.spines[s_].set_visible(False)
            axp.set_xticks([])
            axp.set_yticks([0, 1])
            if first and j == 0:
                axp.set_yticklabels(["0", "1"])
                axp.set_ylabel("mean\ncoverage", fontsize=FSS, labelpad=1)
            else:
                axp.set_yticklabels([])
            axh = poafig.mm_axes(fig, bx, yH, BW, hH)                 # rows: equal height per divergence class
            for b_ in range(4):
                m = bb == b_
                img = np.nan_to_num(M[m], nan=0.0)
                if img.shape[0] == 0:
                    continue
                axh.imshow(img, cmap=cmap, vmin=0, vmax=1, aspect="auto", interpolation="antialiased",
                           extent=(0, NBIN, b_ + 1, b_), origin="upper")
            for b_ in (1, 2, 3):
                axh.axhline(b_, color="white", lw=0.8)
            axh.set_xlim(0, NBIN)
            axh.set_ylim(4, 0)
            if first:
                axh.set_xticks([NB5 / 2, NB5 + NBI / 2, NB5 + NBI + NB3 / 2])
                axh.set_xticklabels(["LTR", "internal", "LTR"], fontsize=FSS)
                axh.tick_params(axis="x", length=0, pad=1.2)
            else:
                axh.set_xticks([])
            for s_ in axh.spines.values():
                s_.set_linewidth(0.4)
            axh.set_yticks([1, 2, 3])
            if first and j == 0:
                axh.set_yticklabels(["0.01", "0.03", "0.1"])
                axh.set_ylabel("LTR divergence (K2P)", fontsize=FS, labelpad=1)
            else:
                axh.set_yticklabels([])
        for b_ in range(4):
            num(f"B {ds} loci in K2P class {BINS[b_]}", int((bb == b_).sum()),
                "report/tables/library_age_bias.tsv n_loci")
    num("B A. thaliana panel", "in fig_benchmark G (main figure, Fig. 1G), not redrawn here",
        "scripts/adj_summary.py frozen rules")
    # key row beneath B: shading scale, then the three profile lines
    yK = yH + hH + 6.2
    axk = poafig.key_axes(fig, x0, yK - 1.6, XR - x0, 3.4)
    gw = 12.0
    for t_ in range(40):
        axk.add_patch(Rectangle((t_ * gw / 40, 0.6), gw / 40 + 0.02, 1.6, fc=cmap((t_ + 0.5) / 40), ec="none"))
    axk.add_patch(Rectangle((0, 0.6), gw, 1.6, fc="none", ec=INK, lw=0.3))
    axk.text(-0.6, 1.4, "0", ha="right", va="center", fontsize=FSS)
    axk.text(gw + 0.6, 1.4, "100% curated annotation", ha="left", va="center", fontsize=FSS)
    items = [(poafig.g_line(LOCIC, lw=1.0), "K2P < 0.01"),
             (poafig.g_line(lighten(LOCIC, 0.55), lw=1.0), "K2P ≥ 0.1"),
             (poafig.g_line(NULLC, lw=1.0), "moved 20–200 kb")]
    poafig.hkey(axk, items, 50.0, 1.4, XR - x0 - 50, fontsize=FSS)
    log(f"figure bottom at {yK + 2:.0f} mm of {H:.0f}")
    short = legend(sd["share"], numbers, old, young)
    res = poafig.write_figure(STEM, fig, short, None, HERE, analysis_dir=analysis_dir,
                              script=str(HERE / f"{STEM}.py"))
    with open(HERE / f"{STEM}_numbers.tsv", "w") as f:
        f.write("item\tvalue\tsource\n")
        for n_, v, s_ in numbers:
            f.write(f"{n_}\t{v}\t{s_}\n")
    print(f"{STEM}: wrote {res['pdf']} ({res['n_pages']} pages; legend {res['short_words']} words, "
          f"{res['short_mm']} mm), {res['png'][0]}, {HERE / (STEM + '_numbers.tsv')}", file=sys.stderr)


def legend(share, numbers, o, y):
    """legend page text (poafig markup); every number in it is recorded in numbers"""
    def n(name, value, source):
        numbers.append((name, value, source))
        return value
    sh = {(d, b): share[(d, b)] for d in DS_B for b in LEGEND_BINS}
    lo_y, hi_y = min(sh[(d, "<0.01")] for d in DS_B), max(sh[(d, "<0.01")] for d in DS_B)
    lo_o, hi_o = min(sh[(d, ">=0.1")] for d in DS_B), max(sh[(d, ">=0.1")] for d in DS_B)
    n("legend supported share, youngest class, range over the four genomes of B", f"{lo_y:.2f}-{hi_y:.2f}",
      "report/tables/library_age_bias.tsv supported_share <0.01 (rice, finch, maize, gorilla)")
    n("legend supported share, oldest class, range over the four genomes of B", f"{lo_o:.2f}-{hi_o:.2f}",
      "report/tables/library_age_bias.tsv supported_share >=0.1 (rice, finch, maize, gorilla)")
    n("legend rule: moved-locus shift (kb)", "20-200", "memo.md 'Approach 4' local null")
    n("legend rule: library-supported = annotated share of the span (%)", 80, "memo.md 'Approach 4' (proxies.py)")
    n("legend genomes (B here + A. thaliana in Fig. 1G)", 5, "report/tables/library_age_bias.tsv")
    hid = [round(100 - h["div"]) for h in o["hits"]]
    yid = [round(100 - h["div"]) for h in y["hits"]]
    n("legend unsupported copy hit identity range (%)", f"{min(hid)}-{max(hid)}", "mask/arab/curated/arab.fa.out")
    n("legend supported copy hit identity range (%)", f"{min(yid)}-{max(yid)}", "mask/arab/curated/arab.fa.out")
    n("legend unsupported copy annotated share (%)", round(100 * o["cov"]), "computed (proxies.ref_ltr rule)")
    return f"""**Fig. S3. Curated repeat libraries annotate old LTR retrotransposons less completely.** In five genomes, \
LTR retrotransposons that carry their own protein domains are less completely annotated by the curated library the more \
their two LTRs have diverged, so masking benchmarks count parts of real elements as false positives. The analysis was \
added after the benchmark design was fixed.

**A**, An *A. thaliana* element called by LTR_peel (unsupported copy) carries the protein domains of a coding \
element in their canonical order, yet the library annotates only {round(100 * o['cov'])}% of it, through scattered hits \
of {min(hid)}–{max(hid)}% identity and none on its LTRs. A clean copy of the same curated family (supported copy) is \
annotated end to end at {min(yid)}–{max(yid)}% identity. Both are drawn 5′ to 3′, with their target-site duplications \
as small open boxes.

**B**, Every protein-confirmed locus found by any caller in rice, zebra finch, maize and gorilla, one row each, in four \
classes of LTR divergence (substitutions per site), young at the top. Shading is the curated annotation along the \
element scaled to one length, beside the same loci moved 20–200 kb along the chromosome (chance), with the means as lines above. The share of loci the library supports (at least 80% annotated) falls from {lo_y:.2f}–{hi_y:.2f} in the \
youngest class to {lo_o:.2f}–{hi_o:.2f} in the oldest. The *A. thaliana* version of B is Fig. 1G."""


def main(argv=None):
    global VERBOSE
    ap = argparse.ArgumentParser(description="Fig. S3 figS3_library_age: figure page, legend page, PNG, numbers.")
    ap.add_argument("--prep-only", action="store_true", help="build the per-genome locus caches in tmp/ and stop")
    ap.add_argument("--redo-prep", action="store_true", help="rebuild the caches even if present")
    ap.add_argument("--from-source-data", action="store_true",
                    help=f"redraw from source_data/{STEM}/ only; needs no BENCH_ROOT")
    ap.add_argument("-v", "--verbose", action="store_true", help="progress messages")
    a = ap.parse_args(argv)
    VERBOSE = a.verbose
    if a.from_source_data:
        build(read_source_data(), str(SD))
        return
    if B is None:
        raise SystemExit(f"{STEM}: BENCH_ROOT is not set (or use --from-source-data)")
    import_frozen()
    TMP.mkdir(exist_ok=True)
    for ds in DS_B:
        if a.redo_prep or not (TMP / f"prep_{ds}.npz").exists():
            prep(ds)
    if a.prep_only:
        return
    write_source_data()
    sd = read_source_data()
    check_read_back(sd)
    build(sd, str(B))


if __name__ == "__main__":
    main()
