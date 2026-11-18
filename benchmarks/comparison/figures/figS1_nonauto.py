#!/usr/bin/env python3
"""Fig. S1 figS1_nonauto: size floors and protein rules remove protein-less LTR retrotransposons.

Panels. Classes: protein-less = no LTR-RT protein domain; coding = at least
two kinds of LTR-RT domain from one superfamily; partly coding = the rest (table key 'partial').
  A  both ends of a typical protein-less rice element (a Dasheng copy) at bp: target-site duplications and
     the TG..CA ends; the element named in-panel
  B  every protein-less and coding simulated element of S1 and S2 (3 seeds
     each; coding class only, partly coding not drawn) as a column sorted by length, protein-less block beside the
     coding block, found (tool color) or missed; 2 kb and 4 kb marked
  C  recall vs element length per tool, coding (solid) and
     protein-less (dashed), family-bootstrap 95% bands; line-style key in the header
  D  (new) scorecard: precision, recall, F1 for coding / partly coding / protein-less elements, S2, mean of 3 genomes
     (n summed), D_TABLE (constant at the top; now experiments/figpub/sim_class_prf_joint_summary.tsv) system
     'protein' (one basis with Fig. 1B, whose coding and protein-less columns these equal); viridis on a fixed 0-1
     scale, column best bold; the legend's basis clause checked by check_basis()
  E  library-supported share of each tool's protein-less calls with block-bootstrap CI, chance (same spans
     moved 20-200 kb) and the same tool's coding calls on the same axis; n of protein-less calls; failed runs dark
The aging series is in the main figure (fig_benchmark C-D).
Numbers: experiments/nonauto/*.tsv (sim_elements, sim_gap, div_d50, real_calls, real_loci, real_locus_table),
  D_TABLE (D; its per-genome rows and the call side of sim_nonauto_prf_ext.tsv for checks), sims/truth/S2_*.truth.tsv
  + calls/S2_*/<tool>.tsv + experiments/nonauto/lib/lib_class.tsv (check_basis), report/tables/resources.tsv (failed
  runs). Drawings:
  data/rice/genome.fa (+ .fai), calls/rice/<tool>.tsv (0-based half-open, original chromosome names ->
  data/rice/names.tsv), calls/rice/ltrquest.ids.tsv (TSD), runs/rice/dante_ltr/dante.gff3 (1-based closed),
  mask/rice/curated/rice.fa.out (1-based closed) for the example checks.
Outputs (this directory): figS1_nonauto.{pdf,png}, figS1_nonauto_numbers.tsv.
Usage: figS1_nonauto.py [-v]
"""
import argparse
import csv
import glob
import os
import re
import sys
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

B = Path(os.environ["BENCH_ROOT"])
NA = B / "experiments" / "nonauto"
FP = B / "experiments" / "figpub"
# D's table (the same table as Fig. 1B): the one place to switch it. Its per-genome rows are read from the same name
# without '_summary'. check_basis() asserts the legend's D basis clause ("every simulated element, including copies of
# library families without annotated LTRs") against it: a variant that leaves those families out (e.g.
# sim_class_prf_joint_hasltr_summary.tsv) stops the build there until that clause is reworded.
D_TABLE = FP / "sim_class_prf_joint_summary.tsv"
D_SRC = f"{D_TABLE.relative_to(B)} (system protein)"
HERE = Path(__file__).resolve().parent
STEM = "figS1_nonauto"
sys.path.append(str(B / "scripts"))          # appended: scripts/queue.py would shadow the stdlib queue
import proxies  # noqa: E402
from score_elements import match, read_calls, read_truth  # noqa: E402  (check_basis: D's matching reproduced)

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
from matplotlib.patches import Rectangle  # noqa: E402
from matplotlib.colors import to_rgb  # noqa: E402

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
CLS2 = ["proteinless", "coding"]
CLS3 = ["coding", "partial", "proteinless"]               # D column groups (table keys)
# names on the figure, in the legend and in numbers.tsv: 'partial' is shown as "partly coding"
CLS3_NAME = {"coding": "coding", "partial": "partly coding", "proteinless": "protein-less"}
METRICS = [("precision", "prec."), ("recall", "rec."), ("f1", "F1")]   # order of Fig. 1B
# missing tool x genome results, set-wide rule: two styles only. Tool side ("tool failed", dark
# gray, white text) = report/tables/resources.tsv status timeout / mem_guard / failed; evaluation side ("not
# evaluated", light hatch) = our evaluation step did not run (no such cell in this figure). Built in load_missing().
FAILED_STATUS = {"timeout", "mem_guard", "failed"}
FAILC = "#5A5A5A"
EXPECT_FAILED = {("finch", "ltr_retriever"), ("maize", "ltr_retriever"), ("gorilla", "ltr_retriever"),
                 ("gorilla", "ltrharvest"), ("gorilla", "ltrpred"), ("gorilla", "look4ltrs"),
                 ("finch", "dante_ltr"), ("gorilla", "dante_ltr")}
MISSING = {}                                  # (ds, tool) -> label; filled by load_missing()
# the example (selection rule and checks in the code below); coordinates 0-based half-open on data/rice/genome.fa ids
EX = [dict(key="pl", cls="proteinless", locus="22", chrom="c0001", s=2007764, e=2014480, title="protein-less",
           name="Dasheng", chrom_name="Chr1")]
FLANK = 450                                   # bp of flank drawn on each side in A
MM = poafig.MM
INK, MUTED = poafig.INK, poafig.MUTED
LTRC, INTC, HITC = "#2E2E2E", "#D4D4D4", "#5E5E5E"
MISS_SIM = "#E4E4E4"                          # B: element missed
UNDEFC = "#EEEEEE"                            # D: precision / F1 undefined (no calls of that class)
NULLC = "#8C8C8C"
FS = FSS = FST = 6.0                           # every label, small label and tick label: 6 pt (brief: >= 6 pt)
# score map of the whole figure set (shared by the figure set): viridis on a fixed 0-1 scale, CMAP(v) with v in [0, 1]
# (cividis_r made 0.97 and 0.67 hard to tell apart); text color by cell luminance (text_color)
CMAP = mpl.colormaps["viridis"]
RNG_SEED = 1                                   # family bootstrap draws (as nonauto_stats.py)
NB = 1000
FIG_H = 144.0                                  # figure height (mm)
VERBOSE = False
NUM = []


def log(msg):
    if VERBOSE:
        print(f"[{STEM}] {msg}", file=sys.stderr, flush=True)


def num(panel, ds, tool, quantity, value, source):
    NUM.append(dict(panel=panel, dataset=ds, tool=tool, quantity=quantity, value=value, source=source))
    return value


def tint(c, a):
    """color c mixed with white; a = share of the color"""
    r, g, b = to_rgb(c)
    return (1 - a + a * r, 1 - a + a * g, 1 - a + a * b)


def tsv(p):
    return pd.read_csv(p, sep="\t")


def text_color(v):
    """white on dark cells, ink on light ones (relative luminance of the cell color; as fig3_real_genomes.py)"""
    r, g, b = (c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4 for c in to_rgb(CMAP(v)))
    return "white" if 0.2126 * r + 0.7152 * g + 0.0722 * b < 0.18 else INK


# ============================================================================ data: simulations (B, C)
D_EDGES = {"proteinless": np.array([350, 1000, 1414, 2000, 2828, 4000, 5657, 8000, 11314, 16000]),
           "coding": np.array([1000, 2000, 2828, 4000, 5657, 8000, 11314, 20000])}


def load_sim():
    e = pd.read_csv(NA / "sim_elements.tsv", sep="\t")
    # B's counting, every class (check_basis): elements per S2 genome and their families (one row per element per tool)
    e2 = e[(e.tool == "ltrquest") & e.ds.str.startswith("S2")]
    b_basis = {"n": e2.groupby("ds").size().to_dict(), "families": set(e2.family)}
    e = e[e.tool.isin(TOOLS) & e.cls.isin(CLS2)].copy()
    e["set"] = e.ds.str[:2]
    gap = tsv(NA / "sim_gap.tsv")
    # cross-check: per set recall from the element table = sim_gap.tsv
    for (st, tool, cls), r in e.groupby(["set", "tool", "cls"]).hit.mean().items():
        row = gap[(gap.set == st) & (gap.tool == tool)].iloc[0]
        ref = row["recall_coding" if cls == "coding" else "recall_proteinless"]
        if abs(r - ref) > 6e-5:
            raise SystemExit(f"{STEM}: sim recall {st} {tool} {cls} {r:.4f} != sim_gap {ref}")
    # size floors in the simulations (legend D): elements found below 2 kb (HiTE, HiRepEat) and 4 kb (LTRpred, DANTE)
    floors = {}
    for tool, kb in (("hirepeat", 2000), ("hite", 2000), ("ltrpred", 4000), ("dante_ltr", 4000)):
        for cls in CLS2:
            s = e[(e.tool == tool) & (e.cls == cls) & (e.len < kb)]
            floors[(tool, cls)] = (int(s.hit.sum()), len(s), int(s[s.hit == 1].len.min()) if s.hit.any() else None)
            num("legend/C", "S1+S2", tool,
                f"{CLS3_NAME[cls]} elements < {kb // 1000} kb found / present (shortest found, bp)",
                f"{floors[(tool, cls)][0]} / {floors[(tool, cls)][1]} ({floors[(tool, cls)][2]})",
                "sim_elements.tsv (hit, len)")
    out = {"floors": floors, "b_basis": b_basis}
    rng = np.random.default_rng(RNG_SEED)
    for cls in CLS2:
        s = e[e.cls == cls]
        p = s.pivot_table(index=["ds", "element", "family", "len"], columns="tool", values="hit",
                          aggfunc="first").reset_index()
        if p[TOOLS].isna().any().any():
            raise SystemExit(f"{STEM}: an element lacks a tool row ({cls})")
        p = p.sort_values(["len", "ds", "element"], kind="mergesort").reset_index(drop=True)
        hits = p[TOOLS].to_numpy().T.astype(np.int8)               # tools x elements
        L = p["len"].to_numpy()
        # D: recall per length bin with a family bootstrap (families = the independent unit)
        ed = D_EDGES[cls]
        if L.min() < ed[0] or L.max() >= ed[-1]:
            raise SystemExit(f"{STEM}: {cls} lengths outside the D bins")
        b = np.searchsorted(ed, L, side="right") - 1
        fams = sorted(p.family.unique())
        fi = p.family.map({f: i for i, f in enumerate(fams)}).to_numpy()
        nF, nb = len(fams), len(ed) - 1
        Cnt = np.zeros((nF, nb))
        np.add.at(Cnt, (fi, b), 1)
        H = np.zeros((len(TOOLS), nF, nb))
        for t in range(len(TOOLS)):
            np.add.at(H[t], (fi, b), hits[t])
        draws = rng.integers(0, nF, size=(NB, nF))
        M = np.zeros((NB, nF))
        for i in range(NB):
            M[i] = np.bincount(draws[i], minlength=nF)
        den = M @ Cnt
        rec = H.sum(1) / Cnt.sum(0)                                 # tools x bins
        lo, hi = np.zeros_like(rec), np.zeros_like(rec)
        for t in range(len(TOOLS)):
            bs = (M @ H[t]) / np.where(den > 0, den, np.nan)
            lo[t], hi[t] = np.nanpercentile(bs, 2.5, axis=0), np.nanpercentile(bs, 97.5, axis=0)
        nfam_bin = (Cnt > 0).sum(0)
        log(f"sim {cls}: {len(p)} elements, {nF} families; per D bin elements {Cnt.sum(0).astype(int).tolist()}, "
            f"families {nfam_bin.tolist()}")
        out[cls] = dict(p=p, hits=hits, L=L, edges=ed, rec=rec, lo=lo, hi=hi, n_bin=Cnt.sum(0), nfam_bin=nfam_bin,
                        n_fam=nF, recall=hits.mean(1))
    out["gap"] = gap
    return out


# ============================================================================ data: scorecard (D)
def d_table_genomes():
    """per-genome rows of D_TABLE: the same file name without '_summary'"""
    p = D_TABLE.with_name(D_TABLE.name.replace("_summary", ""))
    if p == D_TABLE or not p.exists():
        raise SystemExit(f"{STEM}: no per-genome table beside {D_TABLE} (expected {p.name})")
    return p


def load_prf():
    """D: D_TABLE (now experiments/figpub/sim_class_prf_joint_summary.tsv), system 'protein' (the basis of Fig. 1B):
    truth = every truth.tsv element (ext), one-to-one RO >= 0.9; true element class = its library entry's class, call
    class = the REXdb domains inside the call; precision, recall, F1 = mean of the 3 genomes, n_truth / n_calls summed.
    Checked here against the per-genome rows (d_table_genomes()) and, for the call side, against
    sim_nonauto_prf_ext.tsv (same calls, same matching). Returns {(set, tool, cls): dict(precision, recall, f1,
    n_truth, tp_truth, n_calls, tp_calls, n_genomes)}; precision / F1 NaN where a tool makes no call of a class."""
    s = tsv(D_TABLE)
    s = s[(s.system == "protein") & s.tool.isin(TOOLS)]
    g = tsv(d_table_genomes())
    g = g[(g.system == "protein") & g.tool.isin(TOOLS)]
    # basis guard: the legend's D clause ("every simulated element") and the checks below (call side = the every-element
    # matching of sim_nonauto_prf_ext.tsv; check_basis) hold only for a table that counts every truth.tsv row
    for st in ("S1", "S2"):
        for seed in (1, 2, 3):
            n_rows = len(read_truth(str(B / f"sims/truth/{st}_{seed}.truth.tsv")))
            n_tab = int(g[(g["set"] == st) & (g.seed == seed) & (g.tool == "ltrquest")].n_truth.sum())
            if n_tab != n_rows:
                raise SystemExit(f"{STEM}: {D_TABLE.name} counts {n_tab} of the {n_rows} simulated elements of "
                                 f"{st}_{seed}. The legend's D clause 'Like Fig. 1B, D counts every simulated "
                                 "element, including copies of library families without annotated LTRs ...' and the "
                                 "basis checks "
                                 "(load_prf call side, check_basis) assume every element: reword and adapt them first")
    ext = tsv(FP / "sim_nonauto_prf_ext.tsv")
    ext = ext[ext.tool.isin(TOOLS) & (ext.seed.astype(str) != "pooled")]
    if set(ext.precision_basis) != {"ext"}:
        raise SystemExit(f"{STEM}: sim_nonauto_prf_ext.tsv has a basis other than ext")
    ext = ext.assign(seed=ext.seed.astype(int)).set_index(["set", "seed", "tool", "cls"])
    if set(s.cls) != set(CLS3) or set(g.cls) != set(CLS3) or set(s.role) != {"primary"}:
        raise SystemExit(f"{STEM}: unexpected classes or roles in the joint table")
    P = {}
    for st in ("S1", "S2"):
        for cls in CLS3:
            nt = set()
            for tool in TOOLS:
                r = s[(s["set"] == st) & (s.tool == tool) & (s.cls == cls)]
                rg = g[(g["set"] == st) & (g.tool == tool) & (g.cls == cls)].sort_values("seed")
                if len(r) != 1 or sorted(rg.seed) != [1, 2, 3]:
                    raise SystemExit(f"{STEM}: joint table lacks {st} {tool} {cls} (summary {len(r)}, genomes "
                                     f"{sorted(rg.seed)})")
                r = r.iloc[0]
                for _, x in rg.iterrows():             # per genome: definitions, and the call side = earlier table
                    if abs(x.tp_truth / x.n_truth - x.recall) > 6e-5:
                        raise SystemExit(f"{STEM}: D recall {st} {x.seed} {tool} {cls} != tp/n")
                    e = ext.loc[(st, int(x.seed), tool, cls)]
                    if (int(e.n_calls), int(e.tp_calls)) != (int(x.n_calls), int(x.tp_calls)):
                        raise SystemExit(f"{STEM}: D calls {st} {x.seed} {tool} {cls} differ from "
                                         "sim_nonauto_prf_ext.tsv")
                    if x.n_calls == 0:
                        if not (pd.isna(x.precision) and pd.isna(x.f1)):
                            raise SystemExit(f"{STEM}: D {st} {x.seed} {tool} {cls}: 0 calls but precision/F1")
                        continue
                    if abs(x.tp_calls / x.n_calls - x.precision) > 6e-5:
                        raise SystemExit(f"{STEM}: D precision {st} {x.seed} {tool} {cls} != tp_calls/n_calls")
                    if x.precision + x.recall == 0:      # P = R = 0: F1 undefined in the table
                        if not pd.isna(x.f1):
                            raise SystemExit(f"{STEM}: D {st} {x.seed} {tool} {cls}: P = R = 0 but F1 given")
                        continue
                    f1 = 2 * x.precision * x.recall / (x.precision + x.recall)
                    if abs(f1 - x.f1) > 1.2e-4:
                        raise SystemExit(f"{STEM}: D F1 {st} {x.seed} {tool} {cls} {x.f1} != 2PR/(P+R) {f1:.4f}")
                # summary = mean of the defined per-genome values, recomputed here from the integer counts (the
                # per-genome table is rounded to 4 digits); n summed
                ex_ = {"recall": [], "precision": [], "f1": []}
                for _, x in rg.iterrows():
                    rc = x.tp_truth / x.n_truth
                    ex_["recall"].append(rc)
                    if x.n_calls > 0:
                        pc = x.tp_calls / x.n_calls
                        ex_["precision"].append(pc)
                        if pc + rc > 0:
                            ex_["f1"].append(2 * pc * rc / (pc + rc))
                for key in ("precision", "recall", "f1"):
                    v = ex_[key]
                    if int(r[f"{key}_n_genomes"]) != len(v) or len(v) not in (0, 3):
                        raise SystemExit(f"{STEM}: D {st} {tool} {cls} {key}: {r[f'{key}_n_genomes']} genomes "
                                         f"(need 3, or 0 where no calls)")
                    if len(v) and abs(np.mean(v) - r[f"{key}_mean"]) > 5.01e-5:
                        raise SystemExit(f"{STEM}: D {st} {tool} {cls} {key} mean {r[f'{key}_mean']} != "
                                         f"{np.mean(v):.5f}")
                    if not len(v) and not pd.isna(r[f"{key}_mean"]):
                        raise SystemExit(f"{STEM}: D {st} {tool} {cls} {key}: no defined genome but a mean")
                if (int(r.n_truth), int(r.n_calls)) != (int(rg.n_truth.sum()), int(rg.n_calls.sum())):
                    raise SystemExit(f"{STEM}: D {st} {tool} {cls} n_truth / n_calls are not the genome sums")
                nt.add(int(r.n_truth))
                P[(st, tool, cls)] = dict(precision=r.precision_mean, recall=r.recall_mean, f1=r.f1_mean,
                                          n_truth=int(r.n_truth), tp_truth=int(rg.tp_truth.sum()),
                                          n_calls=int(r.n_calls), tp_calls=int(rg.tp_calls.sum()),
                                          n_genomes=int(r.precision_n_genomes))
            if len(nt) != 1:
                raise SystemExit(f"{STEM}: D {st} {cls}: n_truth differs between tools {sorted(nt)}")
    return P


def check_basis(S, P):
    """Legend D: 'Like Fig. 1B, D counts every simulated element, including copies of library families without
    annotated LTRs, which nearly every caller misses, so its recall runs below B's.' D's matching (one-to-one RO >= 0.9
    against every truth.tsv row, score_elements.match, as sim_class_prf_joint.py) is redone for the three S2 genomes
    and must give D_TABLE's per-genome n_truth / tp_truth exactly; then each part of the clause is asserted:
      every element: per genome, n_truth summed over D's classes = truth.tsv rows;
      families without annotated LTRs: copies of lib_class.tsv entries with has_ltr != 1 are in D's truth, B counts
        none of them (B's S2 elements = the std == 1 copies of the other entries, every class);
      nearly every caller misses them: S2 recall on those copies < 1/2 for every caller and > 1/3 for at most one;
      so its recall runs below B's: each caller's recall on those copies is below its recall on B's elements, and D's
        recall is below B's drawn recall (S1 + S2 pooled) and below B's counting in S2 (sim_gap.tsv) for every caller
        and class, except where B's recall is 0 (DANTE_LTR protein-less: both print 0.00)."""
    lib = {r["family"]: r for r in csv.DictReader(open(NA / "lib/lib_class.tsv"), delimiter="\t")}
    g = tsv(d_table_genomes())
    g = g[(g.system == "protein") & g.tool.isin(TOOLS)].set_index(["set", "seed", "tool", "cls"])
    bb = S["b_basis"]
    cnt = {}                       # (tool, origin) -> [n, tp] over the 3 genomes; origin old / ext / untagged
    for seed in (1, 2, 3):
        ds = f"S2_{seed}"
        truth = read_truth(str(B / f"sims/truth/{ds}.truth.tsv"))
        miss = sorted({t["family"] for t in truth if t["family"] not in lib})
        if miss:
            raise SystemExit(f"{STEM}: {ds} truth families missing from lib_class.tsv: {miss[:3]}")
        cls = [lib[t["family"]]["cls"] for t in truth]
        org = ["untagged" if lib[t["family"]]["has_ltr"] != "1" else ("old" if t["std"] == 1 else "ext") for t in truth]
        n_old = org.count("old")
        if n_old != bb["n"][ds]:
            raise SystemExit(f"{STEM}: {ds}: B counts {bb['n'][ds]} elements, std == 1 copies of entries with "
                             f"annotated LTRs number {n_old}")
        for tool in TOOLS:
            _, t2c = match(read_calls(str(B / f"calls/{ds}/{tool}.tsv")), truth, 0.9)
            for c in CLS3:
                idx = [i for i, x in enumerate(cls) if x == c]
                n, tp = len(idx), sum(i in t2c for i in idx)
                r = g.loc[("S2", seed, tool, c)]
                if (n, tp) != (int(r.n_truth), int(r.tp_truth)):
                    raise SystemExit(f"{STEM}: {D_TABLE.name} {ds} {tool} {c}: n/tp {int(r.n_truth)}/{int(r.tp_truth)} "
                                     f"!= every-element matching {n}/{tp}. The legend's D basis clause says D counts "
                                     "every simulated element: reword it for this table")
            for o in ("old", "ext", "untagged"):
                idx = [i for i, x in enumerate(org) if x == o]
                k = cnt.setdefault((tool, o), [0, 0])
                k[0] += len(idx)
                k[1] += sum(i in t2c for i in idx)
        if sum(int(g.loc[("S2", seed, "ltrquest", c)].n_truth) for c in CLS3) != len(truth):
            raise SystemExit(f"{STEM}: {ds}: D's classes do not sum to every truth.tsv row")
        num("legend/D (basis)", ds, "", "elements counted by D (n_truth summed over classes) = truth.tsv rows",
            len(truth), f"{d_table_genomes().relative_to(B)}; sims/truth/{ds}.truth.tsv")
    if bb["families"] & {f for f, r in lib.items() if r["has_ltr"] != "1"}:
        raise SystemExit(f"{STEM}: B's elements include copies of entries without annotated LTRs")
    nu = cnt[("ltrquest", "untagged")][0]
    if nu == 0:
        raise SystemExit(f"{STEM}: no copies of entries without annotated LTRs in S2 (legend D clause)")
    num("legend/D (basis)", "S2", "", "copies of library families without annotated LTRs (has_ltr != 1), 3 genomes; "
        "B counts none of them", nu, "experiments/nonauto/lib/lib_class.tsv; sims/truth/S2_*.truth.tsv")
    num("legend/D (basis)", "S2", "", "elements B counts (std == 1 copies of entries with annotated LTRs, all classes)",
        cnt[("ltrquest", "old")][0], "sim_elements.tsv; sims/truth/S2_*.truth.tsv")
    rec_u = {t: cnt[(t, "untagged")][1] / cnt[(t, "untagged")][0] for t in TOOLS}
    rec_o = {t: cnt[(t, "old")][1] / cnt[(t, "old")][0] for t in TOOLS}
    for t in TOOLS:
        num("legend/D (basis)", "S2", t, "recall on copies of families without annotated LTRs (found / copies)",
            f"{rec_u[t]:.4f} ({cnt[(t, 'untagged')][1]} / {cnt[(t, 'untagged')][0]})",
            "calls/S2_*/<tool>.tsv vs sims/truth (score_elements.match, RO >= 0.9, every truth row)")
        num("legend/D (basis)", "S2", t, "recall on the elements B counts (found / elements)",
            f"{rec_o[t]:.4f} ({cnt[(t, 'old')][1]} / {cnt[(t, 'old')][0]})",
            "calls/S2_*/<tool>.tsv vs sims/truth (score_elements.match, RO >= 0.9, every truth row)")
    if max(rec_u.values()) >= 0.5 or sum(v > 1 / 3 for v in rec_u.values()) > 1:
        raise SystemExit(f"{STEM}: 'nearly every caller misses them' fails: {rec_u}")
    if any(rec_u[t] >= rec_o[t] for t in TOOLS):
        raise SystemExit(f"{STEM}: a caller finds copies without annotated LTRs as often as B's elements: {rec_u}")
    gap = S["gap"].set_index(["set", "tool"])
    for i, t in enumerate(TOOLS):
        for c in CLS2:
            d = P[("S2", t, c)]["recall"]
            b = float(S[c]["recall"][i])                                    # drawn in B (S1 + S2 pooled)
            b2 = float(gap.loc[("S2", t), "recall_coding" if c == "coding" else "recall_proteinless"])   # B's, S2
            num("legend/D (basis)", "S2", t, f"{CLS3_NAME[c]} recall, D (S2) vs B as drawn (S1+S2) vs B's counting "
                "in S2", f"{d:.4f} vs {b:.4f} vs {b2:.4f}", f"{D_SRC}; sim_elements.tsv; sim_gap.tsv")
            if round(d, 2) > round(b, 2) or (b >= 0.01 and not (d < b and d < b2)):
                raise SystemExit(f"{STEM}: 'so its recall runs below B's' fails for {t} {c} ({d}, {b}, {b2})")
    return rec_u


# ============================================================================ data: real genomes (E)


def load_missing():
    """tool x genome cells drawn as missing, from report/tables/resources.tsv (asserted against the expected list)"""
    r = tsv(B / "report/tables/resources.tsv")
    r = r[r.ds.isin(DS) & r.tool.isin(TOOLS)]
    if len(r) != len(DS) * len(TOOLS):
        raise SystemExit(f"{STEM}: resources.tsv lacks real-genome tool rows ({len(r)})")
    bad = set(r.status) - FAILED_STATUS - {"ok"}
    if bad:
        raise SystemExit(f"{STEM}: resources.tsv status not mapped: {bad}")
    failed = {(d, t) for d, t, st in zip(r.ds, r.tool, r.status) if st in FAILED_STATUS}
    if failed != EXPECT_FAILED:
        raise SystemExit(f"{STEM}: failed runs {sorted(failed)} differ from the expected list")
    MISSING.clear()
    for d, t, st in zip(r.ds, r.tool, r.status):
        if st in FAILED_STATUS:
            MISSING[(d, t)] = "tool failed"
            num("E", d, t, "drawn as 'tool failed'", st, "report/tables/resources.tsv (status)")


def load_real():
    """E: support of protein-less calls (real_calls.tsv); gorilla protein-less share of loci for the legend"""
    loci = tsv(NA / "real_loci.tsv")
    calls = tsv(NA / "real_calls.tsv")
    lt = pd.read_csv(NA / "real_locus_table.tsv", sep="\t", usecols=["ds", "cls", "len"])
    out = {"comp": {}, "prec": {}, "coding_share": {}}
    for ds in DS:
        r = loci[(loci.ds == ds) & (loci.tool == "ltrquest") & (loci.scope == "all")].set_index("cls").n_loci
        tot = int(r.sum())
        out["comp"][ds] = {c: int(r[c]) / tot for c in ("coding", "partial", "proteinless")}
        out["comp"][ds]["n"] = tot
        out["comp"][ds]["n_pl"] = int(r["proteinless"])
        if len(lt[(lt.ds == ds) & (lt.cls == "proteinless")]) != int(r["proteinless"]):
            raise SystemExit(f"{STEM}: {ds} protein-less loci differ between real_locus_table and real_loci.tsv")
        for tool in TOOLS:
            p = B / f"eval/adj/{ds}/{tool}.calls.tsv"
            if (ds, tool) in MISSING:
                if p.exists() and p.stat().st_size > 0 and tool != "dante_ltr":
                    raise SystemExit(f"{STEM}: {ds} {tool} has calls but is marked missing")
                continue
            pr = calls[(calls.ds == ds) & (calls.tool == tool)].set_index("cls")
            out["prec"][(ds, tool)] = pr.loc["proteinless"].to_dict()
            out["coding_share"][(ds, tool)] = float(pr.loc["coding", "supported_share"])
    return out


# ============================================================================ data: the example (A)
def names_map():
    return {r["original_id"]: r["new_id"] for r in csv.DictReader(open(B / "data/rice/names.tsv"), delimiter="\t")}


def rm_hits(c, s, e):
    """curated RepeatMasker hits overlapping [s, e): .out is 1-based closed -> 0-based half-open"""
    out = []
    with open(B / "mask/rice/curated/rice.fa.out") as f:
        for line in f:
            t = line.split()
            if len(t) > 10 and t[0].isdigit() and t[4] == c:
                a, b = int(t[5]) - 1, int(t[6])
                if b > s and a < e:
                    out.append(dict(s=a, e=b, strand=t[8], name=t[9], cls=t[10], div=float(t[1])))
    return out


def gff_records(c, s, e):
    """all DANTE protein_domain records overlapping [s, e): GFF3 1-based closed -> 0-based half-open"""
    out = []
    with open(B / "runs/rice/dante_ltr/dante.gff3") as f:
        for line in f:
            if line.startswith("#"):
                continue
            t = line.rstrip("\n").split("\t")
            if len(t) < 9 or t[0] != c or t[2] != "protein_domain":
                continue
            a, b = int(t[3]) - 1, int(t[4])
            if b > s and a < e:
                fc = re.search(r"Final_Classification=([^;]+)", t[8]).group(1)
                out.append(dict(s=a, e=b, strand=t[6], name=re.search(r"Name=([^;]+)", t[8]).group(1), cls=fc))
    return out


def tool_calls(c, s, e):
    """every call of every tool overlapping [s, e) (calls/rice/<tool>.tsv: 0-based half-open, original names)"""
    nm = names_map()
    out = {}
    for tool in TOOLS:
        rows = []
        for r in csv.DictReader(open(B / f"calls/rice/{tool}.tsv"), delimiter="\t"):
            if nm.get(r["chrom"], r["chrom"]) != c:
                continue
            a, b = int(r["start"]), int(r["end"])
            if b > s and a < e:
                lt = int(r["lltr_end"]) if r["lltr_end"] not in (".", "") else None
                rt = int(r["rltr_start"]) if r["rltr_start"] not in (".", "") else None
                rows.append(dict(s=a, e=b, lltr_end=lt, rltr_start=rt, id=r["id"]))
        out[tool] = rows
    return out


def load_example():
    idx = proxies.fai(str(B / "data/rice/genome.fa"))
    fh = open(B / "data/rice/genome.fa", "rb")
    cur = proxies.ref_ltr(str(B / "mask/rice/curated/rice.fa.out"))
    cst = {k: [x[0] for x in v] for k, v in cur.items()}
    lt = pd.read_csv(NA / "real_locus_table.tsv", sep="\t", dtype={"locus": str})
    lt = lt[lt.ds == "rice"]
    res = []
    for ex in EX:
        c, s, e = ex["chrom"], ex["s"], ex["e"]
        row = lt[lt.locus == ex["locus"]].iloc[0]
        calls = tool_calls(c, s - FLANK, e + FLANK)
        lq = [x for x in calls["ltrquest"] if x["s"] == s and x["e"] == e]
        tsd_call = next(r for r in csv.DictReader(open(B / "calls/rice/ltrquest.ids.tsv"), delimiter="\t")
                        if r["chrom"] == c and int(r["start"]) == s and int(r["end"]) == e)["tsd"]
        k = len(tsd_call)
        seq_l = proxies.fetch(fh, idx, c, s - 14, s + 14)
        seq_r = proxies.fetch(fh, idx, c, e - 14, e + 14)
        hits = rm_hits(c, s - FLANK, e + FLANK)
        recs = gff_records(c, s, e)
        ltr_recs = [d for d in recs if d["cls"].startswith("Class_I|LTR")]
        cov = proxies.covered(cur, cst, c, s, e) / (e - s)
        chk = dict(locus_row=row.cls == ex["cls"] and row.precise == 1 and "ltrquest" in row.tools.split(","),
                   ltrquest_call=len(lq) == 1, ltr_bounds=bool(lq) and lq[0]["lltr_end"] is not None
                   and lq[0]["rltr_start"] is not None,
                   tsd=4 <= k <= 6 and seq_l[14 - k:14] == tsd_call and seq_r[14:14 + k] == tsd_call,
                   tgca=seq_l[14:16] == "TG" and seq_r[12:14] == "CA", curated_cov=cov >= 0.8)
        if ex["cls"] == "proteinless":
            chk["no_ltr_domain"] = len(ltr_recs) == 0
        else:
            names = [d["name"] for d in sorted(ltr_recs, key=lambda d: d["s"])]
            ref = ["GAG", "PROT", "RT", "RH", "aRH", "INT", "CHD"]
            chk["one_strand"] = len({d["strand"] for d in ltr_recs}) == 1
            chk["one_superfamily"] = len({d["cls"].split("|")[2] for d in ltr_recs}) == 1
            chk["domain_order"] = all(n in ref for n in names) and \
                [ref.index(n) for n in names] == sorted(ref.index(n) for n in names)
        if not all(chk.values()):
            raise SystemExit(f"{STEM}: example {ex['key']} fails a check: {chk}")
        ex = dict(ex, calls=calls, lq=lq[0], tsd=tsd_call, seq_l=seq_l, seq_r=seq_r,
                  hits=[h for h in hits if proxies.LTR.search(h["cls"])], other_hits=[h for h in hits if not
                                                                                     proxies.LTR.search(h["cls"])],
                  doms=ltr_recs, other_recs=[d for d in recs if not d["cls"].startswith("Class_I|LTR")],
                  cov=cov, checks=chk, row=row)
        log(f"example {ex['key']}: {c}:{s}-{e} ({e - s} bp) checks {chk}; tools {row.tools}; "
            f"K2P {row.k2p}; LTR hits {len(ex['hits'])}; other DANTE records {len(ex['other_recs'])}")
        res.append(ex)
    fh.close()
    return res


# ============================================================================ drawing helpers
def ax_mm(fig, x, y, w, h, xlim=None, ylim=None, off=True):
    ax = poafig.mm_axes(fig, x, y, w, h)
    ax.set_xlim(*(xlim or (0, w)))
    ax.set_ylim(*(ylim or (h, 0)))
    if off:
        ax.axis("off")
    return ax


def failed_box(ax, x, y, w, h, fs=6.0):
    """tool side of the set-wide missing-result rule: dark gray, white 'tool failed'"""
    ax.add_patch(Rectangle((x, y), w, h, fc=FAILC, ec="none", lw=0))
    ax.text(x + w / 2, y + h / 2, "tool failed", ha="center", va="center", fontsize=fs, color="white")


def tool_labels(fig, x_right, y, rows_y, rh):
    """tool names (LTRquest bold) with a color square, right-aligned to x_right; rows_y = row tops (mm)"""
    ax = ax_mm(fig, 0, y, x_right, max(rows_y) + rh - y)
    for t, yy in zip(TOOLS, rows_y):
        yc = yy - y + rh / 2
        ax.text(x_right - 2.6, yc, NAME[t], ha="right", va="center", fontsize=FS,
                fontweight="bold" if t == "ltrquest" else "normal")
        ax.add_patch(Rectangle((x_right - 2.0, yc - 0.6), 1.2, 1.2, fc=COL[t], ec="none"))
    return ax


def text_w_mm(ax, t):
    """width of a text artist in the axes' data units (mm for ax_mm axes)"""
    r = ax.figure.canvas.get_renderer()
    bb = t.get_window_extent(renderer=r)
    inv = ax.transData.inverted()
    (a, _), (b, _) = inv.transform([(bb.x0, bb.y0), (bb.x1, bb.y0)])
    return abs(b - a)


def letter(fig, x, y, L):
    poafig.panel_letter(fig, x, y, L)


def fmt_int(v):
    return f"{int(v):,}"


# ============================================================================ panel A
def mono_row(ax, x, y, segs, cw, fs=6.0):
    """sequence text at fixed pitch cw (mm); segs: [(text, style)] style plain / tsd / ltr"""
    for text, style in segs:
        if style == "tsd":
            ax.add_patch(Rectangle((x - 0.1, y - 1.35), cw * len(text) + 0.2, 2.7, fc="none", ec=INK, lw=0.5,
                                   zorder=3))
        if style == "ltr":
            ax.add_patch(Rectangle((x, y - 1.3), cw * len(text), 2.6, fc=LTRC, ec="none", zorder=1))
        for ch in text:
            ax.text(x + cw / 2, y, ch, ha="center", va="center", fontsize=fs, family="Liberation Mono",
                    color="white" if style == "ltr" else (MUTED if style == "plain" else INK), zorder=4)
            x += cw
    return x


def panel_ends(fig, x0, y0, ex):
    """both ends of the protein-less copy: 9 bp flank + TSD + 12 bp LTR | 12 bp LTR + TSD + 9 bp flank; the element
    is named in-panel (protein-less Dasheng copy, rice Chr1 2.008 Mb)"""
    H = 13.0
    ax = ax_mm(fig, 0, y0, 180, H)
    k = len(ex["tsd"])
    sl, sr = ex["seq_l"], ex["seq_r"]                 # [s-14, s+14), [e-14, e+14)
    cw = 1.34
    yS = 9.6
    nfl = 9
    left = [(sl[14 - k - nfl:14 - k], "plain"), (sl[14 - k:14], "tsd"), (sl[14:26], "ltr")]
    right = [(sr[2:14], "ltr"), (sr[14:14 + k], "tsd"), (sr[14 + k:14 + k + nfl], "plain")]
    t1 = ax.text(x0, 1.6, "protein-less", ha="left", va="center", fontsize=FS, fontweight="bold")
    ax.text(x0 + text_w_mm(ax, t1) + 1.2, 1.6, f"{ex['name']} copy", ha="left", va="center", fontsize=FS,
            color=MUTED)
    ax.text(x0 - 2.6, yS - 1.3, f"rice {ex['chrom_name']}", ha="right", va="center", fontsize=FS)
    ax.text(x0 - 2.6, yS + 1.5, f"{ex['s'] / 1e6:.3f} Mb", ha="right", va="center", fontsize=FS, color=MUTED)
    xa = x0
    xe = mono_row(ax, xa, yS, left, cw)
    ax.text(xa + cw * (nfl + k / 2), yS - 2.2, "TSD", ha="center", va="bottom", fontsize=FSS)
    ax.text(xa + cw * (nfl + k + 6), yS - 2.2, "LTR start", ha="center", va="bottom", fontsize=FSS)
    mid = f"… {ex['e'] - ex['s'] - 24:,} bp …"
    tm = ax.text(xe + 1.4, yS, mid, ha="left", va="center", fontsize=FSS, color=MUTED)
    xr = xe + 1.4 + text_w_mm(ax, tm) + 1.6
    xe2 = mono_row(ax, xr, yS, right, cw)
    ax.text(xr + cw * 6, yS - 2.2, "LTR end", ha="center", va="bottom", fontsize=FSS)
    ax.text(xr + cw * (12 + k / 2), yS - 2.2, "TSD", ha="center", va="bottom", fontsize=FSS)
    num("A", "rice", "", f"{ex['name']} position, {ex['chrom_name']} 1-based start (Mb shown)", ex["s"] + 1,
        "real_locus_table.tsv / calls/rice/ltrquest.tsv (0-based start + 1)")
    num("A", "rice", "", f"{ex['name']} length (bp)", ex["e"] - ex["s"], "calls/rice/ltrquest.tsv")
    num("A", "rice", "", "TSD", ex["tsd"], "data/rice/genome.fa [s-5, s) and [e, e+5), 0-based")
    num("A", "rice", "", "bases between the drawn LTR ends", ex["e"] - ex["s"] - 24, "LTRquest call span - 24")
    num("A", "rice", "", "member tools of the locus", ex["row"].tools, "real_locus_table.tsv (tools)")
    return y0 + H, xe2


# ============================================================================ panels B, C, D (shared tool rows)
def panel_tools(fig, y0, S, P, geo):
    LAB = geo["LAB"]
    RH, rgap = geo["RH"], geo["rgap"]
    hdr = geo["hdr"]
    yR = y0 + hdr
    rows_y = [yR + i * RH for i in range(9)]
    tool_labels(fig, LAB, yR, rows_y, RH - rgap)
    # ---------------- B: every simulated element, found or missed
    bw, bgap, vtxt = geo["B_bw"], geo["B_gap"], geo["B_txt"]
    xB = [LAB + 0.5, LAB + 0.5 + bw + vtxt + bgap]
    for j, cls in enumerate(CLS2):
        d = S[cls]
        N = len(d["L"])
        xb = xB[j]
        title = "protein-less" if cls == "proteinless" else "coding"
        tax = ax_mm(fig, xb, y0, bw, hdr)
        tax.text(bw / 2, 1.6, f"{title} elements", ha="center", va="center", fontsize=FS)
        tax.text(bw / 2, 4.4, f"n = {N:,}", ha="center", va="center", fontsize=FS, color=MUTED)
        num("B", "S1+S2", "", f"{CLS3_NAME[cls]} simulated elements (3 seeds x 2 regimes)", N, "sim_elements.tsv")
        # 2 / 4 kb: dashed rank lines through the rows, labeled above the first row (no key beside B)
        xk = {}
        for kb in (2000, 4000):
            r = int(np.searchsorted(d["L"], kb))
            xk[kb] = r / N * bw
            tax.plot([xk[kb], xk[kb]], [hdr - 2.2, hdr], color=INK, lw=0.5, ls=(0, (2, 1.4)))
            num("B", "S1+S2", "", f"{CLS3_NAME[cls]} elements shorter than {kb // 1000} kb", r,
                "sim_elements.tsv len")
        tax.text(xk[2000], hdr - 2.5, "2", ha="center", va="bottom", fontsize=FST)
        tax.text(xk[4000] - 0.65, hdr - 2.5, "4 kb", ha="left", va="bottom", fontsize=FST)
        cbest = max(round(float(x), 2) for x in d["recall"])     # set-wide rule: column best bold, ties all
        for i, t in enumerate(TOOLS):
            rax = poafig.mm_axes(fig, xb, rows_y[i], bw, RH - rgap)
            rgb = np.array([to_rgb(COL[t]), to_rgb(MISS_SIM)])
            img = rgb[1 - d["hits"][i]][None, :, :]
            rax.imshow(img, aspect="auto", extent=(0, N, 0, 1), interpolation="antialiased",
                       interpolation_stage="rgba", resample=True)
            rax.set_xlim(0, N)
            rax.set_ylim(0, 1)
            rax.axis("off")
            for kb in (2000, 4000):
                r = int(np.searchsorted(d["L"], kb))
                rax.axvline(r, color="white" if t in ("hite",) else INK, lw=0.5, ls=(0, (2, 1.4)))
            v = d["recall"][i]
            num("B", "S1+S2", t, f"{CLS3_NAME[cls]} recall (pooled)", round(float(v), 3), "sim_elements.tsv")
            vax = ax_mm(fig, xb + bw, rows_y[i], vtxt, RH - rgap)
            vax.text(0.7, (RH - rgap) / 2, f"{v:.2f}", ha="left", va="center", fontsize=FSS,
                     color=INK, fontweight="bold" if round(float(v), 2) == cbest else "normal")
        bax = ax_mm(fig, xb, rows_y[-1] + RH, bw, 6)
        bax.text(bw / 2, 1.4, "shortest → longest", ha="center", va="center", fontsize=FSS, color=MUTED)
    # ---------------- C: recall vs length
    xC, wC = geo["xC"], geo["wC"]
    kax = poafig.key_axes(fig, xC, y0, wC, hdr - 1.0)
    poafig.hkey(kax, [(poafig.g_line(INK, w=3.6, lw=0.9), "coding")], 0.4, 1.6, wC, fontsize=FSS)
    poafig.hkey(kax, [(poafig.g_line(INK, w=3.6, lw=0.9, ls=(0, (2.2, 1.2))), "protein-less")], 0.4, 4.4, wC,
                fontsize=FSS)
    poafig.hkey(kax, [(poafig.g_box("#D8D8D8", w=3.6, h=1.3), "95% CI")], 0.4, 7.2, wC, fontsize=FSS)
    for i, t in enumerate(TOOLS):
        ax = poafig.mm_axes(fig, xC, rows_y[i], wC, RH - rgap)
        for cls in CLS2:
            d = S[cls]
            ed = d["edges"] / 1000
            lw = 1.1 if t == "ltrquest" else 0.8
            ls = "-" if cls == "coding" else (0, (2.2, 1.2))
            ax.stairs(d["hi"][i], ed, baseline=d["lo"][i], fill=True, color=tint(COL[t], 0.22), lw=0)
            ax.stairs(d["rec"][i], ed, color=COL[t], lw=lw, ls=ls, baseline=None)
            for b in range(len(ed) - 1):
                num("C", "S1+S2", t, f"{CLS3_NAME[cls]} recall {d['edges'][b]}-{d['edges'][b + 1]} bp "
                    f"[{d['lo'][i][b]:.3f}, {d['hi'][i][b]:.3f}]", round(float(d["rec"][i][b]), 3),
                    "sim_elements.tsv (family bootstrap here)")
        for kb in (2, 4):
            ax.axvline(kb, color=INK, lw=0.5, ls=(0, (2, 1.4)), zorder=0)
        ax.set_xscale("log")
        ax.set_xlim(0.35, 20)
        ax.set_ylim(-0.04, 1.08)
        ax.set_yticks([0, 1])
        ax.set_yticklabels(["0", "1"], fontsize=FST)
        ax.tick_params(axis="both", length=1.5, pad=1)
        ax.xaxis.set_minor_locator(mpl.ticker.NullLocator())
        ax.set_xticks([0.5, 1, 2, 4, 8, 16])
        if i < 8:
            ax.set_xticklabels([])
        else:
            ax.set_xticklabels(["0.5", "1", "2", "4", "8", "16"], fontsize=FST)
            ax.set_xlabel("element length (kb)", fontsize=FS, labelpad=1.5)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    lax = ax_mm(fig, xC - 6.4, yR, 3, 9 * RH)
    lax.text(1.5, 9 * RH / 2, "elements found (recall)", rotation=90, ha="center", va="center", fontsize=FS)
    # ---------------- D: scorecard, S2 (new)
    xD, cwD, ggD = geo["xD"], geo["cwD"], geo["ggD"]
    wD = 9 * cwD + 2 * ggD
    hax = ax_mm(fig, xD, y0, wD, hdr)
    hax.text(wD / 2, 1.6, "S2, mean of three genomes", ha="center", va="center", fontsize=FS)
    src = D_SRC
    for g, cls in enumerate(CLS3):
        xg = g * (3 * cwD + ggD)
        hax.text(xg + 1.5 * cwD, 4.6, CLS3_NAME[cls], ha="center", va="center", fontsize=FS)
        n_el = int(P[("S2", "ltrquest", cls)]["n_truth"])
        hax.text(xg + 1.5 * cwD, 7.1, f"n = {n_el:,}", ha="center", va="center", fontsize=FS, color=MUTED)
        num("D", "S2", "", f"{CLS3_NAME[cls]} simulated elements (n_truth summed over 3 genomes)", n_el, src)
        hax.plot([xg + 0.3, xg + 3 * cwD - 0.3], [8.5, 8.5], color=INK, lw=0.4, solid_capstyle="butt")
        for m, (_, lab) in enumerate(METRICS):
            hax.text(xg + (m + 0.5) * cwD, hdr - 2.3, lab, ha="center", va="center", fontsize=FS, color=MUTED)
    pad = 0.16
    hcell = RH - rgap
    # set-wide rule: bold the best printed (2-decimal) value per column among the nine tools, ties all bold, only
    # where >= 2 tools are scored; undefined cells never count
    best = {}
    for cls in CLS3:
        for key, _ in METRICS:
            vals = [round(float(P[("S2", t, cls)][key]), 2) for t in TOOLS if not pd.isna(P[("S2", t, cls)][key])]
            best[(cls, key)] = max(vals) if len(vals) >= 2 else None
            if best[(cls, key)] is not None:
                num("D", "S2", ", ".join(NAME[t] for t in TOOLS if not pd.isna(P[("S2", t, cls)][key]) and
                                         round(float(P[("S2", t, cls)][key]), 2) == best[(cls, key)]),
                    f"{CLS3_NAME[cls]} {key}: best printed value (bold)", f"{best[(cls, key)]:.2f}", src)
    for i, t in enumerate(TOOLS):
        ax = ax_mm(fig, xD, rows_y[i], wD, hcell)
        for g, cls in enumerate(CLS3):
            r = P[("S2", t, cls)]
            xg = g * (3 * cwD + ggD)
            for m, (key, _) in enumerate(METRICS):
                v = r[key]
                num("D", "S2", t, f"{CLS3_NAME[cls]} {key} (mean of 3 genomes; summed n_truth {r['n_truth']}, "
                    f"tp_truth {r['tp_truth']}, n_calls {r['n_calls']}, tp_calls {r['tp_calls']})",
                    "NA (0 calls)" if pd.isna(v) else round(float(v), 4), src)
                if pd.isna(v):
                    continue
                x = xg + m * cwD
                isbest = best[(cls, key)] is not None and round(float(v), 2) == best[(cls, key)]
                ax.add_patch(Rectangle((x + pad, pad), cwD - 2 * pad, hcell - 2 * pad, fc=CMAP(float(v)), ec="none"))
                ax.text(x + cwD / 2, hcell / 2, f"{v:.2f}", ha="center", va="center", fontsize=FS,
                        color=text_color(float(v)), fontweight="bold" if isbest else "normal")
            if r["n_calls"] == 0:                      # precision and F1 undefined (no calls of that class from
                for m, (key, _) in enumerate(METRICS):     # a finished run): plain light gray with a dash, no hatch
                    if pd.isna(r[key]):
                        x = xg + m * cwD
                        ax.add_patch(Rectangle((x + pad, pad), cwD - 2 * pad, hcell - 2 * pad, fc=UNDEFC, ec="none"))
                        ax.text(x + cwD / 2, hcell / 2, "–", ha="center", va="center", fontsize=FS, color=MUTED)
    # color key under the scorecard: 'score 0 [bar] 1'
    kax = ax_mm(fig, xD, rows_y[-1] + RH + 0.2, wD, 3.2)
    barw = 14.0
    x1 = wD - 2.2
    x0b = x1 - barw
    gg = np.linspace(0, 1, 128)
    for k in range(len(gg) - 1):
        kax.add_patch(Rectangle((x0b + barw * gg[k], 1.0), barw / (len(gg) - 1) + 0.02, 1.2, fc=CMAP(gg[k]),
                                ec="none"))
    kax.add_patch(Rectangle((x0b, 1.0), barw, 1.2, fc="none", ec=INK, lw=0.3))
    kax.text(x1 + 0.6, 1.6, "1", ha="left", va="center", fontsize=FS)
    kax.text(x0b - 0.6, 1.6, "0", ha="right", va="center", fontsize=FS)
    kax.text(x0b - 2.4, 1.6, "score", ha="right", va="center", fontsize=FS)
    return rows_y[-1] + RH + 6.6


# ============================================================================ panel E
def panel_E(fig, y0, R, geo):
    LAB = geo["LAB"]
    cw, cg = geo["cw"], geo["cg"]
    xs = [LAB + 0.5 + j * (cw + cg) for j in range(5)]
    RHh = geo["RHh"]
    wbar = cw - 7.6
    kax = poafig.key_axes(fig, xs[0], y0, 157, 4.0)

    def g_chance(ax, x, y, measure=False):
        if not measure:
            ax.add_patch(Rectangle((x, y - 1.1), 1.2, 2.2, fc=NULLC, ec="none", alpha=0.45))
            ax.plot([x + 0.6, x + 0.6], [y - 1.1, y + 1.1], color="#4D4D4D", lw=0.6)
        return 1.2

    def g_tri(ax, x, y, measure=False):
        if not measure:
            ax.plot([x + 0.6], [y], marker="v", ms=2.3, color=INK, mew=0, ls="none")
        return 1.2

    poafig.hkey(kax, [(poafig.g_box("#5E5E5E", w=2.4, h=1.2), "protein-less calls, library-supported"),
                      (poafig.g_box("#E6E6E6", w=2.4, h=1.2), "unsupported"),
                      (g_chance, "chance (same spans moved)"), (g_tri, "same caller's coding calls")],
                0, 1.6, 157, fontsize=FSS, gap=3.4)
    yh = y0 + 4.8
    hax = ax_mm(fig, 0, yh, 180, 3.0)
    for j, ds in enumerate(DS):
        hax.text(xs[j] + wbar / 2, 1.5, DSNAME[ds], ha="center", va="center", fontsize=FS)
        hax.text(xs[j] + cw - 0.2, 1.5, "calls", ha="right", va="center", fontsize=FSS, color=MUTED)
    hh_rows = [yh + 3.6 + i * RHh for i in range(9)]
    tool_labels(fig, LAB, hh_rows[0], hh_rows, RHh - 0.4)
    for j, ds in enumerate(DS):
        for i, t in enumerate(TOOLS):
            ax = poafig.mm_axes(fig, xs[j], hh_rows[i], wbar, RHh - 0.4)
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.axis("off")
            nax = ax_mm(fig, xs[j] + wbar, hh_rows[i], cw - wbar, RHh - 0.4)
            if (ds, t) in MISSING:
                failed_box(ax, 0, 0.05, 1, 0.9)
                ax.set_xlim(0, 1)
                continue
            p = R["prec"][(ds, t)]
            n = int(p["n_calls"])
            nax.text(cw - wbar - 0.2, (RHh - 0.4) / 2, fmt_int(n), ha="right", va="center", fontsize=FSS,
                     color=INK)
            num("E", ds, t, "protein-less calls", n, "real_calls.tsv (n_calls)")
            if n == 0:
                ax.text(0.02, 0.45, "no protein-less calls", ha="left", va="center", fontsize=FSS, color=MUTED)
                continue
            sh = float(p["supported_share"])
            y0b, y1b = 0.12, 0.62
            ax.add_patch(Rectangle((0, y0b), 1, y1b - y0b, fc=tint(COL[t], 0.14), ec="none"))
            ax.add_patch(Rectangle((0, y0b), sh, y1b - y0b, fc=COL[t], ec="none"))
            if str(p["sup_lo"]) not in ("nan", "NA"):
                ax.plot([float(p["sup_lo"]), float(p["sup_hi"])], [(y0b + y1b) / 2] * 2,
                        color=NULLC if t == "hite" else INK, lw=0.5, solid_capstyle="butt")
            nl, nh, ns_ = float(p["null_lo"]), float(p["null_hi"]), float(p["null_supported_share"])
            ax.add_patch(Rectangle((nl, 0.0), max(nh - nl, 0.006), 0.74, fc=NULLC, ec="none", alpha=0.45, zorder=4))
            ax.plot([ns_, ns_], [0.0, 0.74], color="#4D4D4D", lw=0.6, zorder=5)
            cs = R["coding_share"][(ds, t)]
            ax.plot([cs], [0.86], marker="v", ms=2.2, color=INK, mew=0, ls="none", zorder=6, clip_on=False)
            num("E", ds, t, f"protein-less supported share [{p['sup_lo']}, {p['sup_hi']}]", round(sh, 3),
                "real_calls.tsv")
            num("E", ds, t, f"chance (moved spans) [{nl}, {nh}]", round(ns_, 3), "real_calls.tsv (null)")
            num("E", ds, t, "protein-less library-supported calls", int(p["supported"]), "real_calls.tsv")
            num("E", ds, t, "protein-less supported share / chance", round(sh / ns_, 2) if ns_ > 0 else "NA",
                "real_calls.tsv (ratio computed here)")
            num("E", ds, t, "coding calls supported share", round(cs, 3), "real_calls.tsv")
        bax = poafig.mm_axes(fig, xs[j], hh_rows[-1] + RHh - 0.4, wbar, 0.01)
        bax.set_xlim(0, 1)
        bax.set_xticks([0, 0.5, 1])
        bax.set_xticklabels(["0", "0.5", "1"], fontsize=FST)
        bax.tick_params(axis="x", length=1.5, pad=1)
        bax.set_yticks([])
        for sp in ("top", "right", "left"):
            bax.spines[sp].set_visible(False)
        if j == 2:
            bax.set_xlabel("library-supported share", fontsize=FS, labelpad=1.5)
    return hh_rows[-1] + RHh + 7.0


# ============================================================================ legend
def extra_numbers(S, P, R):
    """numbers the legend rests on that are not drawn as such: per-set recall (sim_gap.tsv), the S1 scorecard,
    LTRquest's d50 in the aging series (div_d50.tsv; drawn in Fig. 1C-D), the precision comparison of D, the
    real-genome support ratios"""
    g = S["gap"]
    for _, r in g.iterrows():
        for cls, k in (("coding", "recall_coding"), ("proteinless", "recall_proteinless")):
            num("legend/B", r["set"], r["tool"], f"{CLS3_NAME[cls]} recall (3 seeds pooled) [{r[k + '_lo']}, "
                f"{r[k + '_hi']}]", r[k], "sim_gap.tsv (family bootstrap CI)")
        num("legend/C", r["set"], r["tool"], f"coding recall at protein-less lengths [{r['recall_coding_lenmatched_lo']}, "
            f"{r['recall_coding_lenmatched_hi']}]", r["recall_coding_lenmatched"], "sim_gap.tsv")
    src = D_SRC
    # S1 scorecard (not drawn; legend: 'the moderate regime (S1) gives the same picture')
    for t in TOOLS:
        for cls in CLS3:
            r = P[("S1", t, cls)]
            for key, _ in METRICS:
                v = r[key]
                num("legend/D (S1, not drawn)", "S1", t, f"{CLS3_NAME[cls]} {key} (mean of 3 genomes; summed "
                    f"n_truth {r['n_truth']}, n_calls {r['n_calls']})",
                    "NA (0 calls)" if pd.isna(v) else round(float(v), 4), src)
    # D legend: LTRquest protein-less precision 0.94 vs LTR_retriever, LTR_Finder, HiTE 0.99-1.00 (S2), so its lead
    # there is recall; DANTE_LTR no partly coding / protein-less calls; 'S1 gives the same picture' = in S1 too,
    # LTRquest's protein-less precision is below all three, LTRquest has the highest recall and F1 in every class
    # (S2's bold pattern) and DANTE_LTR makes no partly coding / protein-less calls
    strict = ("ltr_retriever", "ltr_finder", "hite")
    for st in ("S2", "S1"):
        lq = P[(st, "ltrquest", "proteinless")]["precision"]
        sp = [P[(st, t, "proteinless")]["precision"] for t in strict]
        if st == "S2" and not (round(lq, 2) == 0.94 and all(0.985 <= v <= 1.0 for v in sp)):
            raise SystemExit(f"{STEM}: legend D precision sentence fails ({lq}, {sp})")
        if not lq < min(sp):
            raise SystemExit(f"{STEM}: {st}: LTRquest protein-less precision {lq} not below {sp}")
        for cls in CLS3:
            for key in ("recall", "f1"):
                v = {t: P[(st, t, cls)][key] for t in TOOLS if not pd.isna(P[(st, t, cls)][key])}
                top = sorted(v, key=v.get, reverse=True)
                if top[0] != "ltrquest" or round(v[top[1]], 2) >= round(v["ltrquest"], 2):
                    raise SystemExit(f"{STEM}: {st} {cls} {key}: LTRquest not the clear best ({v})")
                num("legend/D", st, "ltrquest", f"{CLS3_NAME[cls]} {key}, LTRquest vs next best ({NAME[top[1]]})",
                    f"{v['ltrquest']:.4f} vs {v[top[1]]:.4f}", src)
        for cls in ("partial", "proteinless"):
            if P[(st, "dante_ltr", cls)]["n_calls"] != 0:
                raise SystemExit(f"{STEM}: DANTE_LTR makes {CLS3_NAME[cls]} calls in {st}")
            num("legend/D", st, "dante_ltr", f"{CLS3_NAME[cls]} calls (summed over 3 genomes)", 0, src)
        num("legend/D", st, "ltrquest", "protein-less precision vs LTR_retriever, LTR_Finder, HiTE",
            f"{lq:.4f} vs {', '.join(f'{v:.4f}' for v in sp)}", src)
    # aging pointer (Fig. 1C-D, coding and protein-less aging series): LTRquest d50 coding vs protein-less
    d50 = tsv(NA / "div_d50.tsv").set_index(["tool", "cls"])
    for cls in CLS2:
        r = d50.loc[("ltrquest", cls)]
        num("legend/C (Fig. 1C-D pointer)", "div", "ltrquest", f"{CLS3_NAME[cls]} d50, K2P [{r.d50_lo}, {r.d50_hi}]",
            r.d50, "experiments/nonauto/div_d50.tsv")
    if not d50.loc[("ltrquest", "proteinless")].d50 < d50.loc[("ltrquest", "coding")].d50:
        raise SystemExit(f"{STEM}: legend 'LTRquest loses protein-less elements sooner' fails")
    ratios = []
    for ds in ("rice", "finch"):
        for t in ("ltrquest", "ltrharvest", "look4ltrs", "hirepeat", "ltrpred"):
            p = R["prec"][(ds, t)]
            ratios.append(float(p["supported_share"]) / float(p["null_supported_share"]))
    num("legend/E", "rice+finch", "permissive", "protein-less support / chance, range",
        f"{min(ratios):.2f}-{max(ratios):.2f}", "real_calls.tsv (supported_share / null_supported_share; LTRquest, "
        "LTRharvest, Look4LTRs, HiRepEat, LTRpred)")
    # headline: only LTRquest and Look4LTRs lose <= 0.03; every other caller loses more
    for _, r in g.iterrows():
        for k in ("gap", "length_part", "residual_part"):
            num("legend/headline", r["set"], r["tool"], f"protein-less minus coding recall, {k} [{r[k + '_lo']}, "
                f"{r[k + '_hi']}]", r[k], "sim_gap.tsv (family bootstrap CI)")
    gi = g.set_index(["set", "tool"])
    for t in TOOLS:
        lo = min(gi.loc[("S1", t), "gap"], gi.loc[("S2", t), "gap"])
        if (t in ("ltrquest", "look4ltrs")) != (lo >= -0.03):
            raise SystemExit(f"{STEM}: headline 'only LTRquest and Look4LTRs about as often' fails for {t} ({lo})")
    return min(ratios), max(ratios)


def legend_text(S, P, R, exs):
    rmin, rmax = extra_numbers(S, P, R)
    pr = R["prec"]
    c = R["comp"]
    if not (round(rmin) == 2 and round(rmax) == 5):
        raise SystemExit(f"{STEM}: legend says two to five times chance, table gives {rmin:.2f}-{rmax:.2f}")
    fl = S["floors"]
    if not (fl[("hirepeat", "proteinless")][0] == 0 and fl[("hirepeat", "coding")][0] == 0
            and all(fl[(t, k)][0] == 0 for t in ("ltrpred", "dante_ltr") for k in CLS2)):
        raise SystemExit(f"{STEM}: legend C size floors fail: {fl}")
    hite_n, hite_tot, _ = fl[("hite", "proteinless")]
    if not (0 < hite_n <= 0.02 * hite_tot and fl[("hite", "coding")][0] == 0):
        raise SystemExit(f"{STEM}: legend C 'HiTE almost none below 2 kb' fails: {fl[('hite', 'proteinless')]}")
    ex0 = exs[0]
    lqg, hrg = pr[("gorilla", "ltrquest")], pr[("gorilla", "hirepeat")]
    p2 = {t: P[("S2", t, "proteinless")]["precision"] for t in ("ltrquest", "ltr_retriever", "ltr_finder", "hite")}
    lo_s, hi_s = min(v for t, v in p2.items() if t != "ltrquest"), max(v for t, v in p2.items() if t != "ltrquest")
    txt = (
        "**Fig. S1. Size floors and protein rules remove LTR retrotransposons without proteins of their own from "
        "annotation.** Protein-less LTR retrotransposons are lost through rules written for elements that carry their "
        "own proteins, not because they are hard to see. Callers with a minimum element length miss every shorter "
        "copy and callers "
        "that need protein domains miss nearly all. Only LTRquest and Look4LTRs find them about as often as coding "
        "copies, the others losing some, partly at limits on element or LTR length. This analysis was added after "
        "the benchmark design was fixed. The simulated genomes (B to D) were made with PrinTE, which shares an author "
        "with LTRquest, and LTRquest was tuned on PrinTE simulations and on rice.\n\n"
        f"**A**, Both ends of a typical protein-less rice element, a {ex0['name']} copy without protein domains, with "
        "its target-site duplication (boxed) and TG…CA termini.\n\n"
        "**B**, Every protein-less (left) and coding (right) element of six simulated genomes (two regimes, three "
        "genomes each), one column each from shortest to longest, colored where a caller found it, with overall "
        "recall beside each row. **C**, Recall against element length (95% intervals from resampling element "
        "families). HiRepEat finds nothing shorter than 2 kb and HiTE "
        f"almost nothing ({hite_n} of {hite_tot:,} protein-less elements), LTRpred and DANTE_LTR nothing shorter than "
        "4 kb. At equal length, LTRpred and DANTE_LTR still find fewer protein-less than coding elements, and "
        "LTRharvest misses the longest protein-less ones. As elements age, LTRquest loses protein-less copies sooner "
        "than coding ones (Fig. 1C–D).\n\n"
        "**D**, Precision, recall and F1 for coding elements (at least two kinds of LTR retrotransposon domain from "
        "one superfamily), partly coding ones (some domains, short of that) and protein-less ones in the three "
        "LTR-rich, deeply nested genomes (S2, mean of three), recall counted by each element's family and precision "
        "by the domains inside each call. The coding and protein-less columns are those of Fig. 1B, and this panel "
        "adds the partly coding class. Like Fig. 1B, D counts every simulated element, including copies of library "
        "families without annotated LTRs, which nearly every caller misses, so its recall runs below B's. "
        "LTRquest's protein-less calls are less precise "
        f"({p2['ltrquest']:.2f}) than those of LTR_retriever, LTR_Finder and HiTE ({lo_s:.2f} to {hi_s:.2f}), so its "
        "lead in that class is recall. DANTE_LTR makes no partly coding or protein-less calls. –, no calls of that "
        "class (precision undefined). The moderate "
        "regime (S1) gives the same picture. Simulated genomes lack features of real genomes such as tandem "
        "duplicated genes and incomplete curated libraries, so precision in simulations runs higher than on real "
        "genomes.\n\n"
        "**E**, Share of each caller's protein-less calls in five real genomes that the curated library supports "
        "(95% interval), against chance and against the caller's coding calls. Strict callers make few protein-less "
        "calls. Permissive callers' protein-less calls, LTRquest's included, are supported at or below chance in "
        "*A. thaliana* and at two to five times chance in rice and zebra finch. In gorilla, where "
        f"{round(100 * c['gorilla']['proteinless'])}% of library-supported LTR retrotransposon loci are "
        f"protein-less, LTRquest reports {int(lqg['supported']):,} supported protein-less elements and HiRepEat "
        f"{int(hrg['supported']):,} among {int(hrg['n_calls']):,} calls. Maize precision is uninformative because "
        "chance support exceeds every caller's. Dark gray, the caller failed on that genome (it did not finish "
        "within 36 hours, exceeded 200 GB of memory or stopped without calls)."
    )
    num("legend/D", "S2", "", "protein-less precision, LTR_retriever / LTR_Finder / HiTE range (2 decimals)",
        f"{lo_s:.2f} to {hi_s:.2f}", D_SRC)
    num("legend/D", "S2", "ltrquest", "protein-less precision (2 decimals)", f"{p2['ltrquest']:.2f}",
        D_SRC)
    num("legend/E", "gorilla", "", "protein-less share of library-supported loci (%)",
        round(100 * c["gorilla"]["proteinless"]), "real_loci.tsv (n_loci)")
    num("legend/E", "gorilla", "ltrquest", "supported protein-less calls", int(lqg["supported"]), "real_calls.tsv")
    num("legend/E", "gorilla", "hirepeat", "supported protein-less calls of all protein-less calls",
        f"{int(hrg['supported'])} of {int(hrg['n_calls'])}", "real_calls.tsv")
    return txt


# ============================================================================ main
def build():
    load_missing()
    S = load_sim()
    log("simulations loaded")
    P = load_prf()
    log("scorecard table loaded and checked")
    rec_u = check_basis(S, P)
    log("legend D basis clause checked; S2 recall on copies without annotated LTRs: "
        + ", ".join(f"{NAME[t]} {v:.3f}" for t, v in rec_u.items()))
    R = load_real()
    log("real genomes loaded")
    exs = load_example()
    W = 180.0
    geo = dict(LAB=19.0, RH=6.7, rgap=1.3, hdr=13.0, B_bw=23.0, B_gap=1.4, B_txt=5.4, xC=86.5, wC=29.5,
               xD=122.6, cwD=6.0, ggD=0.9, cw=29.0, cg=3.2, RHh=3.2)
    fig = plt.figure(figsize=(W * MM, FIG_H * MM))
    yA = 3.0
    letter(fig, 0, yA - 1.0, "A")
    yAend, _ = panel_ends(fig, geo["LAB"], yA, exs[0])
    yB = yAend + 2.6
    letter(fig, 0, yB, "B")
    letter(fig, geo["xC"] - 10.4, yB, "C")
    letter(fig, geo["xD"] - 4.0, yB, "D")
    yBend = panel_tools(fig, yB, S, P, geo)
    yE = yBend + 1.0
    letter(fig, 0, yE, "E")
    yend = panel_E(fig, yE, R, geo)
    log(f"figure height used: {yend:.1f} mm of {FIG_H}")
    if yend > FIG_H + 0.5:
        raise SystemExit(f"{STEM}: content ({yend:.1f} mm) runs past the page ({FIG_H} mm)")
    return fig, S, P, R, exs


def main(argv=None):
    global VERBOSE
    ap = argparse.ArgumentParser(description="Fig. S1 figS1_nonauto: protein-less LTR-RT benchmark figure "
                                             "(figure + legend PDF, 300-dpi PNG, numbers table).")
    ap.add_argument("-v", "--verbose", action="store_true", help="progress and checks")
    a = ap.parse_args(argv)
    VERBOSE = a.verbose
    fig, S, P, R, exs = build()
    short = legend_text(S, P, R, exs)
    res = poafig.write_figure(STEM, fig, short, None, HERE, str(NA), str(Path(__file__).resolve()))
    log(f"wrote {res['pdf']} ({res['fig_mm']}, legend {res['short_mm']} mm, {res['short_words']} words)")
    with open(HERE / f"{STEM}_numbers.tsv", "w") as f:
        w = csv.DictWriter(f, fieldnames=["panel", "dataset", "tool", "quantity", "value", "source"], delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(NUM)
    log(f"numbers: {len(NUM)} rows")


if __name__ == "__main__":
    main()
