#!/usr/bin/env python3
"""Collect every finished evaluation output into report tables (TSV), one per measure; rerun as results arrive.

Tables (report/tables/):
  sim_elements.tsv  simulations, element level: mean and SD over seeds per set (S1, S2); RO >= 0.9 unless a column
                    says otherwise; std = PrinTE's intact set, ext = + rebuilt hosts
  sim_boundaries.tsv / sim_fp_anatomy.tsv  simulations: boundary error of matched calls; what false calls hit
  sim_crosscheck.tsv  our std F1 vs PrinTE's own evaluator (bedtools.py --printe-intact -r 0.9) on the same calls
  sim_masking.tsv   simulations, base level: each library's masking vs the exact LTR-RT truth (mean over seeds)
  masking.tsv       real genomes, EDTA lib-test LTR row vs the curated library; 95% block-bootstrap CI and the paired
                    difference to LTRquest when eval/lib/<ds>.ci.tsv exists (scripts/libtest_ci.py)
  proxies.tsv       real genomes, per-call sequence proxies
  phasing.tsv       allopolyploids, all copies; k80 = fewest copies (of 500/1k/2k/5k/10k) from which every larger
                    tested subset is correct in >= 80% of 5 seeds and the full set is correct ('<=500' if already at
                    the smallest size; 'all (n)' if no tested subset qualifies but all n copies do; 'never'; 'pending' if the
                    power curve has not run yet)
  phasing_power.tsv allopolyploids, random subsets: P(correct) and root-split strength per copy number
  resources.tsv     every tool run: status, wall and CPU hours, peak memory, load average at start
Column 'role': primary (the pre-registered row), secondary (native library / native phasing path), or post_hoc
(ltrquest_hc, developer diagnostic added after the primary results were seen).
usage: summarize.py [--outdir report/tables] [-v]
"""
import argparse
import csv
import glob
import json
import math
import os
import statistics as st
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import libtest_ci  # noqa: E402

B = os.environ["BENCH_ROOT"]
TOOLS = ["ltrquest", "ltr_retriever", "ltrharvest", "ltr_finder", "look4ltrs", "hirepeat", "hite", "ltrpred",
         "dante_ltr"]
POST_HOC = {"ltrquest_hc"}
SUBS = [500, 1000, 2000, 5000, 10000]
SIM = "S[12]_*"        # approach 1 sets S1, S2; the approach 1b divergence series (SD*) has its own tables (div_score.py)


def role(name):
    base = name.split(".")[0].split("/")[0]
    if base in POST_HOC:
        return "post_hoc"
    if name.endswith(".native") or base == "ltrquest_native":
        return "secondary"
    return "primary"


def rows(path):
    with open(path) as f:
        return [r for r in csv.DictReader(f, delimiter="\t") if r.get("label") != "label"]


def num(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return None if v != v else v


def ms(v, d=4):
    """mean and SD of the non-missing values, formatted; '.' if none."""
    v = [x for x in v if x is not None]
    if not v:
        return [".", "."]
    return [f"{st.mean(v):.{d}f}", f"{st.stdev(v):.{d}f}" if len(v) > 1 else "."]


def write(path, header, data):
    with open(path, "w") as f:
        f.write("\t".join(header) + "\n")
        for r in data:
            f.write("\t".join(str(x) for x in r) + "\n")


def tool_order(t):
    base = t.split(".")[0]
    return (TOOLS.index(base) if base in TOOLS else len(TOOLS), t)


def sim_elements(out):
    d = defaultdict(lambda: defaultdict(list))
    for p in glob.glob(f"{B}/eval/sim/{SIM}/*.metrics.tsv"):
        ds, tool = p.split("/")[-2], os.path.basename(p)[:-len(".metrics.tsv")]
        for r in rows(p):
            k = (r["scope"], r["ro"])
            for m in ("precision", "recall", "f1"):
                d[(ds[:2], tool)][(k, m)].append(num(r[m]))
            d[(ds[:2], tool)]["seeds"].append(ds)
    cols = [(("std", "0.9"), "f1"), (("std", "0.9"), "precision"), (("std", "0.9"), "recall"),
            (("ext", "0.9"), "f1"), (("ext", "0.9"), "precision"), (("ext", "0.9"), "recall"),
            (("ext:host", "0.9"), "recall"), (("ext:tenant", "0.9"), "recall"), (("ext:depth>=2", "0.9"), "recall"),
            (("nest_links", "0.9"), "f1"), (("std", "0.8"), "f1"), (("std", "0.95"), "f1")]
    names = ["f1_std", "prec_std", "rec_std", "f1_ext", "prec_ext", "rec_ext", "rec_host", "rec_tenant",
             "rec_depth2plus", "nest_link_f1", "f1_std_ro0.8", "f1_std_ro0.95"]
    header = ["set", "tool", "role", "n_seeds"] + [x for n in names for x in (n, n + "_sd")]
    data = []
    for (s, tool), v in sorted(d.items(), key=lambda kv: (kv[0][0], tool_order(kv[0][1]))):
        data.append([s, tool, role(tool), len(set(v["seeds"]))] + [x for c in cols for x in ms(v[c])])
    write(f"{out}/sim_elements.tsv", header, data)
    return data


def sim_masking(out):
    d = defaultdict(list)
    for p in glob.glob(f"{B}/eval/simlib/{SIM}/*.tsv"):
        ds, lib = p.split("/")[-2], os.path.basename(p)[:-4]
        if lib.endswith(".uniform"):
            continue
        for r in rows(p):
            d[(ds[:2], lib)].append((num(r["F1"]), num(r["sens"]), num(r["prec"]), num(r["masked_mb"])))
    header = ["set", "lib", "role", "n_seeds", "F1", "F1_sd", "sens", "sens_sd", "prec", "prec_sd", "masked_mb"]
    data = [[s, lib, role(lib), len(v)] + ms([x[0] for x in v]) + ms([x[1] for x in v]) + ms([x[2] for x in v])
            + ms([x[3] for x in v], 1)[:1] for (s, lib), v in sorted(d.items(), key=lambda kv: (kv[0][0], tool_order(kv[0][1])))]
    write(f"{out}/sim_masking.tsv", header, data)


def sim_diagnostics(out):
    """boundary accuracy and false-call anatomy (means over seeds), and the cross-check against PrinTE's evaluator."""
    bd, an, xc = defaultdict(list), defaultdict(list), []
    for p in glob.glob(f"{B}/eval/sim/{SIM}/*.boundaries.tsv"):
        ds, tool = p.split("/")[-2], os.path.basename(p)[:-len(".boundaries.tsv")]
        for r in rows(p):
            bd[(ds[:2], tool, r["end"])].append([num(r[k]) for k in ("exact", "le5", "le20", "median_abs")])
    write(f"{out}/sim_boundaries.tsv", ["set", "tool", "role", "end", "n_seeds", "frac_exact", "frac_le5bp",
                                        "frac_le20bp", "median_abs_bp"],
          [[s_, t, role(t), e, len(v)] + [ms([x[i] for x in v], 3)[0] for i in range(4)]
           for (s_, t, e), v in sorted(bd.items(), key=lambda kv: (kv[0][0], tool_order(kv[0][1]), kv[0][2]))])
    for p in glob.glob(f"{B}/eval/sim/{SIM}/*.anatomy.tsv"):
        ds, tool = p.split("/")[-2], os.path.basename(p)[:-len(".anatomy.tsv")]
        for r in rows(p):
            an[(ds[:2], tool, r["category"])].append(num(r["frac_of_calls"]))
    write(f"{out}/sim_fp_anatomy.tsv", ["set", "tool", "role", "category", "n_seeds", "frac_of_calls"],
          [[s_, t, role(t), c, len(v), ms(v)[0]]
           for (s_, t, c), v in sorted(an.items(), key=lambda kv: (kv[0][0], tool_order(kv[0][1]), kv[0][2]))])
    for p in sorted(glob.glob(f"{B}/eval/sim/{SIM}/*.printe.txt")):
        ds, tool = p.split("/")[-2], os.path.basename(p)[:-len(".printe.txt")]
        f1p = next((line.split()[2] for line in open(p) if line.startswith("F1")), None)
        m = f"{B}/eval/sim/{ds}/{tool}.metrics.tsv"
        f1m = next((r["f1"] for r in rows(m) if r["scope"] == "std" and r["ro"] == "0.9"), None) \
            if os.path.exists(m) else None
        if f1p and f1m:
            xc.append([ds, tool, f1m, f1p, f"{float(f1m) - float(f1p):+.4f}"])
    write(f"{out}/sim_crosscheck.tsv", ["ds", "tool", "f1_std_ours", "f1_printe_bedtools", "diff"], xc)


# masking runs left incomplete (not scored; no library rule added to make them feasible)
INCOMPLETE = {
    ("finch", "hirepeat.intact"): "incomplete: masking projected ~10 days (457/19693 batches in 5.9 h); not run",
    ("gorilla", "hirepeat.intact"): "incomplete: library build stopped at the 36 h cap (near-duplicate removal of "
                                    "43,620 calls, 395 Mb, by cd-hit-est); library not built, so not masked and no "
                                    "BM_RM2 (calls: 43,695, 827 Mb of span, longest 26.8 Mb)",
    ("gorilla", "ltr_finder.intact"): "incomplete: masking stopped at the 36 h cap at batch 26,019 of 61,148 (slow in "
                                      "tandem arrays); not scored",
    ("gorilla", "ltrquest.intact"): "incomplete: masking stopped at the 36 h cap at batch 55,136 of 61,148 (90%; 36 h "
                                    "counted from its 11-14 12:48 restart after the kernel stall); not scored",
}


def lib_mb(path):
    """total sequence length of a FASTA in Mb (None if absent)"""
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        return sum(len(x.strip()) for x in fh if not x.startswith(">")) / 1e6


# maize intact-call libraries are 122-1,492 Mb (<= 91 Mb elsewhere). Masking rates measured with MIMALLOC_ALLOW_THP=0
# at 16 threads, 37,882 batches; rate falls faster than library size grows
MAIZE_RATE = {"ltrpred": 1410, "hite": 593, "dante_ltr": 146}   # batches/h
for _t in TOOLS + ["ltrquest_hc"]:
    if os.path.exists(f"{B}/calls/maize/{_t}.tsv"):
        _mb = lib_mb(f"{B}/libs/maize/{_t}.intact.fa")
        if _t == "ltrpred":       # scored row replaces this once its targeted run finishes
            _why = "masking running (targeted run, 36 h cap; measured 1,410 batches/h -> ~27 h)"
        elif _t in MAIZE_RATE:
            _why = f"masking measured at {MAIZE_RATE[_t]} batches/h -> ~{37882 / MAIZE_RATE[_t]:.0f} h > 36 h; not run"
        else:
            _why = "library as large as or larger than dante_ltr's (445 Mb, ~260 h); masking > 36 h; not run"
        INCOMPLETE[("maize", f"{_t}.intact")] = (
            f"incomplete: {_why}" + (f"; library {_mb:.0f} Mb" if _mb is not None else "; library not built"))


def masking(out):
    header = ["ds", "lib", "role", "F1", "F1_lo", "F1_hi", "prec", "sens", "FP_mb", "FN_mb", "dF1_vs_ltrquest",
              "dF1_lo", "dF1_hi", "ci"]
    data = []
    for dsdir in sorted(glob.glob(f"{B}/eval/lib/*/")):
        ds = os.path.basename(os.path.normpath(dsdir))
        ci = {}
        scored = [p for p in glob.glob(f"{dsdir}*.tsv") if not p.endswith(".uniform.tsv") and os.path.getsize(p) > 0]
        cif = f"{B}/eval/lib/{ds}.ci.tsv"
        if scored and (not os.path.exists(cif) or max(map(os.path.getmtime, scored)) > os.path.getmtime(cif)):
            libtest_ci.main([ds])                         # refresh CIs when any library was (re)scored
        if os.path.exists(f"{B}/eval/lib/{ds}.ci.tsv"):
            ci = {r["lib"]: r for r in csv.DictReader(open(f"{B}/eval/lib/{ds}.ci.tsv"), delimiter="\t")}
        for p in sorted(glob.glob(f"{dsdir}*.tsv"), key=lambda p: tool_order(os.path.basename(p)[:-4])):
            lib = os.path.basename(p)[:-4]
            if lib.endswith(".uniform") or os.path.getsize(p) == 0:
                continue
            r = rows(p)[0]
            c = ci.get(lib)
            ok = c is not None and c["check"] == "ok"
            data.append([ds, lib, role(lib), f"{float(r['F1']):.4f}", c["F1_lo"] if ok else ".",
                         c["F1_hi"] if ok else ".", f"{float(r['prec']):.4f}", f"{float(r['sens']):.4f}",
                         f"{int(r['FP']) / 1e6:.2f}", f"{int(r['FN']) / 1e6:.2f}",
                         c["dF1"] if ok else ".", c["dF1_lo"] if ok else ".", c["dF1_hi"] if ok else ".",
                         f"block {int(c['block_bp']) // 1_000_000} Mb x {c['n_blocks']}" if ok else "none"])
    have = {(x[0], x[1]) for x in data}
    data += [[ds, lib, role(lib)] + ["."] * 10 + [why] for (ds, lib), why in INCOMPLETE.items() if (ds, lib) not in have]
    data.sort(key=lambda x: (x[0], tool_order(x[1])))
    write(f"{out}/masking.tsv", header, data)
    return data


def proxies(out):
    data = []
    for p in sorted(glob.glob(f"{B}/eval/proxy/*/*.tsv")):
        ds, tool = p.split("/")[-2], os.path.basename(p)[:-4]
        for r in rows(p):
            data.append([ds, tool, role(tool), r["n_calls"], r["median_len"], r["frac_tgca"], r["frac_tsd"],
                         r["frac_curated"]])
    data.sort(key=lambda x: (x[0], tool_order(x[1])))
    write(f"{out}/proxies.tsv", ["ds", "tool", "role", "n_calls", "median_len", "frac_tgca", "frac_tsd",
                                 "frac_curated"], data)


def copies_for_80(p_correct, full_ok, n):
    """Fewest tested copies k with P(correct) >= 0.8 at k and every larger tested size, given the full set is correct.
    p_correct: {copies: fraction of seeds correct}."""
    if not full_ok:
        return "never"
    ks = sorted(p_correct)
    ok = [k for k in ks if all(p_correct[j] >= 0.8 for j in ks if j >= k)]
    if ok:
        return f"<={min(ok)}" if min(ok) == SUBS[0] else str(min(ok))
    return f"all ({n})"                                 # no tested subset qualifies; only the full set phases


def phasing(out):
    power = defaultdict(lambda: defaultdict(list))
    for p in glob.glob(f"{B}/eval/phase/*/*/sub/*/score.tsv"):
        parts = p.split("/")
        ds, tool, (k, _) = parts[-5], parts[-4], parts[-2].split("_")
        for r in rows(p):
            power[(ds, tool)][int(k)].append((int(r["correct"] or 0), num(r["strength"])))
    pdata = []
    for (ds, tool), v in sorted(power.items(), key=lambda kv: (kv[0][0], tool_order(kv[0][1]))):
        for k in sorted(v):
            c = [x[0] for x in v[k]]
            pdata.append([ds, tool, role(tool), k, len(c), f"{sum(c) / len(c):.2f}"] + ms([x[1] for x in v[k]], 3))
    write(f"{out}/phasing_power.tsv", ["ds", "tool", "role", "copies", "n_seeds", "p_correct", "strength",
                                       "strength_sd"], pdata)
    data = []
    for p in sorted(glob.glob(f"{B}/eval/phase/*/*/score.tsv")):
        ds, tool = p.split("/")[-3], p.split("/")[-2]
        for r in rows(p):
            pc = {k: sum(x[0] for x in xs) / len(xs) for k, xs in power.get((ds, tool), {}).items()}
            # no subsets yet although n >= the smallest subset size: the power curve has not run
            k80 = "pending" if not pc and int(r["n_ltr"]) >= SUBS[0] else copies_for_80(pc, r["correct"] == "1", r["n_ltr"])
            data.append([ds, tool, role(tool), r["n_sg"], r["correct"], r["accuracy"], r["replicability"],
                         r["strength"], r["margin"], r["n_ltr"], k80])
    data.sort(key=lambda x: (x[0], tool_order(x[1])))
    write(f"{out}/phasing.tsv", ["ds", "tool", "role", "n_subgenomes", "correct", "accuracy", "replicability",
                                 "strength", "margin", "n_copies", "k80"], data)


def resources(out):
    data = []
    for p in sorted(glob.glob(f"{B}/runs/*/*/resources.json")):
        ds, tool = p.split("/")[-3], p.split("/")[-2]
        if ds == "smoke" or ds.startswith("SD"):   # SD*: approach 1b divergence series, not the main jobs
            continue
        j = json.load(open(p))
        # both CPU measures are lower bounds: rusage misses descendants never reaped by the waiter (a job killed at
        # the cap loses its main process), the 2 s /proc sampler misses short-lived processes; report the larger
        cpu = max(j["cpu_s"], j.get("cpu_proc_s") or 0)
        data.append([ds, tool, j["status"], j.get("exit_code", "."), f"{j['wall_s'] / 3600:.2f}",
                     f"{cpu / 3600:.2f}", f"{j['peak_rss_gb']:.1f}", f"{j.get('peak_pss_gb', 0):.1f}",
                     j["threads"], f"{j['load_start'][0]:.0f}",
                     # user vs kernel CPU (rusage; kernel-heavy runs suggest memory stalls);
                     # "." when rusage saw < half the CPU (runs killed at a cap lose their unreaped processes)
                     *([f"{j['cpu_user_s'] / 3600:.2f}", f"{j['cpu_sys_s'] / 3600:.2f}",
                        f"{j['cpu_sys_s'] / j['cpu_s']:.2f}"] if j["cpu_s"] >= 0.5 * cpu > 0 else [".", ".", "."])])
    data.sort(key=lambda x: (x[0], tool_order(x[1])))
    write(f"{out}/resources.tsv", ["ds", "tool", "status", "exit_code", "wall_h", "cpu_h", "peak_rss_gb",
                                   "peak_pss_gb", "threads", "load_at_start", "user_h", "sys_h", "sys_share"], data)
    return data


def bm_rm2(out):
    """BM_RM2 (RepeatModeler2's library benchmark): curated LTR families perfect/good/present."""
    data = []
    for p in sorted(glob.glob(f"{B}/eval/rm2/*/*.tsv")):
        ds, lib = p.split("/")[-2], os.path.basename(p)[:-4]
        for r in rows(p):
            if r["set"] != "LTR":
                continue
            n = int(r["n_curated"])
            data.append([ds, lib, role(lib), n, r["perfect"], r["good"], r["perfect_or_good"], r["present"],
                         f"{int(r['perfect_or_good']) / n:.4f}" if n else ".", f"{int(r['present']) / n:.4f}" if n else "."])
    data.sort(key=lambda x: (x[0], tool_order(x[1])))
    write(f"{out}/bm_rm2.tsv", ["ds", "lib", "role", "n_curated_LTR", "perfect", "good", "perfect_or_good", "present",
                                "frac_perfect_or_good", "frac_present"], data)


def approach4(out):
    """Approach 4 (post hoc; rules frozen before any per-tool number): element_library, unsupported_evidence, library_age_bias.
    Pools use the 9 primary tools present for a genome; complete only when every tool has finished."""
    import adj_summary as su

    def f(x, d=4):
        return "." if x is None or (isinstance(x, float) and math.isnan(x)) else f"{x:.{d}f}"

    def ci(x):
        return f"{f(x[0])}-{f(x[1])}"

    el, un, ab = [], [], []
    for ds in ["arab", "rice", "finch", "maize", "gorilla"]:
        rb = su.load_ds(f"{B}/eval/adj/{ds}", TOOLS + sorted(POST_HOC))
        if not rb:
            continue
        pool = [t for t in TOOLS if t in rb]
        res = su.element_library(rb, pool)
        for t in sorted(rb, key=tool_order):
            d = res[t]
            el.append([ds, t, role(t), d["n_calls"], len(pool), d["n_loci_lib"], f(d["precision_lib"]),
                       ci(d["precision_lib_ci"]), f(d["recall_lib"]), ci(d["recall_lib_ci"]), f(d["recall_lib_ro50"]),
                       f(d["F1_lib"]), ci(d["F1_lib_ci"]), d["n_loci_ext"], f(d["recall_ext"]),
                       f(d["F1_ext"]), ci(d["F1_ext_ci"])])
            u = su.unsupported(rb[t], t)
            un.append([ds, t, role(t), u["n_calls"], u["n_supported"], u["n_unsupported"]] +
                      [f(u[f"{g}_{k}"], 3) for g in ("sup", "unsup", "null") for k in ("coding", "junk")] +
                      [f(u[f"unsup_{k}"], 3) for k in ("gene", "nonltr", "sat", "trf")] +
                      [u["unsup_coding_only"], u["unsup_junk_only"], u["unsup_coding_and_junk"], u["unsup_neither"],
                       f(u["sup_k2p_median"]), f(u["unsup_k2p_median"]), f(u["f_real"], 3), ci(u["f_real_ci"]),
                       u["f_real_bins"] or ".", f(u["precision_lib"]), f(u["precision_upper"]), u["flag"]])
        bins, n_loci, slope, sci = su.age_bias(rb, pool)
        for name, n, share, lo, hi in bins:
            ab.append([ds, len(pool), n_loci, name, n, f(share, 3), f"{f(lo, 3)}-{f(hi, 3)}", f(slope, 3), ci(sci)])
    write(f"{out}/element_library.tsv",
          ["ds", "tool", "role", "n_calls", "n_pool_tools", "n_loci_lib", "precision_lib", "precision_lib_ci",
           "recall_lib", "recall_lib_ci", "recall_lib_ro50", "F1_lib", "F1_lib_ci", "n_loci_ext", "recall_ext",
           "F1_ext", "F1_ext_ci"], el)
    write(f"{out}/unsupported_evidence.tsv",
          ["ds", "tool", "role", "n_calls", "n_supported", "n_unsupported", "sup_coding", "sup_junk", "unsup_coding",
           "unsup_junk", "null_coding", "null_junk", "unsup_gene", "unsup_nonltr", "unsup_sat", "unsup_trf",
           "unsup_coding_only", "unsup_junk_only", "unsup_coding_and_junk", "unsup_neither", "sup_k2p_median",
           "unsup_k2p_median", "f_real", "f_real_ci", "f_real_bins", "precision_lib", "precision_upper", "flag"], un)
    write(f"{out}/library_age_bias.tsv",
          ["ds", "n_pool_tools", "n_coding_loci", "age_bin", "n_loci", "supported_share", "wilson95", "logit_slope",
           "logit_slope_ci"], ab)


def main(argv=None):
    ap = argparse.ArgumentParser(description="Collect finished evaluation outputs into report/tables/*.tsv.")
    ap.add_argument("--outdir", default=f"{B}/report/tables", help="output directory (default report/tables)")
    ap.add_argument("-v", "--verbose", action="store_true", help="print row counts per table")
    a = ap.parse_args(argv)
    os.makedirs(a.outdir, exist_ok=True)
    sim_elements(a.outdir)
    sim_masking(a.outdir)
    sim_diagnostics(a.outdir)
    masking(a.outdir)
    proxies(a.outdir)
    phasing(a.outdir)
    resources(a.outdir)
    bm_rm2(a.outdir)
    approach4(a.outdir)
    if a.verbose:
        for p in sorted(glob.glob(f"{a.outdir}/*.tsv")):
            print(f"{os.path.basename(p)}\t{sum(1 for _ in open(p)) - 1} rows")
    print(f"summarize: tables in {a.outdir}")


if __name__ == "__main__":
    main()
