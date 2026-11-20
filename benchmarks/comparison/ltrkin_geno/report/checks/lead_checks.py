#!/usr/bin/env python3
"""Post hoc checks on the headline: does the leading tool's ltrkin yield survive the pool's known confounds?

Written after the per-tool numbers were seen. It changes no scoring rule: it reuses score.py's per-locus outcomes
(end-to-end ltrkin runs, the same matching, the same block bootstrap and seed) and only restricts the loci.

Subsets of pooled loci (every table below is given per subset where it has a `subset` column):
  all          every pooled locus;
  true_ends    loci pooled through their true element (basis true_confirmed): membership depends on no call's ends;
  not_lq_only  loci that some tool other than LTR_peel nominated (drops the LTR_peel-only loci).

  subsets.tsv  per subset, metric and tool: n (loci with a hit), value = n / n_loci with its 95% CI, and the paired
               difference d = value(tool) - value(LTR_peel) with its CI (same draws). Metrics:
                 yield         genotyped F/E;
                 recall        found (matched at RO >= 0.9);
                 found_not_fe  found but not genotyped F/E. Its d is the part of LTR_peel's yield lead over the tool
                               that arises at genotyping (yield lead - recall lead); `share` = that part / yield lead.
  gap.tsv      per subset and tool: loci where LTR_peel's matched call is genotyped F/E and the tool's is not, split
               by the tool's outcome (missed, partial, found with another genotype), and the reverse count.
  ends.tsv     all loci, per tool: found loci by the larger end error of the matched call (score.end_errors): exact
               (both ends in the microhomology window), le5, le20, gt20 bp, or no_true_ends; n_found and n_fe; a last
               row per dataset, tool "all", sums the primary tools (a locus counts once per tool).
  length.tsv   per subset, length bin (as by_length.tsv) and tool: n_loci and n_found.
  end_distance.tsv
               added later still (after Fig. 1 showed that the >20 bp bin of ends.tsv mixes ends 21 bp
               and 1,000 bp off). All loci, primary tools that finished: found loci with true ends (no_true_ends
               left out) by the larger |end error| of the matched call, in Fig. 1 E/F's half-decade bins: exact (both
               ends in the microhomology window, as score.end_errors), 1-3, 4-10, 11-31, 32-100, 101-316, >=317 bp.
               Columns n_found, n_fe, frac_fe and its Wilson 95% interval lo hi (as Fig. 1); median_far_bp, the
               median larger |end error| of the tool's found loci in ends.tsv's gt20 bin (statistics.median, repeated
               on each of the tool's rows); frac_calls_tsd, the share of the tool's calls (calls TSV, read_calls) that
               carry a TSD on that dataset. group: no_tsd = tools none of whose calls carries a TSD on any dataset
               checked where they finished, else tsd. Group rows (tool "all", group tsd / no_tsd / all) sum their
               tools' rows (a locus counts once per tool); their median_far_bp pools the group's far ends.
  end_distance_contrast.tsv
               per dataset and range of the larger |end error|: exact-31 (bins exact to 11-31), 32-100, 101-316, the
               pooled 32-316 and >=317 bp. F/E / found for the tsd and no_tsd groups with Wilson intervals, and the
               two-sided Fisher exact p of the 2x2 table (group x F/E or not). The pooled 32-316 row mixes distances
               (the per-bin rows compare at closer distance), so it also carries mh_p: the Mantel-Haenszel test
               stratified over the 32-100 and 101-316 bp bins (continuity-corrected chi-square, 1 df). Calls of
               different tools at one locus are not independent, so these p values are descriptive.
  end_distance_paired.tsv
               added in the same pass (post hoc): the within-locus test. Same insertion, same reads: does
               a far-off call fail where a near call succeeds? Per dataset and far range (32-100, 101-316, >=317 and
               all >=32 bp): loci with a near call (found, true ends, larger |end error| <= 31 bp, bins exact to
               11-31) and a far call in the range, from any two primary tools. Locus outcome: near = majority of the
               near calls at the locus genotyped F/E, a tie counting as F/E; far likewise over its far calls
               (n_tie_loci counts loci where either side tied). Columns: n_loci; the four cells near_ok_far_ok,
               near_ok_far_wrong, near_wrong_far_ok, near_wrong_far_wrong; the far calls by group (far_tsd_n,
               far_tsd_fe, far_no_tsd_n, far_no_tsd_fe; calls, not loci); sign_p, the exact two-sided sign test on
               the discordant loci (near_ok_far_wrong vs near_wrong_far_ok).

Self-check: before writing anything, the script dies with a message naming the mismatch unless its `all` rows
reproduce score.py's tables in <root>/report/tables: yield and CIs (ltrkin_yield.tsv, e2e rows), recall and CIs
(recall.tsv) and found per length bin (by_length.tsv); ends.tsv sums per tool must equal found and the e2e F/E count,
and gap.tsv's net count must equal the F/E difference. For end_distance.tsv, per tool the bins must sum to ends.tsv's
exact + le5 + le20 + gt20 (n_found and n_fe) with the exact bins equal, each group row must equal the sum of its
tools, and each contrast row must equal its tools' rows summed over the range's bins. For end_distance_paired.tsv,
the four cells must sum to n_loci, each locus must have a far call, the >=32 row's far calls must equal the sum over
the three far bins (each call is in one bin), its n_loci must lie between the largest bin's and their sum, and the far
calls must be a subset of end_distance.tsv's group counts in those bins. Then, not fatal, the end_distance group rows,
the 32-316 bp contrast and the TSD shares are compared with Fig. 1's numbers table (stale until the figure is rebuilt
after a roll-in).
"""
import argparse
import csv
import math
import os
import re
import statistics
import sys

sys.dont_write_bytecode = True
LG = os.path.join(os.environ["BENCH_ROOT"], "ltrkin_geno")   # data root
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "scripts"))
import score  # noqa: E402
from lg_config import ALL_TOOLS, BLOCK, DATASETS, SECONDARY  # noqa: E402
from lg_io import die, final_calls_path, read_calls  # noqa: E402

REF = "ltr_peel"
SUBSETS = ["all", "true_ends", "not_lq_only"]
METRICS = {"yield": lambda st, cls: cls == "F/E",
           "recall": lambda st, cls: st == "found",
           "found_not_fe": lambda st, cls: st == "found" and cls != "F/E"}


def nominators(locus):
    """Tools with a member (nominating call) at a pooled locus; members are 'tool:element_id:start:end,...'."""
    return {m.split(":")[0] for m in locus["members"].split(",") if m}


def subset_index(pool):
    """Subset name -> indices of its loci, in pool order."""
    keep = dict(all=[True] * len(pool),
                true_ends=[p["basis"] == "true_confirmed" for p in pool],
                not_lq_only=[nominators(p) != {REF} for p in pool])
    return {name: [i for i, x in enumerate(keep[name]) if x] for name in SUBSETS}


def gap_rows(ds, subs, outcomes, ok):
    """Per subset and tool: loci where REF is F/E and the tool is not (by the tool's class), and the reverse."""
    rows = []
    for name, idx in subs.items():
        ref = [outcomes[REF][i] for i in idx]
        for t in ok:
            if t == REF:
                continue
            cls = [outcomes[t][i] for i in idx]
            lost = [c for r, c in zip(ref, cls) if r == "F/E" and c != "F/E"]
            rows.append(dict(ds=ds, subset=name, tool=t, ref_fe_tool_not=len(lost), tool_missed=lost.count("missed"),
                             tool_partial=lost.count("partial"),
                             tool_found_other_genotype=sum(c not in ("missed", "partial") for c in lost),
                             tool_fe_ref_not=sum(c == "F/E" and r != "F/E" for r, c in zip(ref, cls))))
    return rows


def subset_rows(ds, pool, subs, status, outcomes, ok, blocks):
    """Per subset and metric: value, CI and paired difference vs REF (score.block_bootstrap, same seed and blocks)."""
    rows = []
    for name, idx in subs.items():
        if not idx:
            die(f"{ds}: subset {name} is empty")
        ub = [(pool[i]["chrom"], pool[i]["match_start"] // BLOCK) for i in idx]
        val = {}
        for metric, hit in METRICS.items():
            h = {t: [int(hit(status[t][i], outcomes[t][i])) for i in idx] for t in ok}
            ci = score.block_bootstrap([(ub[j], {t: (h[t][j], 1) for t in ok}) for j in range(len(idx))], blocks, ok)
            val[metric] = {t: sum(h[t]) / len(idx) for t in ok}
            for t in ok:
                rows.append(dict(ds=ds, subset=name, metric=metric, tool=t, n_loci=len(idx), n=sum(h[t]),
                                 value=val[metric][t], lo=ci[t]["lo"], hi=ci[t]["hi"],
                                 d=val[metric][t] - val[metric][REF], dlo=ci[t]["dlo"], dhi=ci[t]["dhi"]))
        for r in rows:                      # share of the yield lead that arises at genotyping
            if r["ds"] == ds and r["subset"] == name and r["metric"] == "found_not_fe" and r["tool"] != REF:
                lead = val["yield"][REF] - val["yield"][r["tool"]]
                r["share"] = r["d"] / lead if lead else None
    return rows


END_BINS = ["exact", "le5", "le20", "gt20", "no_true_ends"]


def end_bin(err):
    """Bin of a found locus by its matched call's larger absolute end error (None: the locus has no true ends)."""
    if err is None:
        return "no_true_ends"
    if err["exact"]:
        return "exact"
    worst = max(abs(err["left"]), abs(err["right"]))
    return "le5" if worst <= 5 else "le20" if worst <= 20 else "gt20"


def ends_rows(ds, res, outcomes, ok):
    """Per tool and end bin: found loci and those genotyped F/E; then the sum over primary tools."""
    rows, tot = [], {b: [0, 0] for b in END_BINS}
    for t in ok:
        cnt = {b: [0, 0] for b in END_BINS}
        for (st, _), err, cls in zip(res[t]["status_rows"], res[t]["errs"], outcomes[t]):
            if st != "found":
                continue
            b = end_bin(err)
            cnt[b][0] += 1
            cnt[b][1] += cls == "F/E"
        for b in END_BINS:
            rows.append(dict(ds=ds, tool=t, end_bin=b, n_found=cnt[b][0], n_fe=cnt[b][1],
                             frac_fe=cnt[b][1] / cnt[b][0] if cnt[b][0] else None))
            if t not in SECONDARY:
                tot[b][0] += cnt[b][0]
                tot[b][1] += cnt[b][1]
    rows += [dict(ds=ds, tool="all", end_bin=b, n_found=tot[b][0], n_fe=tot[b][1],
                  frac_fe=tot[b][1] / tot[b][0] if tot[b][0] else None) for b in END_BINS]
    return rows


# Fig. 1 E/F's half-decade bins on the larger |end error| (figures/fig1_ltrkin_yield/fig1_ltrkin_yield.py EDGES:
# edges at 10^(k/2) bp, rounded), after the exact bin.
DIST_EDGES = [1, 4, 11, 32, 101, 317]
DIST_BINS = ["exact", "1-3", "4-10", "11-31", "32-100", "101-316", ">=317"]
CONTRAST = {"exact-31": DIST_BINS[:4], "32-100": ("32-100",), "101-316": ("101-316",),   # end_distance_contrast.tsv
            "32-316": ("32-100", "101-316"), ">=317": (">=317",)}                       # ranges: their bins
GROUPS = ("tsd", "no_tsd")
NEAR_BINS = DIST_BINS[:4]                    # end_distance_paired.tsv: a near call is <= 31 bp off
PAIRED = {"32-100": ("32-100",), "101-316": ("101-316",), ">=317": (">=317",), ">=32": tuple(DIST_BINS[4:])}


def dist_bin(err):
    """Half-decade bin of a found locus by its matched call's larger |end error| (None: no true ends)."""
    if err is None:
        return None
    if err["exact"]:
        return "exact"
    return DIST_BINS[sum(max(abs(err["left"]), abs(err["right"])) >= e for e in DIST_EDGES)]


def wilson(k, n, z=1.959964):
    """Wilson score 95% interval of k / n, the formula of Fig. 1 (fig1_ltrkin_yield.wilson); (None, None) if n is 0."""
    if n == 0:
        return None, None
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - h), min(1.0, c + h)


def fisher_two_sided(a, b, c, d):
    """Two-sided Fisher exact p of the 2x2 table [[a, b], [c, d]]: the summed probability, under fixed margins, of
    every table no more likely than the observed one (exact integer weights, so ties are exact). None if empty."""
    r1, c1, n = a + b, a + c, a + b + c + d
    if n == 0:
        return None
    w = [math.comb(c1, x) * math.comb(n - c1, r1 - x) for x in range(max(0, r1 + c1 - n), min(r1, c1) + 1)]
    obs = math.comb(c1, a) * math.comb(n - c1, b)
    return sum(x for x in w if x <= obs) / math.comb(n, r1)


def sign_test(b, c):
    """Exact two-sided sign test of b against c discordant pairs (binomial, p = 1/2). None if b + c is 0."""
    n = b + c
    if n == 0:
        return None
    return min(1.0, 2 * sum(math.comb(n, i) for i in range(min(b, c) + 1)) / 2 ** n)


def mantel_haenszel(tables):
    """Mantel-Haenszel test of a common odds ratio of 1 over strata [[a, b], [c, d]]: continuity-corrected chi-square
    (|sum a - sum E[a]| - 1/2)^2 / sum Var[a], 1 df, p = erfc(sqrt(chi2 / 2)). Strata with n < 2 add nothing; None if
    no stratum has variance."""
    sa = se = sv = 0.0
    for a, b, c, d in tables:
        n = a + b + c + d
        if n < 2:
            continue
        sa += a
        se += (a + b) * (a + c) / n
        sv += (a + b) * (c + d) * (a + c) * (b + d) / (n * n * (n - 1))
    if sv == 0:
        return None
    chi2 = max(0.0, abs(sa - se) - 0.5) ** 2 / sv
    return math.erfc(math.sqrt(chi2 / 2))


def distance_counts(res, outcomes, ok):
    """Per primary tool: {bin: [n_found, n_fe]} over found loci with true ends, and the larger |end error| of its
    found loci in ends.tsv's gt20 bin."""
    out = {}
    for t in ok:
        if t in SECONDARY:
            continue
        cnt, far = {b: [0, 0] for b in DIST_BINS}, []
        for (st, _), err, cls in zip(res[t]["status_rows"], res[t]["errs"], outcomes[t]):
            b = dist_bin(err) if st == "found" else None
            if b is None:
                continue
            cnt[b][0] += 1
            cnt[b][1] += cls == "F/E"
            if end_bin(err) == "gt20":
                far.append(max(abs(err["left"]), abs(err["right"])))
        out[t] = (cnt, far)
    return out


def tsd_calls(root, ds, ok):
    """Per primary tool: (calls carrying a TSD, calls) in its final calls TSV (read_calls keeps a TSD of 3-8 ACGT)."""
    out = {}
    for t in ok:
        if t not in SECONDARY:
            calls = read_calls(final_calls_path(f"{root}/calls/{ds}", t), t)
            out[t] = (sum(c["tsd"] is not None for c in calls), len(calls))
    return out


def tsd_groups(tsd):
    """tool -> no_tsd if none of its calls carries a TSD on any dataset where it finished, else tsd.
    tsd: {ds: {tool: (n_tsd, n_calls)}}."""
    n = {}
    for per_ds in tsd.values():
        for t, (k, _) in per_ds.items():
            n[t] = n.get(t, 0) + k
    return {t: "no_tsd" if k == 0 else "tsd" for t, k in n.items()}


def distance_rows(ds, counts, group, tsd):
    """end_distance.tsv rows of one dataset: per tool and bin, then group rows (tool "all") tsd, no_tsd and all."""
    def row(tool, g, b, n, k, far, share):
        lo, hi = wilson(k, n)
        return dict(ds=ds, tool=tool, group=g, end_bin=b, n_found=n, n_fe=k, frac_fe=k / n if n else None, lo=lo,
                    hi=hi, median_far_bp=statistics.median(far) if far else None, frac_calls_tsd=share)
    rows = []
    for t, (cnt, far) in counts.items():
        share = tsd[t][0] / tsd[t][1] if tsd[t][1] else None
        rows += [row(t, group[t], b, cnt[b][0], cnt[b][1], far, share) for b in DIST_BINS]
    for g in GROUPS + ("all",):
        mem = [t for t in counts if g in ("all", group[t])]
        far = [x for t in mem for x in counts[t][1]]
        rows += [row("all", g, b, sum(counts[t][0][b][0] for t in mem), sum(counts[t][0][b][1] for t in mem), far,
                     None) for b in DIST_BINS]
    return rows


def range_sum(dist, ds, g, rng):
    """(n_found, n_fe) of group g's rows (tool "all") over the bins of contrast range rng."""
    sel = [x for x in dist if x["ds"] == ds and x["tool"] == "all" and x["group"] == g
           and x["end_bin"] in CONTRAST[rng]]
    return sum(x["n_found"] for x in sel), sum(x["n_fe"] for x in sel)


def contrast_rows(dist):
    """end_distance_contrast.tsv: per dataset and range (CONTRAST), tsd vs no_tsd F/E / found and Fisher's p."""
    rows = []
    for ds, rng in ((d, x) for d in dict.fromkeys(r["ds"] for r in dist) for x in CONTRAST):
        r = dict(ds=ds, range=rng)
        for g in GROUPS:
            n, k = range_sum(dist, ds, g, rng)
            lo, hi = wilson(k, n)
            r.update({f"{g}_n_found": n, f"{g}_n_fe": k, f"{g}_frac_fe": k / n if n else None, f"{g}_lo": lo,
                      f"{g}_hi": hi})
        a, c = r["tsd_n_fe"], r["no_tsd_n_fe"]
        r["fisher_p"] = fisher_two_sided(a, r["tsd_n_found"] - a, c, r["no_tsd_n_found"] - c)
        if rng == "32-316":                 # stratified over its two bins
            st = []
            for b in CONTRAST[rng]:
                (n1, k1), (n2, k2) = (range_sum(dist, ds, g, b) for g in GROUPS)
                st.append((k1, n1 - k1, k2, n2 - k2))
            r["mh_p"] = mantel_haenszel(st)
        rows.append(r)
    return rows


def locus_calls(res, outcomes, ok):
    """Per pooled locus: [(tool, distance bin, genotyped F/E)] of the primary tools' found calls with true ends."""
    prim = [t for t in ok if t not in SECONDARY]
    out = [[] for _ in outcomes[prim[0]]] if prim else []
    for t in prim:
        for i, ((st, _), err, cls) in enumerate(zip(res[t]["status_rows"], res[t]["errs"], outcomes[t])):
            b = dist_bin(err) if st == "found" else None
            if b is not None:
                out[i].append((t, b, cls == "F/E"))
    return out


def paired_rows(ds, loci, group):
    """end_distance_paired.tsv rows of one dataset (see the module docstring)."""
    def ok(calls):
        return 2 * sum(c[2] for c in calls) >= len(calls), 2 * sum(c[2] for c in calls) == len(calls)
    rows = []
    for rng, bins in PAIRED.items():
        r = dict(ds=ds, far_range=rng, n_loci=0, n_tie_loci=0, near_ok_far_ok=0, near_ok_far_wrong=0,
                 near_wrong_far_ok=0, near_wrong_far_wrong=0, far_tsd_n=0, far_tsd_fe=0, far_no_tsd_n=0,
                 far_no_tsd_fe=0)
        for calls in loci:
            near = [c for c in calls if c[1] in NEAR_BINS]
            far = [c for c in calls if c[1] in bins]
            if not near or not far:
                continue
            (a, ta), (b, tb) = ok(near), ok(far)
            r["n_loci"] += 1
            r["n_tie_loci"] += ta or tb
            r[f"near_{'ok' if a else 'wrong'}_far_{'ok' if b else 'wrong'}"] += 1
            for t, _, fe in far:
                r[f"far_{group[t]}_n"] += 1
                r[f"far_{group[t]}_fe"] += fe
        r["sign_p"] = sign_test(r["near_ok_far_wrong"], r["near_wrong_far_ok"])
        rows.append(r)
    return rows


PAIRED_COLS = ["ds", "far_range", "n_loci", "n_tie_loci", "near_ok_far_ok", "near_ok_far_wrong", "near_wrong_far_ok",
               "near_wrong_far_wrong", "far_tsd_n", "far_tsd_fe", "far_no_tsd_n", "far_no_tsd_fe", "sign_p"]
CELLS = ["near_ok_far_ok", "near_ok_far_wrong", "near_wrong_far_ok", "near_wrong_far_wrong"]
FAR = ["far_tsd_n", "far_tsd_fe", "far_no_tsd_n", "far_no_tsd_fe"]


def paired_check(dist, paired):
    """Self-check of end_distance_paired.tsv: list of mismatch messages (see the module docstring)."""
    bad = []
    for d in dict.fromkeys(r["ds"] for r in paired):
        p = {r["far_range"]: r for r in paired if r["ds"] == d}
        for rng, r in p.items():
            if sum(r[c] for c in CELLS) != r["n_loci"] or r["far_tsd_n"] + r["far_no_tsd_n"] < r["n_loci"]:
                bad.append(f"end_distance_paired.tsv {d} {rng}: cells do not sum to n_loci or a locus lacks a far "
                           f"call")
            for g in GROUPS:
                n, k = (sum(x[c] for x in dist if x["ds"] == d and x["tool"] == "all" and x["group"] == g
                            and x["end_bin"] in PAIRED[rng]) for c in ("n_found", "n_fe"))
                if r[f"far_{g}_n"] > n or r[f"far_{g}_fe"] > k:
                    bad.append(f"end_distance_paired.tsv {d} {rng} {g}: more far calls than end_distance.tsv holds")
        bins = [p[b] for b in PAIRED if b != ">=32"]
        if any(p[">=32"][c] != sum(x[c] for x in bins) for c in FAR) or not (
                max(x["n_loci"] for x in bins) <= p[">=32"]["n_loci"] <= sum(x["n_loci"] for x in bins)):
            bad.append(f"end_distance_paired.tsv {d} >=32: far calls or loci disagree with the three far bins")
    return bad


CONTRAST_COLS = ["ds", "range"] + [f"{g}_{c}" for g in GROUPS for c in ("n_found", "n_fe", "frac_fe", "lo", "hi")] \
    + ["fisher_p", "mh_p"]
DIST_COLS = ["ds", "tool", "group", "end_bin", "n_found", "n_fe", "frac_fe", "lo", "hi", "median_far_bp",
             "frac_calls_tsd"]


def distance_check(ends, dist, contrast):
    """Self-check of end_distance(_contrast): list of mismatch messages (see the module docstring)."""
    bad = []
    ekey = {(r["ds"], r["tool"], r["end_bin"]): r for r in ends}
    tools = {(r["ds"], r["tool"]) for r in dist if r["tool"] != "all"}
    for d, t in sorted(tools):
        mine = [r for r in dist if r["ds"] == d and r["tool"] == t]
        for k in ("n_found", "n_fe"):
            want = sum(ekey[(d, t, b)][k] for b in ("exact", "le5", "le20", "gt20"))
            if sum(r[k] for r in mine) != want or [r[k] for r in mine if r["end_bin"] == "exact"] != [
                    ekey[(d, t, "exact")][k]]:
                bad.append(f"end_distance.tsv {d} {t} {k}: bins do not sum to ends.tsv's loci with true ends")
    for r in dist:
        if r["tool"] != "all":
            continue
        mem = [x for x in dist if x["ds"] == r["ds"] and x["tool"] != "all" and x["end_bin"] == r["end_bin"]
               and r["group"] in ("all", x["group"])]
        if (r["n_found"], r["n_fe"]) != (sum(x["n_found"] for x in mem), sum(x["n_fe"] for x in mem)):
            bad.append(f"end_distance.tsv {r['ds']} group {r['group']} {r['end_bin']}: not the sum of its tools")
    for c in contrast:                      # against the tool rows, not the group rows contrast_rows sums
        for g in GROUPS:
            sel = [x for x in dist if x["ds"] == c["ds"] and x["tool"] != "all" and x["group"] == g
                   and x["end_bin"] in CONTRAST[c["range"]]]
            if (c[f"{g}_n_found"], c[f"{g}_n_fe"]) != (sum(x["n_found"] for x in sel), sum(x["n_fe"] for x in sel)):
                bad.append(f"end_distance_contrast.tsv {c['ds']} {c['range']} {g}: not the sum of its tools' bins")
    return bad


FIG1_ITEM = re.compile(r"genotyped correctly, end error (\S+) bp \((no TSDs|TSDs)\)$")
FIG1_VALUE = re.compile(r"(\d+)/(\d+) = [\d.]+% \(95% CI ([\d.]+)-([\d.]+)\)$")


def fig1_compare(path, dist, contrast):
    """Compare with Fig. 1's numbers table where they overlap (not fatal): the group bins of panels E/F (k, n and the
    interval to 0.1 points), the legend's 32-316 bp counts and each tool's share of calls with a TSD (4 decimals).
    Returns (n compared, list of differences), or None when the table is absent."""
    if not os.path.isfile(path):
        return None
    grp = {(r["ds"], r["group"], r["end_bin"]): r for r in dist if r["tool"] == "all"}
    share = {(r["ds"], r["tool"]): r["frac_calls_tsd"] for r in dist if r["tool"] != "all"}
    con = {(c["ds"], c["range"]): c for c in contrast}
    n, diff = 0, []
    for f in read_table(path):
        ds, item, val = f["dataset"], f["item"], f["value"]
        m = FIG1_ITEM.match(item)
        if f["panel"] in ("E", "F") and m:
            r = grp.get((ds, "no_tsd" if m.group(2) == "no TSDs" else "tsd", m.group(1)))
            if r is None:
                continue
            v = FIG1_VALUE.match(val)
            ci = tuple(None if x is None else f"{100 * x:.1f}" for x in (r["lo"], r["hi"]))
            mine = (r["n_fe"], r["n_found"]) + ci
            theirs = v and (int(v.group(1)), int(v.group(2)), v.group(3), v.group(4))
        elif f["panel"] in ("E", "F") and item == "share of calls with a TSD":
            if (ds, f["caller"]) not in share:
                continue
            mine, theirs = share[(ds, f["caller"])], float(val)
            mine = None if mine is None else round(mine, 4)
        elif f["panel"] == "legend" and item.startswith("genotyped correctly, end error 32-316 bp"):
            if (ds, "32-316") not in con:
                continue
            c = con[(ds, "32-316")]
            mine = f"{c['tsd_n_fe']}/{c['tsd_n_found']} | {c['no_tsd_n_fe']}/{c['no_tsd_n_found']}"
            theirs = val
        else:
            continue
        n += 1
        if mine != theirs:
            diff.append(f"{f['panel']} {ds} {f['caller'] or '-'} {item}: Fig. 1 {theirs}, end_distance {mine}")
    return n, diff


def length_rows(ds, pool, subs, status, ok):
    """Per subset, length bin (score.bin_of on the match span, as by_length.tsv) and tool: n_loci, n_found."""
    rows = []
    for name, idx in subs.items():
        for j, lb in enumerate(score.LEN_LABELS):
            ix = [i for i in idx if score.bin_of(pool[i]["match_end"] - pool[i]["match_start"]) == j]
            rows += [dict(ds=ds, subset=name, len_bin=lb, tool=t, n_loci=len(ix),
                          n_found=sum(status[t][i] == "found" for i in ix)) for t in ok]
    return rows


def read_table(path):
    if not os.path.isfile(path):
        die(f"self-check: table not found: {path}")
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def self_check(root, dss, subs_rows, gaps, ends, lengths, dist, contrast, paired):
    """Die unless the `all` rows reproduce score.py's tables and end_distance agrees with ends (see the module
    docstring); return a summary line."""
    tab = f"{root}/report/tables"
    bad = distance_check(ends, dist, contrast) + paired_check(dist, paired)
    mine = {(r["ds"], r["metric"], r["tool"]): r for r in subs_rows if r["subset"] == "all"}
    for metric, fname, col, run in (("yield", "ltrkin_yield.tsv", "yield", "e2e"), ("recall", "recall.tsv", "recall", None)):
        rows = [r for r in read_table(f"{tab}/{fname}") if r["ds"] in dss and r["status"] == "ok"
                and (run is None or r["run"] == run)]
        theirs = {(r["ds"], r["tool"]) for r in rows}
        ours = {(d, t) for (d, m, t) in mine if m == metric}
        if theirs != ours:
            bad.append(f"{fname}: tools with status ok differ: table only {sorted(theirs - ours)}, "
                       f"lead_checks only {sorted(ours - theirs)}")
        for r in rows:
            m = mine.get((r["ds"], metric, r["tool"]))
            if m is None:
                continue
            for mk, tk in (("value", col), ("lo", f"{col}_lo"), ("hi", f"{col}_hi"), ("dlo", f"{col}_dlo"),
                           ("dhi", f"{col}_dhi")):
                if score.fmt(m[mk]) != r[tk]:
                    bad.append(f"{fname} {r['ds']} {r['tool']} {tk}: table {r[tk]}, lead_checks {score.fmt(m[mk])}")
    blen = {(r["ds"], r["tool"], r["len_bin"]): r for r in read_table(f"{tab}/by_length.tsv")
            if r["ds"] in dss and r["status"] == "ok"}
    for r in lengths:
        if r["subset"] != "all":
            continue
        b = blen.get((r["ds"], r["tool"], r["len_bin"]))
        if b is None or (b["n_loci"], b["found"]) != (str(r["n_loci"]), str(r["n_found"])):
            bad.append(f"by_length.tsv {r['ds']} {r['tool']} {r['len_bin']}: table {b and (b['n_loci'], b['found'])}, "
                       f"lead_checks {(r['n_loci'], r['n_found'])}")
    for (d, m, t), r in mine.items():
        if m != "yield" or t == "all":
            continue
        e = [x for x in ends if x["ds"] == d and x["tool"] == t]
        if (sum(x["n_found"] for x in e), sum(x["n_fe"] for x in e)) != (mine[(d, "recall", t)]["n"], r["n"]):
            bad.append(f"ends.tsv {d} {t}: found/F/E sums differ from recall/yield counts")
        if t != REF:
            g = [x for x in gaps if x["ds"] == d and x["subset"] == "all" and x["tool"] == t][0]
            parts = g["tool_missed"] + g["tool_partial"] + g["tool_found_other_genotype"]
            if g["ref_fe_tool_not"] - g["tool_fe_ref_not"] != mine[(d, "yield", REF)]["n"] - r["n"] \
                    or parts != g["ref_fe_tool_not"]:
                bad.append(f"gap.tsv {d} {t}: net count or parts disagree with the F/E counts")
    if bad:
        die("self-check failed: lead_checks no longer reproduces score.py's tables in " + tab + ":\n  "
            + "\n  ".join(bad[:20]))
    n_yield = sum(1 for k in mine if k[1] == "yield")
    return (f"self-check passed: {n_yield} yield and {n_yield} recall rows equal ltrkin_yield.tsv (e2e) and recall.tsv; "
            f"{sum(r['subset'] == 'all' for r in lengths)} length rows equal by_length.tsv; ends and gap sums agree; "
            f"end_distance bins of {sum(r['tool'] != 'all' for r in dist) // len(DIST_BINS)} tool runs sum to "
            f"ends.tsv, group and contrast rows to their tools; paired cells and far calls consistent")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ds", nargs="+", default=["camara_hap2", "camara_hap1"], help="datasets to check")
    ap.add_argument("--root", default=LG, help="ltrkin_geno root, as score.py --root (default: production LG)")
    ap.add_argument("--running", action="append", default=[], metavar="DS:TOOL",
                    help="a tool run still in progress, as score.py --running (repeatable)")
    ap.add_argument("--out-dir", default=f"{LG}/report/checks",
                    help="output directory for subsets.tsv, gap.tsv, ends.tsv, length.tsv, end_distance.tsv, "
                         "end_distance_contrast.tsv and end_distance_paired.tsv")
    ap.add_argument("--fig1-numbers", default=None, metavar="TSV",
                    help="Fig. 1 numbers table for the non-fatal comparison (default: <root>/figures/"
                         "fig1_ltrkin_yield/fig1_ltrkin_yield_numbers.tsv; skipped if absent)")
    ap.add_argument("-v", "--verbose", action="store_true", help="print subset sizes per dataset on stderr")
    a = ap.parse_args()
    running = set()
    for x in a.running:
        ds, sep, t = x.partition(":")
        if not sep or ds not in a.ds or t not in ALL_TOOLS:
            die(f"--running expects DS:TOOL with DS among --ds and a known tool, got {x!r}")
        running.add((ds, REF if t == "ltr_peel_hc" else t))
    gaps, subs_rows, ends, lengths, counts, tsd, loci = [], [], [], [], {}, {}, {}
    for ds in a.ds:
        if ds not in DATASETS:
            die(f"unknown dataset {ds!r}")
        fai = str(DATASETS[ds]["ref"]) + ".fai"
        pool, _, res, ok, n = score.score_dataset(a.root, ds, fai, running)
        if REF not in ok:
            die(f"{ds}: {REF} did not finish; nothing to compare against")
        outcomes = {t: [o[0] for o in res[t]["outcome"]["e2e"]] for t in ok}
        status = {t: [s[0] for s in res[t]["status_rows"]] for t in ok}
        subs = subset_index(pool)
        gaps += gap_rows(ds, subs, outcomes, ok)
        subs_rows += subset_rows(ds, pool, subs, status, outcomes, ok, score.read_blocks(fai))
        ends += ends_rows(ds, res, outcomes, ok)
        lengths += length_rows(ds, pool, subs, status, ok)
        counts[ds] = distance_counts(res, outcomes, ok)
        tsd[ds] = tsd_calls(a.root, ds, ok)
        loci[ds] = locus_calls(res, outcomes, ok)
        if a.verbose:
            print(f"{ds}: {n} pooled loci; subsets {({k: len(v) for k, v in subs.items()})}; tools {len(ok)}",
                  file=sys.stderr)
    group = tsd_groups(tsd)
    dist = [r for ds in a.ds for r in distance_rows(ds, counts[ds], group, tsd[ds])]
    contrast = contrast_rows(dist)
    paired = [r for ds in a.ds for r in paired_rows(ds, loci[ds], group)]
    if a.verbose:
        print(f"no_tsd tools: {', '.join(sorted(t for t, g in group.items() if g == 'no_tsd'))}", file=sys.stderr)
    print(f"lead_checks: {self_check(a.root, a.ds, subs_rows, gaps, ends, lengths, dist, contrast, paired)}",
          file=sys.stderr)
    fig = fig1_compare(a.fig1_numbers or f"{a.root}/figures/fig1_ltrkin_yield/fig1_ltrkin_yield_numbers.tsv", dist,
                       contrast)
    if fig is None:
        print("lead_checks: Fig. 1 comparison skipped (numbers table not found)", file=sys.stderr)
    elif fig[0] == 0:
        print("lead_checks: Fig. 1 comparison: no overlapping values found (item text changed?)", file=sys.stderr)
    else:
        print(f"lead_checks: Fig. 1 comparison (not fatal): {fig[0] - len(fig[1])} of {fig[0]} overlapping values "
              f"equal" + "".join(f"\n  differs: {x}" for x in fig[1]), file=sys.stderr)
    os.makedirs(a.out_dir, exist_ok=True)
    score.write_table(f"{a.out_dir}/subsets.tsv", ["ds", "subset", "metric", "tool", "n_loci", "n", "value", "lo",
                                                   "hi", "d", "dlo", "dhi", "share"], subs_rows)
    score.write_table(f"{a.out_dir}/gap.tsv", ["ds", "subset", "tool", "ref_fe_tool_not", "tool_missed",
                                               "tool_partial", "tool_found_other_genotype", "tool_fe_ref_not"], gaps)
    score.write_table(f"{a.out_dir}/ends.tsv", ["ds", "tool", "end_bin", "n_found", "n_fe", "frac_fe"], ends)
    score.write_table(f"{a.out_dir}/length.tsv", ["ds", "subset", "len_bin", "tool", "n_loci", "n_found"], lengths)
    score.write_table(f"{a.out_dir}/end_distance.tsv", DIST_COLS, dist)
    score.write_table(f"{a.out_dir}/end_distance_contrast.tsv", CONTRAST_COLS, contrast)
    score.write_table(f"{a.out_dir}/end_distance_paired.tsv", PAIRED_COLS, paired)
    print(f"lead_checks: wrote subsets.tsv, gap.tsv, ends.tsv, length.tsv, end_distance.tsv, end_distance_contrast.tsv "
          f"and end_distance_paired.tsv to {a.out_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
