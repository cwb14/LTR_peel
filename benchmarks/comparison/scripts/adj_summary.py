"""Approach 4 (post hoc) summaries from eval/adj/<ds>/<tool>.calls.tsv (rules frozen before any per-tool number).

element_library   precision_lib, pooled recall (RO >= 0.9; 0.5 secondary), F1; ext pool adds coding-unit loci;
                  95% CI by 1,000 bootstrap resamples of 5 Mb blocks
unsupported       per tool: evidence rates (supported / unsupported / null), f_real with bootstrap CI, bracket
age_bias          coding-unit loci pooled over tools: supported share per age bin (Wilson CI) + logistic slope
"""
import bisect
import csv
import math
import os
import statistics
from collections import defaultdict

import numpy as np

AGES = ["<0.01", "0.01-0.03", "0.03-0.1", ">=0.1", "NA"]
BLOCK = 5_000_000


def load(path):
    rows = []
    for r in csv.DictReader(open(path), delimiter="\t"):
        r["start"], r["end"] = int(r["start"]), int(r["end"])
        rows.append(r)
    return rows


def ro_ok(s1, e1, s2, e2, t):
    """reciprocal overlap >= t: overlap >= t x the longer span"""
    ov = min(e1, e2) - max(s1, s2)
    return ov > 0 and ov >= t * max(e1 - s1, e2 - s2)


class Spans:
    """per-chrom sorted spans for overlap queries"""

    def __init__(self, spans):                    # spans: [(chrom, s, e, payload)]
        self.by = defaultdict(list)
        for c, s, e, p in spans:
            self.by[c].append((s, e, p))
        self.starts, self.maxlen = {}, {}
        for c, v in self.by.items():
            v.sort()
            self.starts[c] = [x[0] for x in v]
            self.maxlen[c] = max(e - s for s, e, _ in v)

    def hits(self, c, s, e, t):
        v = self.by.get(c)
        if not v:
            return []
        i = bisect.bisect_left(self.starts[c], s - self.maxlen[c])
        out = []
        while i < len(v) and v[i][0] < e:
            if ro_ok(s, e, v[i][0], v[i][1], t):
                out.append(v[i][2])
            i += 1
        return out


def loci(spans, t=0.9):
    """union-find of spans [(chrom, s, e, key)] at RO >= t -> {key: locus id}, n loci"""
    parent = list(range(len(spans)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    idx = Spans([(c, s, e, i) for i, (c, s, e, _) in enumerate(spans)])
    for i, (c, s, e, _) in enumerate(spans):
        for j in idx.hits(c, s, e, t):
            if j != i:
                parent[find(i)] = find(j)
    roots = {}
    lab = {}
    for i, (_, _, _, key) in enumerate(spans):
        lab[key] = roots.setdefault(find(i), len(roots))
    return lab, len(roots)


def found_loci(member_spans, tool_rows, t):
    """set of locus ids that a tool's calls match (RO >= t to any member)"""
    idx = Spans(member_spans)                     # payload = locus id
    hit = set()
    for r in tool_rows:
        hit.update(idx.hits(r["chrom"], r["start"], r["end"], t))
    return hit


def block(c, s):
    return (c, s // BLOCK)


def element_library(rows_by_tool, pool_tools, n_boot=1000, seed=1):
    """rows_by_tool: {tool: rows}. Returns {tool: dict of metrics}"""
    res = {}
    for ext in (False, True):
        mem = []                                  # (chrom, s, e, key)
        for t in pool_tools:
            for k, r in enumerate(rows_by_tool.get(t, [])):
                if r["supported"] == "1" or (ext and r["coding"] == "1"):
                    mem.append((r["chrom"], r["start"], r["end"], (t, k)))
        lab, n_loci = loci(mem)
        loc_block = {}
        member_spans = []
        for c, s, e, key in mem:
            loc_block.setdefault(lab[key], block(c, s))
            member_spans.append((c, s, e, lab[key]))
        blocks = sorted(set(loc_block.values()) | {block(r["chrom"], r["start"])
                                                   for rows in rows_by_tool.values() for r in rows})
        bpos = {b: i for i, b in enumerate(blocks)}
        loc_in_block = np.bincount([bpos[b] for b in loc_block.values()], minlength=len(blocks))
        rng = np.random.default_rng(seed)
        draws = rng.integers(0, len(blocks), size=(n_boot, len(blocks)))
        w = np.stack([np.bincount(d, minlength=len(blocks)) for d in draws])      # block multiplicities
        for t, rows in rows_by_tool.items():
            d = res.setdefault(t, {"n_calls": len(rows)})
            sup = np.array([r["supported"] == "1" for r in rows])
            cb = np.array([bpos[block(r["chrom"], r["start"])] for r in rows], dtype=int)
            d["precision_lib"] = float(sup.mean()) if len(rows) else math.nan
            calls_b = np.bincount(cb, minlength=len(blocks))
            sup_b = np.bincount(cb, weights=sup, minlength=len(blocks))
            for thr, tag in ((0.9, ""), (0.5, "_ro50")):
                hit = found_loci(member_spans, rows, thr)
                fb = np.bincount([bpos[loc_block[x]] for x in hit], minlength=len(blocks))
                name = ("recall_ext" if ext else "recall_lib") + tag
                d[name] = len(hit) / n_loci if n_loci else math.nan
                if tag:
                    continue
                P = (w @ sup_b) / np.maximum(w @ calls_b, 1)
                R = (w @ fb) / np.maximum(w @ loc_in_block, 1)
                F = np.where(P + R > 0, 2 * P * R / np.maximum(P + R, 1e-12), 0)
                f1 = "F1_ext" if ext else "F1_lib"
                p, r_ = d["precision_lib"], d[name]
                d[f1] = 2 * p * r_ / (p + r_) if p + r_ > 0 else 0.0
                d[name + "_ci"] = tuple(np.percentile(R, [2.5, 97.5]))
                d[f1 + "_ci"] = tuple(np.percentile(F, [2.5, 97.5]))
                if not ext:
                    d["precision_lib_ci"] = tuple(np.percentile(P, [2.5, 97.5]))
            d["n_loci_ext" if ext else "n_loci_lib"] = n_loci
    return res


def rate(rows, key):
    v = [r[key] == "1" for r in rows if r[key] in ("0", "1")]
    return (sum(v) / len(v), len(v)) if v else (math.nan, 0)


def f_real(rows, n_boot=1000, seed=1):
    """age-matched (r_unsup - r_null) / (r_sup - r_null) for 'coding'; weighted by unsupported counts"""
    groups = {}
    for a in AGES:
        s = np.array([r["coding"] == "1" for r in rows if r["age"] == a and r["supported"] == "1"], float)
        u = np.array([r["coding"] == "1" for r in rows if r["age"] == a and r["supported"] == "0"], float)
        n = np.array([r["null_coding"] == "1" for r in rows
                      if r["age"] == a and r["supported"] == "0" and r["null_coding"] in ("0", "1")], float)
        if len(s) >= 10 and len(u) >= 10 and len(n):
            groups[a] = (s, u, n)

    def est(gs):
        num = den = 0.0
        for s, u, n in gs:
            d = s.mean() - n.mean()
            if d <= 0:
                continue
            num += len(u) * (u.mean() - n.mean()) / d
            den += len(u)
        return num / den if den else math.nan

    point = est(groups.values())
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(n_boot):
        gs = [(rng.choice(s, len(s)), rng.choice(u, len(u)), rng.choice(n, len(n))) for s, u, n in groups.values()]
        boots.append(est(gs))
    boots = [b for b in boots if not math.isnan(b)]
    ci = tuple(np.percentile(boots, [2.5, 97.5])) if boots else (math.nan, math.nan)
    return point, ci, sorted(groups)


def unsupported(rows, tool, n_boot=1000):
    sup = [r for r in rows if r["supported"] == "1"]
    uns = [r for r in rows if r["supported"] == "0"]
    nul = [{k[5:]: v for k, v in r.items() if k.startswith("null_")} for r in uns if r["null_shift"] != "NA"]
    d = {"n_calls": len(rows), "n_supported": len(sup), "n_unsupported": len(uns), "n_null": len(nul)}
    for grp, rs in (("sup", sup), ("unsup", uns), ("null", nul)):
        for k in ("coding", "junk", "gene", "nonltr", "sat", "trf"):
            d[f"{grp}_{k}"] = rate(rs, k)[0]
    for grp, rs in (("sup", sup), ("unsup", uns)):
        for k in ("tgca", "tsd"):
            d[f"{grp}_{k}"] = rate(rs, k)[0]
        ks = [float(r["k2p"]) for r in rs if r["k2p"] != "NA"]
        d[f"{grp}_k2p_median"] = statistics.median(ks) if ks else math.nan
    cats = defaultdict(int)
    for r in uns:
        cats[("coding" if r["coding"] == "1" else "") + ("+junk" if r["junk"] == "1" else "")] += 1
    d["unsup_coding_only"] = cats["coding"]
    d["unsup_junk_only"] = cats["+junk"]
    d["unsup_coding_and_junk"] = cats["coding+junk"]
    d["unsup_neither"] = cats[""]
    fr, ci, bins = f_real(rows, n_boot)
    d["f_real"], d["f_real_ci"], d["f_real_bins"] = fr, ci, ",".join(bins)
    gap = d["sup_coding"] - d["null_coding"]
    flags = []
    if tool == "dante_ltr":
        flags.append("circular")
    if not (gap >= 0.1):
        flags.append("uninformative")
    d["flag"] = ";".join(flags) or "ok"
    p = len(sup) / len(rows) if rows else math.nan
    d["precision_lib"] = p
    d["precision_upper"] = p + min(max(fr, 0.0), 1.0) * (1 - p) if not math.isnan(fr) else math.nan
    return d


def wilson(k, n, z=1.96):
    if n == 0:
        return (math.nan, math.nan)
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return (c - h, c + h)


def logistic(x, y, iters=50):
    """slope and its Wald 95% CI of y ~ 1 + x (Newton-Raphson)"""
    X = np.column_stack([np.ones(len(x)), x])
    b = np.zeros(2)
    for _ in range(iters):
        p = 1 / (1 + np.exp(-X @ b))
        W = p * (1 - p)
        H = X.T @ (X * W[:, None])
        step = np.linalg.solve(H + 1e-9 * np.eye(2), X.T @ (y - p))
        b += step
        if np.abs(step).max() < 1e-8:
            break
    p = 1 / (1 + np.exp(-X @ b))
    cov = np.linalg.inv(X.T @ (X * (p * (1 - p))[:, None]))
    se = math.sqrt(cov[1, 1])
    return b[1], (b[1] - 1.96 * se, b[1] + 1.96 * se)


def age_bias(rows_by_tool, pool_tools):
    mem = []
    for t in pool_tools:
        for k, r in enumerate(rows_by_tool.get(t, [])):
            if r["coding"] == "1":
                mem.append((r["chrom"], r["start"], r["end"], (t, k)))
    lab, n = loci(mem)
    by = defaultdict(list)
    for c, s, e, (t, k) in mem:
        by[lab[(t, k)]].append(rows_by_tool[t][k])
    loc = []
    for rs in by.values():
        sup = sum(r["supported"] == "1" for r in rs) >= len(rs) / 2
        ks = [float(r["k2p"]) for r in rs if r["k2p"] != "NA"]
        loc.append((sup, statistics.median(ks) if ks else None))
    out = []
    for name, lo, hi in (("<0.01", 0, 0.01), ("0.01-0.03", 0.01, 0.03), ("0.03-0.1", 0.03, 0.1),
                         (">=0.1", 0.1, math.inf), ("NA", None, None)):
        sel = [s for s, k in loc if (k is None if lo is None else (k is not None and lo <= k < hi))]
        kk = sum(sel)
        out.append((name, len(sel), kk / len(sel) if sel else math.nan, *wilson(kk, len(sel))))
    xy = [(math.log10(k + 0.001), float(s)) for s, k in loc if k is not None]
    slope, ci = logistic(np.array([x for x, _ in xy]), np.array([y for _, y in xy])) if len(xy) >= 20 else \
        (math.nan, (math.nan, math.nan))
    return out, n, slope, ci


def load_ds(adj_dir, tools):
    out = {}
    for t in tools:
        p = os.path.join(adj_dir, f"{t}.calls.tsv")
        if os.path.exists(p):
            out[t] = load(p)
    return out
