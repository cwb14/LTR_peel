#!/usr/bin/env python3
"""Score LTR-RT element calls against PrinTE truth (printe_truth.py outputs).

Matching: one-to-one, greedy on min(reciprocal overlap) among call-truth pairs with both overlaps >= --ro.
Reports (TSV, one row per metric scope):
  std       PrinTE standard set (OUT.std.bed: intact + tenants; PrinTE counts host pieces as fragments)
  ext       extended set (OUT.truth.tsv: intact, tenant, host, host_tenant; hosts rebuilt from their pieces)
  ext:<class> / ext:depth>=2 / ext:host_tenant_in_ltr   recall only, per stratum
Precision is over all calls; recall over the scope's truth. Also: boundary errors of matched pairs (left and right
genome ends; exact / <=5 / <=20 bp), false-positive anatomy (what each unmatched call hits most), and, if calls carry a
parent column, nesting-link precision/recall against OUT.pairs.tsv (LTR-RT tenants only).

Calls TSV (calls_to_bed.py): header with at least chrom, start (0-based), end; optional id, parent (id of the call's
host call, '.' if none).
"""
import argparse
import bisect
import csv
import sys
from collections import Counter, defaultdict


def read_calls(path):
    out = []
    with open(path) as f:
        for i, r in enumerate(csv.DictReader(f, delimiter="\t")):
            out.append(dict(chrom=r["chrom"], start=int(r["start"]), end=int(r["end"]),
                            id=r.get("id") or f"call{i}", parent=r.get("parent", ".") or "."))
    return out


def read_truth(path):
    with open(path) as f:
        return [dict(r, start=int(r["start"]), end=int(r["end"]), depth=int(r["depth"]),
                     tenant_in_ltr=int(r["tenant_in_ltr"]), std=int(r["std"])) for r in csv.DictReader(f, delimiter="\t")]


def read_bed3(path):
    out = []
    with open(path) as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) >= 3:
                out.append(dict(chrom=c[0], start=int(c[1]), end=int(c[2]), cat=c[3] if len(c) > 3 else ""))
    return out


class Index:
    def __init__(self, items):
        self.by = defaultdict(list)
        for k, x in enumerate(items):
            self.by[x["chrom"]].append((x["start"], x["end"], k))
        self.maxlen = {}
        for c, v in self.by.items():
            v.sort()
            self.maxlen[c] = max(e - s for s, e, _ in v)
        self.starts = {c: [s for s, _, _ in v] for c, v in self.by.items()}

    def overlaps(self, chrom, s, e):
        if chrom not in self.by:
            return []
        v, st = self.by[chrom], self.starts[chrom]
        i = bisect.bisect_left(st, s - self.maxlen[chrom])
        out = []
        while i < len(v) and v[i][0] < e:
            a, b, k = v[i]
            if b > s:
                out.append(k)
            i += 1
        return out


def match(calls, truth, ro):
    idx = Index(truth)
    cand = []
    for ci, c in enumerate(calls):
        lc = c["end"] - c["start"]
        for ti in idx.overlaps(c["chrom"], c["start"], c["end"]):
            t = truth[ti]
            ov = min(c["end"], t["end"]) - max(c["start"], t["start"])
            r1, r2 = ov / lc, ov / (t["end"] - t["start"])
            if r1 >= ro and r2 >= ro:
                cand.append((min(r1, r2), ci, ti))
    cand.sort(reverse=True)
    c2t, t2c = {}, {}
    for _, ci, ti in cand:
        if ci not in c2t and ti not in t2c:
            c2t[ci], t2c[ti] = ti, ci
    return c2t, t2c


def prf(tp_c, n_c, tp_t, n_t):
    p = tp_c / n_c if n_c else float("nan")
    r = tp_t / n_t if n_t else float("nan")
    f = 2 * p * r / (p + r) if n_c and n_t and (p + r) > 0 else float("nan")
    return p, r, f


def anatomy(calls, unmatched, truth, feats, ro):
    """Category of what each unmatched call overlaps most (bp)."""
    tidx, fidx = Index(truth), Index(feats)
    cnt = Counter()
    for ci in unmatched:
        c = calls[ci]
        best, cat = 0, "no_feature"
        for ti in tidx.overlaps(c["chrom"], c["start"], c["end"]):
            t = truth[ti]
            ov = min(c["end"], t["end"]) - max(c["start"], t["start"])
            if ov > best:
                best, cat = ov, f"partial_{t['class']}"
        if best < 0.5 * (c["end"] - c["start"]):           # mostly outside intact LTR-RTs
            bf, fcat = 0, None
            for fi in fidx.overlaps(c["chrom"], c["start"], c["end"]):
                f = feats[fi]
                if f["cat"] in ("LTR_intact", "LTR_piece"):
                    continue
                ov = min(c["end"], f["end"]) - max(c["start"], f["start"])
                if ov > bf:
                    bf, fcat = ov, f["cat"]
            if fcat and bf > best:
                cat = fcat
        cnt[cat] += 1
    return cnt


def main(argv=None):
    ap = argparse.ArgumentParser(description="Score LTR-RT calls against PrinTE truth (one-to-one matching).")
    ap.add_argument("--calls", required=True, help="calls TSV (chrom, start 0-based, end[, id, parent])")
    ap.add_argument("--truth", required=True, help="printe_truth.py prefix (reads .truth.tsv .std.bed .features.bed "
                                                     ".pairs.tsv)")
    ap.add_argument("--ro", type=float, nargs="+", default=[0.9, 0.8, 0.95], help="reciprocal-overlap thresholds; "
                    "the first is primary (default 0.9 0.8 0.95)")
    ap.add_argument("--label", default="", help="label for the 'label' column (e.g. tool/dataset)")
    ap.add_argument("--out", required=True, help="output prefix (.metrics.tsv, .boundaries.tsv, .anatomy.tsv)")
    ap.add_argument("-v", "--verbose", action="store_true", help="print the primary metrics")
    a = ap.parse_args(argv)
    calls = read_calls(a.calls)
    truth = read_truth(a.truth + ".truth.tsv")
    std = read_bed3(a.truth + ".std.bed")
    feats = read_bed3(a.truth + ".features.bed")
    pairs = [r for r in csv.DictReader(open(a.truth + ".pairs.tsv"), delimiter="\t") if r["tenant_class"] == "LTR"]
    rows = []
    for ro in a.ro:
        c2t_s, t2c_s = match(calls, std, ro)
        p, r, f = prf(len(c2t_s), len(calls), len(t2c_s), len(std))
        rows.append(("std", f"{ro:g}", len(calls), len(std), len(c2t_s), p, r, f))
        c2t, t2c = match(calls, truth, ro)
        p, r, f = prf(len(c2t), len(calls), len(t2c), len(truth))
        rows.append(("ext", f"{ro:g}", len(calls), len(truth), len(c2t), p, r, f))
        strata = defaultdict(list)
        for ti, t in enumerate(truth):
            strata[f"ext:{t['class']}"].append(ti)
            if t["depth"] >= 2:
                strata["ext:depth>=2"].append(ti)
            if t["tenant_in_ltr"]:
                strata["ext:host_tenant_in_ltr"].append(ti)
            if t["class"] in ("host", "host_tenant") and not t["tenant_in_ltr"]:
                strata["ext:host_tenant_internal"].append(ti)
        for k in sorted(strata):
            n_t = len(strata[k])
            tp = sum(1 for ti in strata[k] if ti in t2c)
            rows.append((k, f"{ro:g}", len(calls), n_t, tp, float("nan"), tp / n_t if n_t else float("nan"), float("nan")))
        if ro == a.ro[0]:
            prim = (c2t, t2c)
    c2t, t2c = prim
    with open(a.out + ".metrics.tsv", "w") as f:
        f.write("label\tscope\tro\tn_calls\tn_truth\ttp\tprecision\trecall\tf1\n")
        for row in rows:
            f.write(a.label + "\t" + "\t".join(f"{x:.4f}" if isinstance(x, float) else str(x) for x in row) + "\n")
    # boundaries (primary threshold, extended truth)
    dl = [calls[ci]["start"] - truth[ti]["start"] for ci, ti in c2t.items()]
    dr = [calls[ci]["end"] - truth[ti]["end"] for ci, ti in c2t.items()]
    with open(a.out + ".boundaries.tsv", "w") as f:
        f.write("label\tend\tn\texact\tle5\tle20\tmedian_abs\n")
        for name, d in (("left", dl), ("right", dr)):
            n = len(d)
            ab = sorted(abs(x) for x in d)
            fr = (lambda k: sum(1 for x in ab if x <= k) / n) if n else (lambda k: float("nan"))
            med = ab[n // 2] if n else float("nan")
            f.write(f"{a.label}\t{name}\t{n}\t{fr(0):.4f}\t{fr(5):.4f}\t{fr(20):.4f}\t{med}\n")
    # false-positive anatomy
    unmatched = [ci for ci in range(len(calls)) if ci not in c2t]
    cnt = anatomy(calls, unmatched, truth, feats, a.ro[0])
    with open(a.out + ".anatomy.tsv", "w") as f:
        f.write("label\tcategory\tn\tfrac_of_calls\n")
        for k, v in cnt.most_common():
            f.write(f"{a.label}\t{k}\t{v}\t{v / max(len(calls), 1):.4f}\n")
    # nesting links (only if the caller reports parents)
    has_links = any(c["parent"] != "." for c in calls)
    if has_links:
        id2ci = {c["id"]: ci for ci, c in enumerate(calls)}
        tid = {t["id"]: ti for ti, t in enumerate(truth)}
        truth_pairs = {(tid[p["host_id"]], tid[p["tenant_id"]]) for p in pairs
                       if p["host_id"] in tid and p["tenant_id"] in tid}
        called = [(id2ci.get(c["parent"]), ci) for ci, c in enumerate(calls) if c["parent"] != "."]
        called = [(h, t) for h, t in called if h is not None]
        ok = sum(1 for h, t in called if h in c2t and t in c2t and (c2t[h], c2t[t]) in truth_pairs)
        found = sum(1 for h, t in truth_pairs if h in t2c and t in t2c and
                    calls[t2c[t]]["parent"] == calls[t2c[h]]["id"])
        with open(a.out + ".metrics.tsv", "a") as f:
            p = ok / len(called) if called else float("nan")
            r = found / len(truth_pairs) if truth_pairs else float("nan")
            f.write(f"{a.label}\tnest_links\t{a.ro[0]:g}\t{len(called)}\t{len(truth_pairs)}\t{ok}\t{p:.4f}\t{r:.4f}\t"
                    f"{(2 * p * r / (p + r)) if p + r > 0 else float('nan'):.4f}\n")
    if a.verbose:
        for row in rows[:2]:
            print(a.label, row, file=sys.stderr)


if __name__ == "__main__":
    main()
