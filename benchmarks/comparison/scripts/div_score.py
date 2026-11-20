#!/usr/bin/env python3
"""Approach 1b (post hoc): divergence series. Per tool and timepoint: recall, precision and F1 against the PrinTE
standard intact set; recall per element LTR-divergence bin; d50, the LTR K2P at which recall falls to half its burn-in
value.

Truth: sims/truth/<ds>.std.bed (all whole, un-nested elements in this design) with LTR lengths from <ds>.truth.tsv.
Recall (primary) counts only elements with an LTR pair (ltr_len > 0): ~4% of PrinTE's standard set come from library
entries with no detectable LTR (e.g. *_INT internal-only entries), which have no LTR divergence and no LTR structure to
find. recall_std (all standard elements, as in approach 1) is kept beside it; precision uses the whole standard set.
Element divergence: the two LTRs compared position by position in that timepoint's genome. PrinTE mutates by
substitution only, so the LTRs keep their length and register. Matching: score_elements.match (one-to-one, reciprocal
overlap >= 0.9, the primary rule of approach 1).
Empty-candidate exits: a run that exits non-zero because its own candidate set is empty (verified from the tool's
intermediate files, see empty_candidates) is scored as 0 calls (recall 0, precision NA), status no_candidates; the
same reading as approach 2's dante_ltr vertebrate runs ("0 elements under defaults"). Any other failure stays failed
(no data for that timepoint).

TSD artifact (found after the design froze): in the first aging step PrinTE's mutator gives TE bases a
random share of the step but all other bases, the TSDs included, the full step. An element's two TSD copies therefore
run ~0.02 substitutions per site ahead of its LTR pair from timepoint 1 on, which real elements (TSDs and LTRs the same
age) do not show. Secondary outputs: sim_divergence_tsd.tsv splits recall by whether the two TSD copies are still
identical; recall_tsd_adj (and d50_k2p_tsd_adj) reweights those two strata to the share of identical TSD pairs expected
if each element's TSD copies had diverged as much as its LTR copies (see tsd_adjusted). The primary metrics are unchanged.

Writes <out>/sim_divergence.tsv, sim_divergence_bins.tsv, sim_divergence_d50.tsv, sim_divergence_tsd.tsv.
"""
import argparse
import csv
import glob
import math
import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from score_elements import match, read_bed3, read_calls  # noqa: E402

TOOLS = ["ltr_peel", "ltrharvest", "ltr_finder", "ltr_retriever", "look4ltrs", "hirepeat", "hite", "ltrpred",
         "dante_ltr", "ltr_peel_hc"]
STEP = 1540000          # generations per snapshot (sims/div/params.tsv)
BINW = 0.02             # K2P bin width
RO = 0.9
TS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}


def snap(k):
    return "burnin" if k == 0 else f"gen{k * STEP}_final"


def fasta_records(path):
    """Stream (name, sequence) one record at a time."""
    name, chunks = None, []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip().upper())
    if name is not None:
        yield name, "".join(chunks)


def pdist(a, b):
    n = d = 0
    for x, y in zip(a, b):
        if x in "ACGT" and y in "ACGT":
            n += 1
            d += x != y
    return d / n if n else float("nan")


def k2p(a, b):
    n = ts = tv = 0
    for x, y in zip(a, b):
        if x not in "ACGT" or y not in "ACGT":
            continue
        n += 1
        if x != y:
            if (x, y) in TS:
                ts += 1
            else:
                tv += 1
    if n == 0:
        return float("nan")
    p, q = ts / n, tv / n
    if 1 - 2 * p - q <= 0 or 1 - 2 * q <= 0:
        return float("inf")
    return -0.5 * math.log(1 - 2 * p - q) - 0.25 * math.log(1 - 2 * q) + 0.0   # + 0.0: no -0.0 for identical LTRs


def element_stats(fasta, elems):
    """Per element: (K2P between its two LTRs, are its two TSD copies identical, LTR p-distance). elems: (chrom, start,
    end, ltr_len,
    tsd_len), 0-based half-open element span without the TSDs, which flank it: [start - tsd_len, start) and
    [end, end + tsd_len). TSD flag None when unknown (no TSD, or a copy runs off the sequence)."""
    by = {}
    for e in elems:
        by.setdefault(e[0], []).append(e)
    out = {}
    for chrom, seq in fasta_records(fasta):
        for c, s, e, L, t in by.get(chrom, []):
            if L > 0 and e - s >= 2 * L:
                d, p = k2p(seq[s:s + L], seq[e - L:e]), pdist(seq[s:s + L], seq[e - L:e])
            else:
                d = p = float("nan")
            same = seq[s - t:s] == seq[e:e + t] if t > 0 and s - t >= 0 and e + t <= len(seq) else None
            out[(c, s, e)] = (d, same, p)
    return out


def fmt(x):
    return f"{x:.4f}" if math.isfinite(x) else "NA"


def tsd_adjusted(strata, q):
    """Recall with its two TSD strata reweighted to q. strata: [n_identical, detected_identical, n_differ,
    detected_differ]; q: share of elements whose TSD copies would still be identical had they diverged only as much as
    the element's LTR copies (mean of (1 - p_LTR) ** tsd_len), as for real elements, whose TSDs and LTRs are the same
    age."""
    n_id, d_id, n_df, d_df = strata
    if n_id and n_df:
        return d_id / n_id * q + d_df / n_df * (1 - q)
    if n_id and q > 0.999:       # burn-in: every TSD pair identical, nothing to reweight
        return d_id / n_id
    return float("nan")


def wilson(k, n, z=1.96):
    if n == 0:
        return float("nan"), float("nan")
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - h), min(1.0, c + h)


def _text(path):
    try:
        with open(path, errors="replace") as f:
            return f.read()
    except OSError:
        return ""


def empty_candidates(ds, tool):
    """True if a failed run's own files show it stopped on an empty candidate set:
    ltr_retriever: <ds>.fa.rawLTR.scn (LTRharvest + LTR_FINDER_parallel output) has no candidate line and LTR_retriever
      exits with "No candidate is found";
    look4ltrs: out/db.fasta (its element database) is empty and the run aborts in "Training Identity"."""
    rd = f"runs/{ds}/{tool}"
    if tool == "ltr_retriever":
        scn = _text(f"{rd}/{ds}.fa.rawLTR.scn")
        return (os.path.exists(f"{rd}/{ds}.fa.rawLTR.scn")
                and not any(l.strip() and not l.startswith("#") for l in scn.splitlines())
                and "No candidate is found" in _text(f"{rd}/run.log"))
    if tool == "look4ltrs":
        db = f"{rd}/out/db.fasta"
        return os.path.exists(db) and os.path.getsize(db) == 0 and "Training Identity" in _text(f"{rd}/run.log")
    return False


def status(ds, tool):
    rd = f"runs/{ds}/{'ltr_peel' if tool == 'ltr_peel_hc' else tool}"
    if os.path.exists(f"calls/{ds}/{tool}.tsv"):
        return "ok"
    if os.path.exists(f"{rd}/.failed"):
        return "no_candidates" if empty_candidates(ds, tool) else "failed"
    return "pending"


def d50(points):
    """points: [(k, median_k2p, recall)] sorted by k, starting at k = 0. Linear interpolation between timepoints."""
    if not points or points[0][0] != 0 or not points[0][2] > 0:
        return "NA"
    half = points[0][2] / 2
    for (k0, x0, r0), (k1, x1, r1) in zip(points, points[1:]):
        if k1 != k0 + 1:
            return "incomplete"
        if r1 <= half:
            return f"{x0 + (x1 - x0) * (r0 - half) / (r0 - r1):.4f}" if r0 != r1 else f"{x1:.4f}"
    return f">{points[-1][1]:.4f}" if points[-1][0] == 10 else "incomplete"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--bench", default=os.environ.get("BENCH_ROOT"), required="BENCH_ROOT" not in os.environ,
                    help="benchmark root (default: $BENCH_ROOT)")
    ap.add_argument("--out", default="report/tables", help="output directory, relative to --bench (default: %(default)s)")
    ap.add_argument("-v", "--verbose", action="store_true", help="per-dataset progress")
    a = ap.parse_args(argv)
    os.chdir(a.bench)
    dss = sorted(os.path.basename(p)[:-len(".truth.tsv")] for p in glob.glob("sims/truth/SD*_*.truth.tsv"))
    if not dss:
        sys.exit("div_score: no sims/truth/SD*_*.truth.tsv found")
    rows, bins, series, series_adj, tsd = [], {}, {}, {}, {}
    nan3 = (float("nan"), None, float("nan"))
    for ds in dss:
        rep, k = int(ds[2]), int(ds.split("_")[1])
        with open(f"sims/truth/{ds}.truth.tsv") as f:
            tr = {(r["chrom"], int(r["start"]), int(r["end"])): r for r in csv.DictReader(f, delimiter="\t")
                  if r["std"] == "1"}
        std = read_bed3(f"sims/truth/{ds}.std.bed")
        key = [(t["chrom"], t["start"], t["end"]) for t in std]
        tl = [len(tr[x]["tsd"]) if x in tr and tr[x]["tsd"] not in ("", "NA") else 0 for x in key]
        es = element_stats(f"sims/div/D{rep}/{snap(k)}.fasta",
                           [x + (int(tr[x]["ltr_len"]) if x in tr else 0, t) for x, t in zip(key, tl)])
        tk = [es.get(x, nan3)[0] for x in key]
        ts = [es.get(x, nan3)[1] for x in key]
        known = [(same, (1 - es[x][2]) ** t) for x, d, same, t in zip(key, tk, ts, tl)
                 if math.isfinite(d) and same is not None]
        q = sum(w for _, w in known) / len(known) if known else float("nan")
        share = sum(1 for same, _ in known if same) / len(known) if known else float("nan")
        finite = [x for x in tk if math.isfinite(x)]
        med = statistics.median(finite) if finite else float("nan")
        n_ltr = len(finite)
        if a.verbose:
            print(f"div_score: {ds} truth {len(std)}, K2P median {med:.4f} ({len(std) - len(finite)} without K2P), "
                  f"TSD pairs identical {share:.3f} (expected at LTR divergence {q:.3f})", file=sys.stderr)
        for tool in TOOLS:
            st = status(ds, tool)
            base = dict(ds=ds, seed=rep, timepoint=k, generations=k * STEP, median_k2p=f"{med:.4f}", tool=tool,
                        status=st, n_truth=n_ltr, n_truth_std=len(std), tsd_identical_share=fmt(share),
                        tsd_identical_expected=fmt(q))
            if st not in ("ok", "no_candidates"):
                rows.append(base)
                continue
            calls, c2t, t2c = [], {}, {}
            if st == "ok":
                calls = read_calls(f"calls/{ds}/{tool}.tsv")
                c2t, t2c = match(calls, std, RO)
            tp, nc = len(t2c), len(calls)
            tp_ltr = sum(1 for ti in t2c if math.isfinite(tk[ti]))
            p = tp / nc if nc else float("nan")
            r = tp_ltr / n_ltr if n_ltr else float("nan")
            f1 = 2 * p * r / (p + r) if nc and (p + r) > 0 else float("nan")
            exact = sum(calls[ci]["start"] == std[ti]["start"] and calls[ci]["end"] == std[ti]["end"]
                        for ci, ti in c2t.items())
            strata = [0, 0, 0, 0]
            for ti, x in enumerate(tk):
                if not math.isfinite(x):
                    continue
                hit = ti in t2c
                n_d = bins.setdefault((tool, int(x / BINW)), [0, 0])
                n_d[0] += 1
                n_d[1] += hit
                if ts[ti] is not None:
                    n_d = tsd.setdefault((tool, rep, k, med, ts[ti]), [0, 0])
                    n_d[0] += 1
                    n_d[1] += hit
                    strata[0 if ts[ti] else 2] += 1
                    strata[1 if ts[ti] else 3] += hit
            adj = tsd_adjusted(strata, q)
            rows.append(dict(base, n_calls=nc, tp=tp_ltr, precision=fmt(p), recall=fmt(r), f1=fmt(f1),
                             recall_std=fmt(tp / len(std)), exact_both=fmt(exact / tp) if tp else "NA",
                             recall_tsd_adj=fmt(adj)))
            series.setdefault((tool, rep), []).append((k, med, r))
            if math.isfinite(adj):
                series_adj.setdefault((tool, rep), []).append((k, med, adj))
    os.makedirs(a.out, exist_ok=True)
    cols = ["ds", "seed", "timepoint", "generations", "median_k2p", "tool", "status", "n_truth", "n_calls", "tp",
            "precision", "recall", "f1", "n_truth_std", "recall_std", "exact_both", "recall_tsd_adj",
            "tsd_identical_share", "tsd_identical_expected"]
    with open(f"{a.out}/sim_divergence.tsv", "w") as f:
        f.write("\t".join(cols) + "\tlabel\n")
        for r in rows:
            f.write("\t".join(str(r.get(c, "NA")) for c in cols) + "\tpost_hoc\n")
    with open(f"{a.out}/sim_divergence_bins.tsv", "w") as f:
        f.write("tool\tk2p_lo\tk2p_hi\tn_element_timepoints\tdetected\trecall\tci_lo\tci_hi\tlabel\n")
        for (tool, b), (n, d) in sorted(bins.items(), key=lambda x: (TOOLS.index(x[0][0]), x[0][1])):
            lo, hi = wilson(d, n)
            f.write(f"{tool}\t{b * BINW:.2f}\t{(b + 1) * BINW:.2f}\t{n}\t{d}\t{d / n:.4f}\t{lo:.4f}\t{hi:.4f}\t"
                    f"post_hoc\n")
    with open(f"{a.out}/sim_divergence_d50.tsv", "w") as f:
        f.write("tool\tseed\trecall_burnin\td50_k2p\td50_k2p_tsd_adj\ttimepoints_scored\tlabel\n")
        for (tool, rep), pts in sorted(series.items(), key=lambda x: (TOOLS.index(x[0][0]), x[0][1])):
            pts.sort()
            adj = sorted(series_adj.get((tool, rep), []))
            rb = f"{pts[0][2]:.4f}" if pts[0][0] == 0 else "NA"
            f.write(f"{tool}\t{rep}\t{rb}\t{d50(pts)}\t{d50(adj)}\t{len(pts)}\tpost_hoc\n")
    with open(f"{a.out}/sim_divergence_tsd.tsv", "w") as f:
        f.write("tool\tseed\ttimepoint\tmedian_k2p\ttsd_copies\tn_elements\tdetected\trecall\tci_lo\tci_hi\tlabel\n")
        for (tool, rep, k, med, same), (n, d) in sorted(tsd.items(), key=lambda x: (TOOLS.index(x[0][0]),) + x[0][1:]):
            lo, hi = wilson(d, n)
            f.write(f"{tool}\t{rep}\t{k}\t{med:.4f}\t{'identical' if same else 'differ'}\t{n}\t{d}\t{d / n:.4f}\t"
                    f"{lo:.4f}\t{hi:.4f}\tpost_hoc\n")
    n_nc = sum(r["status"] == "no_candidates" for r in rows)
    print(f"div_score: {len(dss)} datasets, {sum(r['status'] == 'ok' for r in rows)} tool runs scored, "
          f"{n_nc} empty-candidate exits scored as 0 calls -> {a.out}/", file=sys.stderr)


if __name__ == "__main__":
    main()
