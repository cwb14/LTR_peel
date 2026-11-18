#!/usr/bin/env python3
"""Precision, recall and F1 by nesting class and by protein class on one basis in the nested simulations (S1, S2).

Basis ext (default), shared by both class systems: truth = every row of sims/truth/<S>.truth.tsv (intact, tenant, host,
host_tenant; hosts rebuilt from their pieces); one-to-one matching at reciprocal overlap >= 0.9 (score_elements.match),
run once per genome x tool inside sim_class_prf.class_counts and used by both systems.
Basis hasltr (--basis hasltr): the ext basis (same truth, matching and classes), then in both systems the copies of the
library entries without an ~LTRlen tag (has_ltr = 0) leave the truth and every call matched to one of them leaves the
calls (neither a true positive nor in the precision denominator); unmatched calls stay with their ext class.
system = nesting: sim_class_prf.py as is (host = host + host_tenant, tenant, unnested = intact; an unmatched call takes
  the class of its anatomy category). Rows identical to sim_class_prf.tsv (checked).
system = protein:
  true element: class of its PrinTE library entry, experiments/nonauto/lib/lib_class.tsv (classify_lib.py: coding = >= 2
    distinct REXdb Class_I|LTR domain names of one (strand, superfamily) inside the internal region; proteinless = no
    Class_I|LTR domain record in the entry; partial = the rest). The 68 entries without an ~LTRlen tag (has_ltr = 0)
    are classed by the same rule on the whole entry (classify_lib.py); their elements are kept and counted on stderr
    (--lib-check tests that class against the middle 80% and against LTRs found by self-alignment).
  call: nonauto_sim.call_class (coding unit in the call's middle 80%; proteinless = no Class_I|LTR domain fully inside
    the call; partial = the rest) on runs/<S>/dante_ltr/dante.gff3 with its IDs mapped back to the PrinTE names
    (sim_nonauto_prf.rename_doms); a call is true when matched to any ext element (as sim_nonauto_prf_ext.tsv).
Per genome and class: precision = tp_calls / n_calls, recall = tp_truth / n_truth, F1 = 2PR / (P + R); NA when undefined
(sim_class_prf.prf). Summary: mean and sample SD (n - 1) over the 3 genomes of a set, over finite values (n_genomes
beside), n_truth and n_calls summed (sim_class_prf_summary.tsv convention).
Coordinates: calls and truth TSV 0-based half-open; DANTE GFF3 1-based closed -> 0-based half-open in
evidence.dante_domains (start - 1, end); no conversion here.
Checks (fatal; both bases are computed in every run): ext: per genome x tool and per system, n_truth, tp_truth,
n_calls and tp_calls summed over classes equal the ext ro 0.9 row of eval/sim/<S>/<tool>.metrics.tsv; nesting rows =
sim_class_prf.tsv and nesting summary = sim_class_prf_summary.tsv (printed strings); protein n_calls, tp_calls,
precision per seed = sim_nonauto_prf_ext.tsv; old-basis protein recall recomputed here (std == 1, has_ltr, matched
against std == 1 rows) = its n_truth / tp_truth; cross-tab margins = class n_truth of both systems. hasltr: per genome
x tool the two systems' sums agree and equal the ext row minus the removed elements (n_truth) and removed calls (tp,
n_calls); removed calls = removed elements that were matched.
Outputs (--out): ext: sim_class_prf_joint.tsv, sim_class_prf_joint_summary.tsv, sim_class_joint_crosstab.tsv;
hasltr: sim_class_prf_joint_hasltr.tsv, sim_class_prf_joint_hasltr_summary.tsv (same schema).
"""
import argparse
import csv
import math
import os
import queue  # noqa: F401  (stdlib, cached before div_score, via nonauto_sim, puts scripts/ first: scripts/queue.py)
import multiprocessing.pool  # noqa: F401
import statistics
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from multiprocessing import Pool

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sim_class_prf as SC  # noqa: E402
import sim_nonauto_prf as SN  # noqa: E402
from evidence import coding_unit, dante_domains  # noqa: E402
from score_elements import match, read_bed3, read_calls, read_truth  # noqa: E402

N = SN.N                       # experiments/nonauto/nonauto_sim.py
B, NA = SC.B, SN.NA
TOOLS, DIAG, SETS, RO = SC.TOOLS, SC.DIAG, SC.SETS, SC.RO
SYSTEMS = {"nesting": SC.CLASSES, "protein": SN.CLS}
KEYS = ("n_truth", "tp_truth", "n_calls", "tp_calls")
METRICS = ("precision", "recall", "f1")
ORIGINS = ("old_set", "ext_only", "no_ltr_tag")   # -v only: std == 1 tagged / std == 0 tagged / entry without ~LTRlen
BASES = ("ext", "hasltr")
OUTS = {"ext": ("sim_class_prf_joint.tsv", "sim_class_prf_joint_summary.tsv"),
        "hasltr": ("sim_class_prf_joint_hasltr.tsv", "sim_class_prf_joint_hasltr_summary.tsv")}


def die(msg):
    sys.exit(f"error: {msg}")


def element_classes(truth, lib):
    """Protein class of every ext truth element = class of its library entry (lib_class.tsv); fails on an unknown family."""
    miss = sorted({t["family"] for t in truth if t["family"] not in lib})
    if miss:
        die(f"truth families missing from lib_class.tsv: {miss[:5]} ({len(miss)})")
    return [lib[t["family"]]["cls"] for t in truth]


def tally(truth, ecls, ncls, ccls, c2t, t2c, keep=None):
    """Per system and class: n_truth, tp_truth, n_calls, tp_calls. ecls: protein class per truth element; ncls / ccls:
    nesting / protein class per call. keep: truth mask (None = all); an element not kept leaves the truth and the call
    matched to it leaves the calls. Returns (counts, number of calls removed)."""
    k = {s: {c: dict.fromkeys(KEYS, 0) for c in cl} for s, cl in SYSTEMS.items()}
    for ti, t in enumerate(truth):
        if keep is None or keep[ti]:
            for s, c in (("nesting", SC.truth_class(t)), ("protein", ecls[ti])):
                k[s][c]["n_truth"] += 1
                k[s][c]["tp_truth"] += ti in t2c
    removed = 0
    for ci, (nc, pc) in enumerate(zip(ncls, ccls)):
        ti = c2t.get(ci)
        if ti is not None and keep is not None and not keep[ti]:
            removed += 1
            continue
        for s, c in (("nesting", nc), ("protein", pc)):
            k[s][c]["n_calls"] += 1
            k[s][c]["tp_calls"] += ti is not None
    return k, removed


def joint_counts(calls, truth, feats, ecls, doms, keep=None, ro=RO):
    """Both systems on one matching (that of sim_class_prf.class_counts). Nesting class of a call: that of its matched
    element, else of its anatomy category (sim_class_prf); protein class of a call: nonauto_sim.call_class on doms
    (evidence.dante_domains layout keyed by the calls' sequence names). Returns the counts on all elements, the counts
    with keep applied (= the first when keep is None), the number of calls removed, c2t and t2c."""
    nest, cats, c2t, t2c, _ = SC.class_counts(calls, truth, feats, ro)
    ncls = [SC.truth_class(truth[c2t[ci]]) if ci in c2t else SC.fp_class(cats[ci]) for ci in range(len(calls))]
    ccls = [N.call_class(doms, call) for call in calls]
    k, _ = tally(truth, ecls, ncls, ccls, c2t, t2c)
    if k["nesting"] != nest:
        raise RuntimeError("nesting tally differs from sim_class_prf.class_counts")
    kk, removed = (k, 0) if keep is None else tally(truth, ecls, ncls, ccls, c2t, t2c, keep)
    return k, kk, removed, c2t, t2c


def crosstab(truth, ecls, mask=None):
    """Counts of truth elements by (nesting class, protein class); mask selects elements (default all)."""
    return Counter((SC.truth_class(t), c) for i, (t, c) in enumerate(zip(truth, ecls)) if mask is None or mask[i])


def run_ds(args):
    """Pool task. die() raises SystemExit, which would kill a pool worker and leave pool.map waiting forever: hand it
    to the parent as an ordinary error instead."""
    try:
        return genome_counts(*args)
    except SystemExit as e:
        raise RuntimeError(f"{args[0]}: {e}") from None


def genome_counts(ds, lib):
    truth = read_truth(f"{B}/sims/truth/{ds}.truth.tsv")
    feats = read_bed3(f"{B}/sims/truth/{ds}.features.bed")
    ecls = element_classes(truth, lib)
    noltr = [lib[t["family"]]["has_ltr"] != "1" for t in truth]
    doms = SN.rename_doms(dante_domains(f"{B}/runs/{ds}/dante_ltr/dante.gff3"), f"{B}/data/{ds}/names.tsv")
    old_truth = N.truth(ds, lib)                     # old recall basis: std == 1 rows (has_ltr flag, library class)
    std_idx = [i for i, t in enumerate(truth) if t["std"] == 1]
    if len(std_idx) != len(old_truth) or any((truth[i]["chrom"], truth[i]["start"], truth[i]["end"]) !=
                                             (x["chrom"], x["start"], x["end"]) for i, x in zip(std_idx, old_truth)):
        die(f"{ds}: std == 1 rows of truth.tsv and nonauto_sim.truth differ")
    keep = [not z for z in noltr]                    # hasltr basis
    rows, bad, cmp, removed = {b: [] for b in BASES}, [], {}, {}
    for tool in TOOLS + DIAG:
        cp = f"{B}/calls/{ds}/{tool}.tsv"
        if not os.path.isfile(cp):
            die(f"calls not found: {cp}")
        calls = read_calls(cp)
        k, kh, n_rm_calls, c2t, t2c = joint_counts(calls, truth, feats, ecls, doms, keep)
        ext = SC.read_metrics(f"{B}/eval/sim/{ds}/{tool}.metrics.tsv")["ext"]
        exp = dict(n_truth=int(ext["n_truth"]), tp_truth=int(ext["tp"]), n_calls=int(ext["n_calls"]),
                   tp_calls=int(ext["tp"]))
        n_rm_truth = sum(noltr)
        if n_rm_calls != sum(1 for ti in t2c if noltr[ti]):
            bad.append(f"{ds} {tool}: removed calls {n_rm_calls} != matched removed elements")
        exp_h = dict(n_truth=exp["n_truth"] - n_rm_truth, tp_truth=exp["tp_truth"] - n_rm_calls,
                     n_calls=exp["n_calls"] - n_rm_calls, tp_calls=exp["tp_calls"] - n_rm_calls)
        removed[tool] = (n_rm_truth, n_rm_calls)
        for basis, kb, e in (("ext", k, exp), ("hasltr", kh, exp_h)):
            tot = {s: {x: sum(kk[c][x] for c in kk) for x in KEYS} for s, kk in kb.items()}
            if tot["nesting"] != tot["protein"]:
                bad.append(f"{ds} {tool} {basis}: nesting sums {tot['nesting']} != protein sums {tot['protein']}")
            for s in tot:
                if tot[s] != e:
                    bad.append(f"{ds} {tool} {basis} {s}: summed over classes {tot[s]} vs expected {e}")
        # old basis (sim_nonauto_prf_ext.tsv recall), the same elements under the ext matching, and the ext elements
        # by origin (old set, ext-only rows = rebuilt hosts etc. from tagged entries, entries without ~LTRlen), for -v
        _, t2c_old = match(calls, old_truth, RO)
        old = {c: [0, 0] for c in SN.CLS}
        same = {c: [0, 0] for c in SN.CLS}
        for j, x in enumerate(old_truth):
            if x["has_ltr"]:
                old[x["cls"]][0] += 1
                old[x["cls"]][1] += j in t2c_old
                same[x["cls"]][0] += 1
                same[x["cls"]][1] += std_idx[j] in t2c
        grp = {c: {g: [0, 0] for g in ORIGINS} for c in SN.CLS}
        for ti, t in enumerate(truth):
            g = "no_ltr_tag" if noltr[ti] else ("old_set" if t["std"] == 1 else "ext_only")
            grp[ecls[ti]][g][0] += 1
            grp[ecls[ti]][g][1] += ti in t2c
        cmp[tool] = (old, same, grp)
        for basis, kb in (("ext", k), ("hasltr", kh)):
            for system, cl in SYSTEMS.items():
                for c in cl:
                    p, r, f = SC.prf(kb[system][c])
                    rows[basis].append(dict(set=ds[:2], seed=ds[3], tool=tool, system=system, cls=c, **kb[system][c],
                                            precision=p, recall=r, f1=f, role="diagnostic" if tool in DIAG else "primary"))
    return ds, rows, bad, crosstab(truth, ecls), crosstab(truth, ecls, noltr), cmp, removed


def summarize(rows):
    g = defaultdict(list)
    for r in rows:
        g[(r["set"], r["tool"], r["system"], r["cls"])].append(r)
    out = []
    for (st, tool, system, c), rr in g.items():
        d = dict(set=st, tool=tool, system=system, cls=c)
        for x in METRICS:
            v = [r[x] for r in rr if math.isfinite(r[x])]
            d[f"{x}_mean"] = statistics.mean(v) if v else float("nan")
            d[f"{x}_sd"] = statistics.stdev(v) if len(v) > 1 else float("nan")
            d[f"{x}_n_genomes"] = len(v)
        d["n_truth"] = sum(r["n_truth"] for r in rr)
        d["n_calls"] = sum(r["n_calls"] for r in rr)
        d["role"] = rr[0]["role"]
        out.append(d)
    return out


def read_tsv(path):
    return list(csv.DictReader(open(path), delimiter="\t"))


def compare(mine, ref, key, cols, label):
    """Mismatches between my rows (formatted with sim_class_prf.fmt) and reference TSV rows on `cols`."""
    bad = []
    refd = {tuple(r[k] for k in key): r for r in ref}
    mined = {tuple(str(r[k]) for k in key): r for r in mine}
    if set(refd) != set(mined):
        bad.append(f"{label}: keys differ ({len(set(mined) - set(refd))} only here, {len(set(refd) - set(mined))} only "
                   f"in the reference)")
    for kk in sorted(set(refd) & set(mined)):
        for c in cols:
            if SC.fmt(mined[kk][c]) != refd[kk][c]:
                bad.append(f"{label} {kk} {c}: {SC.fmt(mined[kk][c])} vs {refd[kk][c]}")
    return bad


def terminal_repeat(hsps, L, tol):
    """Outer LTR length from self-alignment HSPs of one sequence of length L (1-based closed qs, qe, ss, se, plus
    strand, diagonal removed): the longest HSP whose arms start within tol bp of the 5' end and end within tol bp of the
    3' end without overlapping; length = the longer arm measured from its end. None if there is no such HSP."""
    term = [h for h in hsps if h[0] - 1 <= tol and L - h[3] <= tol and h[1] < h[2]]
    if not term:
        return None
    qs, qe, ss, se = max(term, key=lambda h: (h[1] - h[0], -h[0]))
    return max(qe, L - ss + 1)


def lib_class_of(doms, sid, s, e, n_dom):
    """classify_lib.py rule with the internal region [s, e)."""
    return "coding" if coding_unit(doms, sid, s, e) else ("proteinless" if n_dom == 0 else "partial")


def lib_check(blastn, tol, verbose):
    """Library entries without an ~LTRlen tag: class on the whole entry (lib_class.tsv) vs the middle 80% (call_class
    window) vs between the terminal direct repeats found by blastn self-alignment."""
    d = f"{NA}/lib"
    names = dict(l.rstrip("\n").split("\t") for l in open(f"{d}/names.tsv"))
    rows = {r["header"]: r for r in read_tsv(f"{d}/lib_class.tsv")}
    doms = dante_domains(f"{d}/dante.gff3")
    bad = []
    for sid, h in names.items():           # lib_class.tsv reproduces from lib/dante.gff3
        r = rows[h]
        L, ll = int(r["len"]), int(r["ltr_len"])
        s, e = (ll, L - ll) if ll else (0, L)
        if lib_class_of(doms, sid, s, e, int(r["n_dom"])) != r["cls"]:
            bad.append(h)
    if bad:
        die(f"lib_class.tsv does not reproduce for {len(bad)} entries: {bad[:3]}")
    sel = {sid: h for sid, h in names.items() if rows[h]["has_ltr"] != "1"}
    with tempfile.TemporaryDirectory() as td:
        fa = f"{td}/noltr.fa"
        with open(fa, "w") as f:
            keep = False
            for line in open(f"{d}/ltr_entries.fa"):
                if line.startswith(">"):
                    keep = line[1:].split()[0] in sel
                if keep:
                    f.write(line)
        out = subprocess.run([blastn, "-task", "blastn", "-dust", "no", "-strand", "plus", "-evalue", "1e-10",
                              "-query", fa, "-subject", fa, "-outfmt", "6 qseqid sseqid qstart qend sstart send"],
                             check=True, capture_output=True, text=True).stdout
    hsps = defaultdict(list)
    for line in out.splitlines():
        q, s, qs, qe, ss, se = line.split("\t")
        qs, qe, ss, se = int(qs), int(qe), int(ss), int(se)
        if q == s and ss > qs:             # off-diagonal self-hits, 5' arm as query
            hsps[q].append((qs, qe, ss, se))
    differ, n_rep, n_mid = [], 0, 0
    for sid, h in sel.items():
        r = rows[h]
        L, nd = int(r["len"]), int(r["n_dom"])
        m = L // 10
        mid = lib_class_of(doms, sid, m, L - m, nd)
        n_mid += mid == r["cls"]
        ltr = terminal_repeat(hsps[sid], L, tol)
        rep = None if ltr is None else lib_class_of(doms, sid, ltr, L - ltr, nd)
        n_rep += ltr is not None
        if mid != r["cls"] or (rep is not None and rep != r["cls"]):
            differ.append((r["family"], L, ltr, r["cls"], mid, rep))
        if verbose:
            print(f"{r['family']}\tlen {L}\tLTR {ltr}\twhole {r['cls']}\tmiddle80 {mid}\tbetween_repeats {rep}",
                  file=sys.stderr)
    fams = {x[0] for x in differ}
    n_el = Counter()
    for ds in SETS:
        n_el[ds[:2]] += sum(1 for r in read_tsv(f"{B}/sims/truth/{ds}.truth.tsv") if r["family"] in fams)
    print(f"lib_class.tsv reproduces from lib/dante.gff3 for all {len(names)} entries; {len(sel)} without an ~LTRlen tag:"
          f" middle-80% class = whole-entry class for {n_mid}; terminal direct repeat (blastn self-alignment, ends <= "
          f"{tol} bp) in {n_rep}, class between them differs for {sum(1 for x in differ if x[5] not in (None, x[3]))}",
          file=sys.stderr)
    for x in differ:
        print("  differs: %s len %d LTR %s whole %s middle80 %s between_repeats %s" % x, file=sys.stderr)
    print(f"  ext elements of the differing entries: S1 {n_el['S1']}, S2 {n_el['S2']} (3 genomes each)", file=sys.stderr)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", default=HERE, help="output dir (default: script dir)")
    ap.add_argument("--basis", choices=BASES, default="ext", help="ext (default): every truth.tsv element; hasltr: "
                    "without the copies of library entries lacking ~LTRlen and the calls matched to them "
                    "(writes sim_class_prf_joint_hasltr*.tsv); both bases are computed and checked in every run")
    ap.add_argument("-t", "--threads", type=int, default=4, help="worker processes, one genome each (default 4)")
    ap.add_argument("--lib-check", action="store_true", help="only test the class of library entries without an "
                    "~LTRlen tag (whole entry vs middle 80%% vs between self-aligned terminal repeats) and exit")
    ap.add_argument("--blastn", default=os.path.join(os.path.dirname(sys.executable), "blastn"),
                    help="blastn for --lib-check (default: next to this python)")
    ap.add_argument("--tol", type=int, default=50, help="--lib-check: max bp between a repeat arm and its end (50)")
    ap.add_argument("-v", "--verbose", action="store_true", help="per-genome counts and the old vs new protein recall "
                    "table on stderr")
    a = ap.parse_args(argv)
    if a.lib_check:
        lib_check(a.blastn, a.tol, a.verbose)
        return
    lib = N.load_lib(NA)
    print(f"sim_class_prf_joint: basis {a.basis}, {len(SETS)} genomes x {len(TOOLS + DIAG)} tools", file=sys.stderr)
    try:
        with Pool(min(a.threads, len(SETS))) as pool:
            res = pool.map(run_ds, [(ds, lib) for ds in SETS])
    except RuntimeError as e:
        die(str(e))
    bad = [b for r in res for b in r[2]]
    order = {t: i for i, t in enumerate(TOOLS + DIAG)}
    sysord = {s: i for i, s in enumerate(SYSTEMS)}
    rows_b = {b: sorted((x for r in res for x in r[1][b]), key=lambda r: (
        r["set"], r["seed"], order[r["tool"]], sysord[r["system"]], SYSTEMS[r["system"]].index(r["cls"]))) for b in BASES}
    summ_b = {b: sorted(summarize(rows_b[b]), key=lambda r: (
        r["set"], order[r["tool"]], sysord[r["system"]], SYSTEMS[r["system"]].index(r["cls"]))) for b in BASES}
    rows, summ = rows_b["ext"], summ_b["ext"]       # every reference check is on the ext basis
    nest = [r for r in rows if r["system"] == "nesting"]
    prot = [r for r in rows if r["system"] == "protein"]
    bad += compare(nest, read_tsv(f"{HERE}/sim_class_prf.tsv"), ("set", "seed", "tool", "cls"),
                   ("n_truth", "tp_truth", "recall", "n_calls", "tp_calls", "precision", "f1", "role"), "sim_class_prf.tsv")
    ref = read_tsv(f"{HERE}/sim_nonauto_prf_ext.tsv")
    ref_seed = [r for r in ref if r["seed"] != "pooled"]
    if any(r["precision_basis"] != "ext" for r in ref):
        die("sim_nonauto_prf_ext.tsv: precision_basis is not ext")
    bad += compare(prot, ref_seed, ("set", "seed", "tool", "cls"), ("n_calls", "tp_calls", "precision", "role"),
                   "sim_nonauto_prf_ext.tsv")
    cmp = {(ds, tool): v for ds, _, _, _, _, c, _ in res for tool, v in c.items()}
    removed = {(ds, tool): v for ds, _, _, _, _, _, rm in res for tool, v in rm.items()}
    old = [dict(set=ds[:2], seed=ds[3], tool=tool, cls=c, n_truth=v[0][c][0], tp_truth=v[0][c][1])
           for (ds, tool), v in cmp.items() for c in SN.CLS]
    bad += compare(old, ref_seed, ("set", "seed", "tool", "cls"), ("n_truth", "tp_truth"), "old recall basis")
    bad += compare([r for r in summ if r["system"] == "nesting"], read_tsv(f"{HERE}/sim_class_prf_summary.tsv"),
                   ("set", "tool", "cls"), [f"{x}_{y}" for x in METRICS for y in ("mean", "sd", "n_genomes")] +
                   ["n_truth", "n_calls", "role"], "sim_class_prf_summary.tsv")
    ct, ct0 = defaultdict(Counter), defaultdict(Counter)
    for ds, _, _, c, c0, _, _ in res:
        ct[ds[:2]].update(c)
        ct0[ds[:2]].update(c0)
    for st in ct:
        for system, j in (("nesting", 0), ("protein", 1)):
            for c in SYSTEMS[system]:
                n = sum(v for kk, v in ct[st].items() if kk[j] == c)
                ref_n = sum(r["n_truth"] for r in rows if r["set"] == st and r["tool"] == TOOLS[0] and
                            r["system"] == system and r["cls"] == c)
                if n != ref_n:
                    bad.append(f"cross-tab {st} {system} {c}: {n} vs n_truth {ref_n}")
    if bad:
        die("checks failed:\n  " + "\n  ".join(bad[:50]) + (f"\n  ... {len(bad) - 50} more" if len(bad) > 50 else ""))
    print(f"checks passed: {len(SETS) * len(TOOLS + DIAG)} genome x tool runs; ext: both systems sum to the ext ro 0.9 "
          "row, nesting rows and summary = sim_class_prf*.tsv, protein precision and old recall = sim_nonauto_prf_ext.tsv; "
          "hasltr: both systems sum to the ext row minus the removed elements and calls", file=sys.stderr)
    if a.basis == "ext":
        for st in sorted(ct0):
            print(f"{st}: {sum(ct0[st].values())} of {sum(ct[st].values())} ext elements come from library entries "
                  f"without an ~LTRlen tag (class on the whole entry): " +
                  ", ".join(f"{c} {sum(v for kk, v in ct0[st].items() if kk[1] == c)}" for c in SN.CLS), file=sys.stderr)
    else:
        print("hasltr: removed elements (copies of entries without ~LTRlen, any tool) per genome: " +
              ", ".join(f"{d} {removed[(d, TOOLS[0])][0]}" for d in SETS), file=sys.stderr)
        print("hasltr: removed calls (matched to a removed element) per genome " + " ".join(SETS) + ":", file=sys.stderr)
        for tool in TOOLS + DIAG:
            print(f"  {tool}: " + " ".join(str(removed[(d, tool)][1]) for d in SETS), file=sys.stderr)
    cols = ["set", "seed", "tool", "system", "cls", "n_truth", "tp_truth", "n_calls", "tp_calls", "precision", "recall",
            "f1", "role"]
    f_rows, f_summ = OUTS[a.basis]
    SC.write(f"{a.out}/{f_rows}", rows_b[a.basis], cols)
    SC.write(f"{a.out}/{f_summ}", summ_b[a.basis],
             ["set", "tool", "system", "cls"] + [f"{x}_{y}" for x in METRICS for y in ("mean", "sd", "n_genomes")] +
             ["n_truth", "n_calls", "role"])
    if a.basis == "hasltr":
        print(f"wrote {f_rows} ({len(rows_b[a.basis])} rows), {f_summ} ({len(summ_b[a.basis])}) to {a.out}",
              file=sys.stderr)
        return
    xrows = [dict(set=st, nesting_cls=nc, protein_cls=pc, n=ct[st][(nc, pc)], n_lib_has_ltr0=ct0[st][(nc, pc)])
             for st in sorted(ct) for nc in SC.CLASSES for pc in SN.CLS]
    SC.write(f"{a.out}/sim_class_joint_crosstab.tsv", xrows, ["set", "nesting_cls", "protein_cls", "n", "n_lib_has_ltr0"])
    if a.verbose:
        for ds, _, _, c, c0, _, _ in res:
            print(ds, "nesting x protein:", dict(c), "| from entries without ~LTRlen:", dict(c0), file=sys.stderr)
        print("protein recall: old = std == 1 has_ltr elements matched against std == 1 rows (sim_nonauto_prf_ext.tsv), "
              "new = all ext elements; old, new, shift and matching (same old elements under the ext matching - old) "
              "are means over genomes; then pooled tp/n of the ext elements by origin: old_set (std == 1, tagged "
              "entry), ext_only (std == 0: rebuilt hosts, 1 tenant in S1_2), no_ltr_tag (entry without ~LTRlen)",
              file=sys.stderr)
        for st in ("S1", "S2"):
            for tool in TOOLS + DIAG:
                for c in SN.CLS:
                    v = [[cmp[(f"{st}_{s}", tool)][i][c] for s in (1, 2, 3)] for i in (0, 1)]
                    o, sm = (statistics.mean(x[1] / x[0] for x in vv) for vv in v)
                    nw = statistics.mean(r["recall"] for r in prot if r["set"] == st and r["tool"] == tool and
                                         r["cls"] == c)
                    g = {o_: [sum(cmp[(f"{st}_{s}", tool)][2][c][o_][j] for s in (1, 2, 3)) for j in (0, 1)]
                         for o_ in ORIGINS}
                    print(f"{st}\t{tool}\t{c}\told {o:.4f}\tnew {nw:.4f}\tshift {nw - o:+.4f}\tmatching {sm - o:+.4f}\t" +
                          "\t".join(f"{o_} {g[o_][1]}/{g[o_][0]} = {SC.fmt(g[o_][1] / g[o_][0] if g[o_][0] else math.nan)}"
                                    for o_ in ORIGINS), file=sys.stderr)
    print(f"wrote sim_class_prf_joint.tsv ({len(rows)} rows), sim_class_prf_joint_summary.tsv ({len(summ)}), "
          f"sim_class_joint_crosstab.tsv ({len(xrows)}) to {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
