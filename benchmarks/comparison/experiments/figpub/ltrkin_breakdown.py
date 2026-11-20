#!/usr/bin/env python3
"""ltrkin genotyping outcomes of every caller's call at the pooled real insertion polymorphisms (C. amara), broken down
by the locus's nesting and protein content.

Inputs (read-only, ltrkin_geno/): pool/<ds>/pool.tsv (0-based half-open); per-locus, per-caller outcomes from
scripts/score.py score_dataset(LG, ds, <ref>.fai, running) exactly as figures/fig1_ltrkin_yield/fig1_ltrkin_yield.py
calls it (running = rows of report/tables/ltrkin_yield.tsv with status "DNF (running)"); calls via lg_io.final_calls_path
+ read_calls (0-based half-open); DANTE domains runs/<ds>/dante_ltr/dante.gff3 (GFF3, 1-based closed, sequence IDs of
data/<ds>/genome.fa: c0001..) -> evidence.dante_domains (start - 1, end: 0-based half-open) -> original IDs via
data/<ds>/names.tsv (new_id -> original_id; lengths checked against the reference .fai).
Locus span = pool match_start..match_end (the span score.py matches calls against).
Outcome of a caller at a locus (fig1 cat_of): missed; partial (overlapped, never matched: wrong extent); correct
(found, genotyped F/E); wrong (found, any other ltrkin result: F/F, F/N, F/S, no call, other, not genotypable).
A caller whose run did not finish has status "DNF (<status>)" and no outcome.
Protein class (nonauto_sim.call_class on the locus span): coding = >= 2
distinct Class_I|LTR domain names of one (strand, superfamily) inside the middle 80% of the span; proteinless = no
Class_I|LTR domain record fully inside the span; partial = the rest.
Nesting (caller consensus): calls of the primary callers that finished, on that haplotype, linked when their reciprocal
overlap is >= 0.9 (score_elements-style: overlap / each length), clusters by union-find; clusters with calls from >= 2
callers are consensus elements, span = median start, median end of the members. A consensus element "matches" the
locus when its reciprocal overlap with the locus span is >= 0.9. tenant = the locus lies inside a consensus element at
least 1 kb longer than the locus; host = the locus contains a consensus element that does not match it, or another
pooled locus; both -> host; neither -> unnested.
CIs: score.block_bootstrap (5 Mb blocks of the reference, 1,000 draws, seed 1, the same draws as the yield table),
share = loci of the class in the outcome / loci of the class; classes with < 20 loci get counts only.
Check (fatal): over all loci, every caller's counts of the 9 ltrkin classes, the yield and its CI reproduce
report/tables/ltrkin_yield.tsv (e2e rows) to the printed digits; the correct-genotype split (LTR_peel-only loci + the
rest) sums to the yield count.
Outputs (--out): ltrkin_locus_class.tsv, ltrkin_breakdown.tsv, ltrkin_outcomes_all.tsv.
"""
import argparse
import csv
import os
import statistics
import sys
from collections import Counter, defaultdict

sys.dont_write_bytecode = True
B = os.environ["BENCH_ROOT"]
LG = f"{B}/ltrkin_geno"                  # data root of the genotyping benchmark
LGM = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "ltrkin_geno")   # its modules
sys.path.insert(0, f"{LGM}/scripts")
sys.path.insert(0, f"{LGM}/report/checks")
import score  # noqa: E402  (ltrkin modules first: nonauto_sim later puts bench/scripts at the front of sys.path)
import lg_config  # noqa: E402
from lg_io import final_calls_path, read_calls  # noqa: E402
from lead_checks import nominators  # noqa: E402
from score_elements import Index  # noqa: E402
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sim_nonauto_prf import rename_doms  # noqa: E402
import nonauto_sim  # noqa: E402
from evidence import dante_domains  # noqa: E402

DS = ["camara_hap2", "camara_hap1"]
PRIMARY = list(lg_config.TOOLS)
CATS = ["correct", "wrong", "partial", "missed"]
PROT = ["coding", "partial", "proteinless"]
NEST = ["host", "tenant", "unnested"]
RO = 0.9
MIN_LONGER = 1000       # bp: an enclosing consensus element must be at least this much longer than the locus
MIN_CI = 20             # classes with fewer loci: counts only


def die(msg):
    sys.exit(f"error: {msg}")


def cat_of(st, o):
    """fig1_ltrkin_yield.cat_of: missed / partial from the match status, else correct (F/E) or wrong."""
    if st[0] in ("missed", "partial"):
        return st[0]
    return "correct" if o[0] == "F/E" else "wrong"


def num(x):
    """Coordinate or length that may be a median of two integers: '12' or '12.5'."""
    return str(int(x)) if float(x).is_integer() else f"{x:.1f}"


def ro(a, b):
    ov = min(a["end"], b["end"]) - max(a["start"], b["start"])
    if ov <= 0:
        return 0.0
    return min(ov / (a["end"] - a["start"]), ov / (b["end"] - b["start"]))


def consensus_elements(calls, ro_min=RO, min_callers=2):
    """calls: dicts with chrom, start, end, tool. Union-find over pairs with reciprocal overlap >= ro_min; clusters
    with calls from >= min_callers tools -> dict(chrom, start, end (medians of members), callers, n_calls)."""
    par = list(range(len(calls)))

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x

    idx = Index(calls)
    for i, c in enumerate(calls):
        for j in idx.overlaps(c["chrom"], c["start"], c["end"]):
            if j > i and ro(c, calls[j]) >= ro_min:
                a, b = find(i), find(j)
                if a != b:
                    par[a] = b
    cl = defaultdict(list)
    for i in range(len(calls)):
        cl[find(i)].append(calls[i])
    out = []
    for mem in cl.values():
        tools = sorted({m["tool"] for m in mem})
        if len(tools) < min_callers:
            continue
        out.append(dict(chrom=mem[0]["chrom"], start=statistics.median(m["start"] for m in mem),
                        end=statistics.median(m["end"] for m in mem), callers=",".join(tools), n_calls=len(mem)))
    out.sort(key=lambda e: (e["chrom"], e["start"], e["end"]))
    return out


def nesting(loci, cons, ro_min=RO, min_longer=MIN_LONGER):
    """loci: dicts with locus_id, chrom, start, end. Returns {locus_id: dict(cls, evidence...)}."""
    cidx, lidx = Index(cons), Index(loci)
    out = {}
    for L in loci:
        ll = L["end"] - L["start"]
        enc, inner, selfm = [], [], []
        for k in cidx.overlaps(L["chrom"], L["start"], L["end"]):
            e = cons[k]
            el = e["end"] - e["start"]
            r = ro(L, e)
            if r >= ro_min:
                selfm.append(e)
            if e["start"] <= L["start"] and L["end"] <= e["end"] and el - ll >= min_longer:
                enc.append(e)
            if L["start"] <= e["start"] and e["end"] <= L["end"] and r < ro_min:
                inner.append((e, r))
        inner_loci = [loci[k]["locus_id"] for k in lidx.overlaps(L["chrom"], L["start"], L["end"])
                      if loci[k]["locus_id"] != L["locus_id"] and L["start"] <= loci[k]["start"]
                      and loci[k]["end"] <= L["end"]]
        host = bool(inner or inner_loci)
        cls = "host" if host else ("tenant" if enc else "unnested")
        sp = lambda e: f"{e['chrom']}:{num(e['start'])}-{num(e['end'])}"
        e0 = min(enc, key=lambda e: e["end"] - e["start"]) if enc else None      # closest enclosing element
        i0 = max(inner, key=lambda x: x[0]["end"] - x[0]["start"]) if inner else None   # largest contained element
        out[L["locus_id"]] = dict(
            nesting=cls, is_tenant=int(bool(enc)), is_host=int(host), n_self_cons=len(selfm),
            enclosing=sp(e0) if e0 else ".", enclosing_len=num(e0["end"] - e0["start"]) if e0 else ".",
            enclosing_callers=e0["callers"] if e0 else ".", n_enclosing=len(enc),
            contained=sp(i0[0]) if i0 else ".", contained_len=num(i0[0]["end"] - i0[0]["start"]) if i0 else ".",
            contained_ro=f"{i0[1]:.3f}" if i0 else ".", contained_callers=i0[0]["callers"] if i0 else ".",
            n_contained=len(inner), contained_loci=",".join(inner_loci) or ".")
    return out


def load_domains(gff, names_tsv, fai):
    """DANTE GFF3 (1-based closed, prep_genome IDs) -> evidence.dante_domains (0-based half-open) keyed by original
    IDs; each mapped sequence's length in names.tsv must equal the reference .fai."""
    doms = rename_doms(dante_domains(gff), names_tsv)
    flen = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(fai)}
    for r in csv.DictReader(open(names_tsv), delimiter="\t"):
        if flen.get(r["original_id"]) != int(r["length"]):
            die(f"{names_tsv}: {r['original_id']} length {r['length']} != {fai} {flen.get(r['original_id'])}")
    return doms


def domain_evidence(doms, chrom, s, e):
    """Class_I|LTR domain records fully inside [s, e) (the records call_class counts): n and sorted names."""
    d = doms.get(chrom)
    if not d:
        return 0, "."
    starts, recs = d
    import bisect
    i = bisect.bisect_left(starts, s - 20000)
    names = []
    while i < len(recs) and recs[i][0] < e:
        if recs[i][0] >= s and recs[i][1] <= e:
            names.append(f"{recs[i][3]}:{recs[i][4]}:{recs[i][2]}")
        i += 1
    return len(names), ",".join(sorted(set(names))) or "."


def fmt6(v):
    return score.fmt(v)


def shares_ci(cats_by_tool, idx, ub, blocks, tools):
    """Per tool and category: count, share and block-bootstrap CI over the loci idx (CI only if len(idx) >= MIN_CI)."""
    out = {t: {} for t in tools}
    n = len(idx)
    for k in CATS:
        if n >= MIN_CI:
            units = [(ub[i], {t: (int(cats_by_tool[t][i] == k), 1) for t in tools}) for i in idx]
            ci = score.block_bootstrap(units, blocks, tools)
        for t in tools:
            c = sum(cats_by_tool[t][i] == k for i in idx)
            out[t][k] = dict(n=c, share=c / n if n else None,
                             lo=ci[t]["lo"] if n >= MIN_CI else None, hi=ci[t]["hi"] if n >= MIN_CI else None)
    return out


def write(path, rows, cols):
    with open(path, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join("." if r.get(c) is None else (f"{r[c]:.6g}" if isinstance(r[c], float) else str(r[c]))
                              for c in cols) + "\n")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", default=os.path.dirname(os.path.abspath(__file__)), help="output dir (default: script dir)")
    ap.add_argument("-v", "--verbose", action="store_true", help="per-dataset counts and checks on stderr")
    a = ap.parse_args(argv)
    yt = list(csv.DictReader(open(f"{LG}/report/tables/ltrkin_yield.tsv"), delimiter="\t"))
    Y = {(r["ds"], r["tool"], r["run"]): r for r in yt}
    running = {(r["ds"], r["tool"]) for r in yt if r["status"] == "DNF (running)"}
    loc_rows, brk_rows, all_rows, bad = [], [], [], []
    for ds in DS:
        fai = str(lg_config.DATASETS[ds]["ref"]) + ".fai"
        pool, _, res, ok, N = score.score_dataset(LG, ds, fai, running)
        blocks = score.read_blocks(fai)
        ub = [(p["chrom"], p["match_start"] // lg_config.BLOCK) for p in pool]
        cats = {t: [cat_of(s, o) for s, o in zip(res[t]["status_rows"], res[t]["outcome"]["e2e"])] for t in ok}
        fine = {t: [s[0] if s[0] != "found" else o[0] for s, o in zip(res[t]["status_rows"], res[t]["outcome"]["e2e"])]
                for t in ok}
        # protein class
        doms = load_domains(f"{LG}/runs/{ds}/dante_ltr/dante.gff3", f"{LG}/data/{ds}/names.tsv", fai)
        loci = [dict(locus_id=p["locus_id"], chrom=p["chrom"], start=p["match_start"], end=p["match_end"]) for p in pool]
        pcls = [nonauto_sim.call_class(doms, L) for L in loci]
        dev = [domain_evidence(doms, L["chrom"], L["start"], L["end"]) for L in loci]
        for c, (nd, _) in zip(pcls, dev):
            if (c == "proteinless") != (nd == 0 and c != "coding"):
                die(f"{ds}: domain evidence disagrees with call_class ({c}, {nd} domains)")
        dl = read_calls(final_calls_path(f"{LG}/calls/{ds}", "dante_ltr"), "dante_ltr") if "dante_ltr" in ok else []
        dlc = Counter(nonauto_sim.call_class(doms, c) for c in dl)
        # nesting from caller consensus
        calls = []
        for t in PRIMARY:
            if t in ok:
                calls += read_calls(final_calls_path(f"{LG}/calls/{ds}", t), t)
        cons = consensus_elements(calls)
        nest = nesting(loci, cons)
        lq_only = [nominators(p) == {"ltr_peel"} for p in pool]
        if a.verbose:
            print(f"{ds}: {N} loci; callers ok {','.join(ok)}; {len(calls)} calls -> {len(cons)} consensus elements; "
                  f"DANTE_LTR calls by class {dict(dlc)}; protein {dict(Counter(pcls))}; nesting "
                  f"{dict(Counter(v['nesting'] for v in nest.values()))}; LTR_peel-only loci {sum(lq_only)}",
                  file=sys.stderr)
        for i, (p, L) in enumerate(zip(pool, loci)):
            nv = nest[L["locus_id"]]
            loc_rows.append(dict(ds=ds, locus_id=L["locus_id"], chrom=L["chrom"], match_start=L["start"],
                                 match_end=L["end"], length=L["end"] - L["start"], basis=p["basis"],
                                 nominators=",".join(sorted(nominators(p))), ltr_peel_only=int(lq_only[i]),
                                 protein_class=pcls[i], n_domains=dev[i][0], domains=dev[i][1], **nv,
                                 **{f"out_{t}": cats[t][i] for t in PRIMARY if t in ok}))
        # check: all-loci fine classes, yield and CI vs ltrkin_yield.tsv (e2e)
        for t in ok:
            y = Y[(ds, t, "e2e")]
            fc = Counter(fine[t])
            for k in score.CLASSES:
                if fmt6(fc[k] / N) != y[f"share_{k}"] and not (fc[k] == 0 and y[f"share_{k}"] == "0"):
                    bad.append(f"{ds} {t} share_{k}: {fmt6(fc[k] / N)} vs {y[f'share_{k}']}")
            if int(y["n_loci"]) != N:
                bad.append(f"{ds} {t}: n_loci {N} vs {y['n_loci']}")
        tools_all = list(ok)
        allc = shares_ci(cats, list(range(N)), ub, blocks, tools_all)
        for t in lg_config.ALL_TOOLS:
            role = "secondary" if t in lg_config.SECONDARY else "primary"
            base = dict(ds=ds, tool=t, role=role, n_loci=N)
            if t not in ok:
                st = res[t]["status"]
                all_rows.append(dict(base, status=f"DNF ({st})"))
                for dim, labels, vals in (("protein", PROT, pcls), ("nesting", NEST, [nest[L["locus_id"]]["nesting"]
                                                                                       for L in loci])):
                    for c in labels:
                        brk_rows.append(dict(base, status=f"DNF ({st})", dimension=dim, cls=c,
                                             n_loci=sum(v == c for v in vals)))
                continue
            y = Y[(ds, t, "e2e")]
            ac = allc[t]
            if fmt6(ac["correct"]["share"]) != y["yield"] or fmt6(ac["correct"]["lo"]) != y["yield_lo"] or \
                    fmt6(ac["correct"]["hi"]) != y["yield_hi"]:
                bad.append(f"{ds} {t}: correct {fmt6(ac['correct']['share'])} [{fmt6(ac['correct']['lo'])}, "
                           f"{fmt6(ac['correct']['hi'])}] vs yield {y['yield']} [{y['yield_lo']}, {y['yield_hi']}]")
            n_lq = sum(cats[t][i] == "correct" and lq_only[i] for i in range(N))
            n_rest = sum(cats[t][i] == "correct" and not lq_only[i] for i in range(N))
            if n_lq + n_rest != ac["correct"]["n"]:
                bad.append(f"{ds} {t}: split {n_lq}+{n_rest} != {ac['correct']['n']}")
            row = dict(base, status="ok")
            for k in CATS:
                row.update({f"n_{k}": ac[k]["n"], f"share_{k}": ac[k]["share"], f"{k}_lo": ac[k]["lo"],
                            f"{k}_hi": ac[k]["hi"]})
            row.update(n_ltr_peel_only_loci=sum(lq_only), n_correct_ltr_peel_only=n_lq, n_correct_rest=n_rest,
                       share_correct_ltr_peel_only=n_lq / N, share_correct_rest=n_rest / N)
            all_rows.append(row)
        for dim, labels, vals in (("protein", PROT, pcls),
                                  ("nesting", NEST, [nest[L["locus_id"]]["nesting"] for L in loci])):
            for c in labels:
                idx = [i for i, v in enumerate(vals) if v == c]
                sc = shares_ci(cats, idx, ub, blocks, tools_all)
                for t in tools_all:
                    row = dict(ds=ds, tool=t, role="secondary" if t in lg_config.SECONDARY else "primary",
                               status="ok", dimension=dim, cls=c, n_loci=len(idx))
                    for k in CATS:
                        row.update({f"n_{k}": sc[t][k]["n"], f"share_{k}": sc[t][k]["share"],
                                    f"{k}_lo": sc[t][k]["lo"], f"{k}_hi": sc[t][k]["hi"]})
                    brk_rows.append(row)
    if bad:
        die("checks failed:\n  " + "\n  ".join(bad))
    print("checks passed: all-loci shares of the 9 ltrkin classes, n_loci, yield and its CI equal ltrkin_yield.tsv "
          "(e2e) for every finished caller; correct = LTR_peel-only + rest", file=sys.stderr)
    order = {t: i for i, t in enumerate(lg_config.ALL_TOOLS)}
    brk_rows.sort(key=lambda r: (DS.index(r["ds"]), r["dimension"] != "nesting",
                                 (NEST + PROT).index(r["cls"]), order[r["tool"]]))
    all_rows.sort(key=lambda r: (DS.index(r["ds"]), order[r["tool"]]))
    cc = [f"{x}" for k in CATS for x in (f"n_{k}", f"share_{k}", f"{k}_lo", f"{k}_hi")]
    write(f"{a.out}/ltrkin_breakdown.tsv", brk_rows, ["ds", "tool", "role", "status", "dimension", "cls", "n_loci"] + cc)
    write(f"{a.out}/ltrkin_outcomes_all.tsv", all_rows,
          ["ds", "tool", "role", "status", "n_loci"] + cc + ["n_ltr_peel_only_loci", "n_correct_ltr_peel_only",
                                                             "n_correct_rest", "share_correct_ltr_peel_only",
                                                             "share_correct_rest"])
    lc = ["ds", "locus_id", "chrom", "match_start", "match_end", "length", "basis", "nominators", "ltr_peel_only",
          "protein_class", "n_domains", "domains", "nesting", "is_tenant", "is_host", "n_self_cons", "enclosing",
          "enclosing_len", "enclosing_callers", "n_enclosing", "contained", "contained_len", "contained_ro",
          "contained_callers", "n_contained", "contained_loci"] + [f"out_{t}" for t in PRIMARY]
    write(f"{a.out}/ltrkin_locus_class.tsv", loc_rows, lc)


if __name__ == "__main__":
    main()
