#!/usr/bin/env python3
"""Ground truth for LTR-RT callers from a PrinTE BED (0-based, half-open), including nested elements.

PrinTE records an insertion into an existing feature by splitting the host into pieces tagged ';CUT_BY:<tenant>' and
tagging the tenant ';NESTED_IN:<host>'. A host's pieces share feature ID, TSD and strand and are separated only by its
tenants (and the tenants' TSD copies). This script rebuilds every intact element from its pieces and classifies LTR-RTs:

  intact       uncut, not nested                       (in PrinTE's standard --printe-intact set)
  tenant       uncut, nested inside another feature     (in the standard set)
  host         rebuilt from >= 2 pieces, not nested     (NOT in the standard set: PrinTE calls host pieces fragments)
  host_tenant  rebuilt host that is itself nested

Outputs (prefix OUT):
  OUT.truth.tsv     one row per intact LTR-RT (columns below), span = rebuilt genome span
  OUT.pairs.tsv     host_id, tenant_id for every element directly nested in an LTR-RT host (tenant may be any TE)
  OUT.features.bed  every BED record with a category (gene, LTR_frag, LTR_solo, nonLTR_TE, nonLTR_frag, LTR_piece,
                    LTR_intact) for false-positive anatomy
"""
import argparse
import bisect
import re
import sys
from collections import defaultdict

LTRLEN = re.compile(r"~LTRlen:(\d+)")
MAX_GAP = 30          # bp of a host's span not tiled by records (a tenant's TSD copy is 2-20 bp)


def parse_fid(fid):
    """'name#Class/Super~LTRlen:N[_SOLO]' -> (family, cls, superfamily, ltr_len, kind)."""
    if fid.startswith("gene"):
        return fid, "gene", "", 0, "gene"
    fam, _, rest = fid.partition("#")
    cls, _, sf = rest.partition("/")
    m = LTRLEN.search(fid)
    sf = sf.split("~")[0]
    kind = "frag" if "_FRAG" in fid else "solo" if "_SOLO" in fid else "intact"
    return fam, cls, sf, int(m.group(1)) if m else 0, kind


def read_bed(path):
    recs = []
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 6:
                continue
            name = c[3]
            fid, *attrs = name.split(";")
            fam, cls, sf, ll, kind = parse_fid(fid)
            recs.append(dict(chrom=c[0], start=int(c[1]), end=int(c[2]), name=name, fid=fid, tsd=c[4],
                             strand=c[5], attrs=attrs, fam=fam, cls=cls, sf=sf, ltr_len=ll, kind=kind,
                             cut=any(a.startswith("CUT_BY:") for a in attrs),
                             nested=any(a.startswith("NESTED_IN:") for a in attrs)))
    return recs


class Coverage:
    """Per-chromosome record intervals for 'is [a, b) tiled by records with gaps <= MAX_GAP?'."""

    def __init__(self, recs):
        self.iv = defaultdict(list)
        for r in recs:
            self.iv[r["chrom"]].append((r["start"], r["end"]))
        for v in self.iv.values():
            v.sort()
        self.starts = {k: [s for s, _ in v] for k, v in self.iv.items()}

    def tiled(self, chrom, a, b, max_gap=MAX_GAP):
        if b - a <= max_gap:
            return True
        iv, st = self.iv[chrom], self.starts[chrom]
        i = max(0, bisect.bisect_right(st, a) - 1)
        pos = a
        while i < len(iv) and iv[i][0] < b:
            s, e = iv[i]
            if e > pos:
                if s - pos > max_gap:
                    return False
                pos = max(pos, e)
            i += 1
        return b - pos <= max_gap


def printe_standard(recs):
    """PrinTE's own intact-LTR-RT set (util/bedtools.py parse_bed_intact, ported verbatim in logic): first supplemental
    attribute containing CUT_BY -> fragment; rule 3 (+-100 records, same TSD and strand, prefix-related name) ->
    fragment; keep class LTR. Returns {(chrom, start, end)}."""
    cat = []
    for r in recs:
        add = r["attrs"]
        if r["fid"].startswith("gene"):
            c = "gene"
        elif "_FRAG" in r["fid"] or (add and "CUT_BY" in add[0]):
            c = "frag"
        elif "_SOLO" in r["fid"]:
            c = "solo"
        else:
            c = "pot"
        cat.append(c)
    for i, r in enumerate(recs):
        if cat[i] == "pot" and r["attrs"]:
            for j in range(max(0, i - 100), min(len(recs), i + 101)):
                o = recs[j]
                if i == j or o["fid"].startswith("gene"):
                    continue
                if o["tsd"] == r["tsd"] and o["strand"] == r["strand"] and \
                        (r["name"].startswith(o["name"]) or o["name"].startswith(r["name"])):
                    cat[i] = "frag"
                    break
    return {(r["chrom"], r["start"], r["end"]) for r, c in zip(recs, cat) if c == "pot" and r["cls"] == "LTR"}


def build(recs):
    """Group records into element instances. Pieces of one element share feature ID, TSD and strand and the space
    between them is tiled by other records (its tenants). Cut status comes from structure (>= 2 pieces), never from
    tags: a tenant's NESTED_IN value is its host's full name, which may itself carry CUT_BY."""
    cov = Coverage(recs)
    groups = defaultdict(list)
    singles = []
    for r in recs:
        if r["kind"] in ("intact", "gene"):
            groups[(r["chrom"], r["fid"], r["tsd"], r["strand"])].append(r)
        else:
            singles.append([r])
    inst = singles
    for key, pieces in groups.items():
        pieces.sort(key=lambda r: r["start"])
        cur = [pieces[0]]
        for p in pieces[1:]:
            if cov.tiled(key[0], cur[-1]["end"], p["start"]):
                cur.append(p)
            else:
                inst.append(cur)
                cur = [p]
        inst.append(cur)
    out = []
    for ps in inst:
        r0 = ps[0]
        for p in ps:
            p["piece"] = len(ps) > 1
        lens = [p["end"] - p["start"] for p in ps]
        cum = [sum(lens[:i + 1]) for i in range(len(lens) - 1)]        # insertion offsets in host coordinates
        out.append(dict(chrom=r0["chrom"], start=min(p["start"] for p in ps), end=max(p["end"] for p in ps),
                        strand=r0["strand"], fid=r0["fid"], fam=r0["fam"], cls=r0["cls"], sf=r0["sf"],
                        ltr_len=r0["ltr_len"], kind=r0["kind"], tsd=r0["tsd"], n_pieces=len(ps),
                        seq_len=sum(lens), cut=len(ps) > 1, nested=any(p["nested"] for p in ps),
                        ins_offsets=cum, pieces=[(p["start"], p["end"]) for p in ps]))
    return out


def classify(inst, std):
    elems = [x for x in inst if x["kind"] == "intact" and x["cls"] == "LTR"]
    for x in elems:
        host = x["n_pieces"] > 1
        x["class"] = ("host_tenant" if x["nested"] else "host") if host else ("tenant" if x["nested"] else "intact")
        x["std"] = int((x["chrom"], x["start"], x["end"]) in std)
        ll, L = x["ltr_len"], x["seq_len"]
        x["tenant_in_ltr"] = int(host and any(o < ll or o > L - ll for o in x["ins_offsets"]))
    return elems


def nest_links(inst, elems):
    """Direct host for every nested element (any intact TE inside an LTR-RT host); depth of each LTR-RT element."""
    hosts = sorted([x for x in elems if x["n_pieces"] > 1], key=lambda x: x["end"] - x["start"])
    by_chrom = defaultdict(list)
    for h in hosts:
        by_chrom[h["chrom"]].append(h)
    cand = [x for x in inst if x["kind"] == "intact" and x["nested"]]
    pairs = []
    for x in cand:
        best = None
        for h in by_chrom[x["chrom"]]:                      # smallest containing host = direct host
            if h is not x and h["start"] <= x["start"] and x["end"] <= h["end"] and \
                    (h["end"] - h["start"]) > (x["end"] - x["start"]):
                best = h
                break
        if best is not None:
            pairs.append((best, x))
            x["parent"] = best
    for x in elems:
        x["depth"] = 0
    for h in hosts:                                          # increasing size: inner hosts first
        kids = [t for hh, t in pairs if hh is h]
        h["depth"] = 1 + max((t.get("depth", 0) for t in kids), default=0)
    return pairs


def eid(x):
    return f"{x['chrom']}:{x['start']}-{x['end']}"


def main(argv=None):
    ap = argparse.ArgumentParser(description="LTR-RT truth (incl. rebuilt nested hosts) from a PrinTE BED.")
    ap.add_argument("bed", help="PrinTE BED (gen<N>_final.bed)")
    ap.add_argument("out", help="output prefix")
    ap.add_argument("-v", "--verbose", action="store_true", help="print class counts and checks")
    a = ap.parse_args(argv)
    recs = read_bed(a.bed)
    if not recs:
        sys.exit(f"error: no records in {a.bed}")
    std = printe_standard(recs)
    inst = build(recs)
    elems = classify(inst, std)
    pairs = nest_links(inst, elems)
    elems.sort(key=lambda x: (x["chrom"], x["start"]))
    with open(a.out + ".truth.tsv", "w") as f:
        f.write("id\tchrom\tstart\tend\tstrand\tclass\tstd\tfamily\tsuperfamily\tltr_len\ttsd\tn_pieces\tseq_len"
                "\tdepth\ttenant_in_ltr\tparent\tpieces\n")
        for x in elems:
            par = eid(x["parent"]) if "parent" in x else "."
            pcs = ",".join(f"{s}-{e}" for s, e in x["pieces"])
            f.write(f"{eid(x)}\t{x['chrom']}\t{x['start']}\t{x['end']}\t{x['strand']}\t{x['class']}\t{x['std']}\t"
                    f"{x['fam']}\t{x['sf']}\t{x['ltr_len']}\t{x['tsd']}\t{x['n_pieces']}\t{x['seq_len']}\t"
                    f"{x['depth']}\t{x['tenant_in_ltr']}\t{par}\t{pcs}\n")
    with open(a.out + ".pairs.tsv", "w") as f:
        f.write("host_id\ttenant_id\ttenant_class\n")
        for h, t in pairs:
            f.write(f"{eid(h)}\t{eid(t)}\t{t['cls']}\n")
    with open(a.out + ".features.bed", "w") as f:
        for r in recs:
            if r["kind"] == "gene":
                cat = "gene"
            elif r["cls"] == "LTR":
                cat = {"frag": "LTR_frag", "solo": "LTR_solo"}.get(r["kind"], "LTR_piece" if r.get("piece") else "LTR_intact")
            else:
                cat = "nonLTR_frag" if r["kind"] == "frag" else "nonLTR_TE"
            f.write(f"{r['chrom']}\t{r['start']}\t{r['end']}\t{cat}\t{r['fid']}\t{r['strand']}\n")
    with open(a.out + ".std.bed", "w") as f:            # PrinTE's standard set, exactly
        for c, s0, e0 in sorted(std):
            f.write(f"{c}\t{s0}\t{e0}\n")
    n = defaultdict(int)
    for x in elems:
        n[x["class"]] += 1
    deep = sum(1 for x in elems if x["depth"] >= 2)
    inltr = sum(x["tenant_in_ltr"] for x in elems)
    print(f"printe_truth: {len(elems)} intact LTR-RTs: " + ", ".join(f"{k} {v}" for k, v in sorted(n.items()))
          + f"; PrinTE standard set {len(std)} ({sum(x['std'] for x in elems)} are whole elements here); depth>=2 hosts {deep}; hosts with tenant in an LTR "
          f"{inltr}; nest pairs {len(pairs)}", file=sys.stderr)


if __name__ == "__main__":
    main()
