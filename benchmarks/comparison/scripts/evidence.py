"""Shared per-base evidence tracks for calls on a real genome (used by adjudicate.py).

Curated classes come from the curated-library RepeatMasker .out (LTR wins on any overlap, as in lib-test); CDS from a
reference gene annotation with transposon-derived genes removed; LTR-RT protein domains from dante (REXdb).
"""
import bisect
import csv
import gzip
import os
import re
from collections import defaultdict

import numpy as np

NONE, LTR, DNA, RC, LINE, SINE, SAT, OTHER = range(8)
# First keyword set (TE_NOTE). Its bare 'hAT' also matches inside words ('phosphatase', 'that'): 22,546 of
# 286,264 TAIR10 CDS lines dropped as TE that way. TE_NOTE_V2 (approach 4) fixes it with word boundaries and adds
# retroviral products for vertebrate RefSeq annotations.
TE_NOTE = re.compile(r"transpos|retrotrans|Ty1|Ty3|copia|gypsy|helitron|MuDR|Mutator|hAT|CACTA|En/Spm", re.I)
TE_NOTE_V2 = re.compile(r"transpos|retrotrans|retrovir|endogenous retro|\bTy1\b|\bTy3\b|copia|gypsy|helitron|"
                        r"\bMuDR\b|\bMutator\b|\bhAT\b|CACTA|En/Spm", re.I)
_PRE = {TE_NOTE: ("transpos", "retrotrans", "ty1", "ty3", "copia", "gypsy", "helitron", "mudr", "mutator", "hat",
                  "cacta", "en/spm"),
        TE_NOTE_V2: ("transpos", "retrotrans", "retrovir", "endogenous retro", "ty1", "ty3", "copia", "gypsy",
                     "helitron", "mudr", "mutator", "hat", "cacta", "en/spm")}


def te_hit(te_note, text):
    """te_note.search(text), with a fast substring prefilter (each alternative is a literal, so it is necessary)"""
    low = text.lower()
    return any(k in low for k in _PRE[te_note]) and te_note.search(text) is not None


def code(cl):
    """curated class code of a RepeatMasker class/family string"""
    c0 = cl.split("/")[0]
    if c0 == "LTR":
        return LTR
    if c0 in ("DNA", "MITE", "DNAnona", "DNAauto"):
        return DNA
    if c0 == "RC" or "Helitron" in cl:
        return RC
    if c0 == "LINE":
        return LINE
    if c0 == "SINE":
        return SINE
    if c0 in ("Satellite", "Simple_repeat", "Low_complexity", "Cent", "Centro", "Centromeric", "satellite"):
        return SAT
    return OTHER


def rm_hits(path):
    with open(path) as f:
        for line in f:
            t = line.split()
            if len(t) > 10 and t[0].isdigit():
                yield t


def curated_classes(out_path, lens):
    """per-base curated class (uint8) from a RepeatMasker .out"""
    arr = {c: np.zeros(n, np.uint8) for c, n in lens.items()}
    for t in rm_hits(out_path):
        if t[4] not in arr:
            continue
        k = code(t[10])
        a = arr[t[4]][int(t[5]) - 1:int(t[6])]
        if k == LTR:
            a[:] = LTR                                  # lib-test: LTR on any overlap
        else:
            a[a == NONE] = k
    return arr


def gene_gff(gdir):
    """the one gene annotation file directly in data/<ds>/genes (None if absent)"""
    if not os.path.isdir(gdir):
        return None
    hits = sorted(os.path.join(gdir, f) for f in os.listdir(gdir) if re.search(r"\.gff3?(\.gz)?$", f))
    return hits[0] if hits else None


def cds_mask(gff_path, names_tsv, lens, te_note=TE_NOTE):
    """non-TE protein-coding CDS as per-base bool arrays, renamed to the benchmark's sequence ids"""
    ren = {r["original_id"]: r["new_id"] for r in csv.DictReader(open(names_tsv), delimiter="\t")}
    arr = {c: np.zeros(n, bool) for c, n in lens.items()}
    op = gzip.open if gff_path.endswith(".gz") else open
    te_genes = set()
    rows = []
    with op(gff_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            t = line.rstrip("\n").split("\t")
            if len(t) < 9:
                continue
            if t[2] == "gene" and te_hit(te_note, t[8]):
                m = re.search(r"ID=([^;]+)", t[8])
                if m:
                    te_genes.add(m.group(1).replace("gene-", ""))
            elif t[2] == "CDS":
                rows.append(t)
    for t in rows:
        if te_hit(te_note, t[8]):
            continue
        m = re.search(r"Parent=([^;]+)", t[8])
        gene = m.group(1).split(".")[0] if m else ""       # MSU: LOC_Os01g01010.1 -> LOC_Os01g01010
        if gene in te_genes:
            continue
        c = ren.get(t[0])
        if c in arr:
            arr[c][int(t[3]) - 1:int(t[4])] = True
    return arr


def dante_domains(path):
    """REXdb Class_I|LTR protein domains: {chrom: (start0, end, strand, name, superfamily) array lists}"""
    out = defaultdict(list)
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            t = line.rstrip("\n").split("\t")
            if len(t) < 9 or t[2] != "protein_domain":
                continue
            m = re.search(r"Final_Classification=([^;]+)", t[8])
            if not m or not m.group(1).startswith("Class_I|LTR"):
                continue
            parts = m.group(1).split("|")
            sf = parts[2] if len(parts) > 2 else "NA"
            name = re.search(r"Name=([^;]+)", t[8]).group(1)
            out[t[0]].append((int(t[3]) - 1, int(t[4]), t[6], name, sf))
    res = {}
    for c, v in out.items():
        v.sort()
        res[c] = ([x[0] for x in v], v)
    return res


def coding_unit(doms, chrom, s, e):
    """True if >= 2 distinct domain names of one (strand, superfamily) lie fully inside [s, e)"""
    d = doms.get(chrom)
    if not d or e <= s:
        return False
    starts, recs = d
    i = bisect.bisect_left(starts, s)
    groups = defaultdict(set)
    while i < len(recs) and recs[i][0] < e:
        a, b, st, name, sf = recs[i]
        if b <= e:
            groups[(st, sf)].add(name)
        i += 1
    return any(len(v) >= 2 for v in groups.values())
