#!/usr/bin/env python3
"""Classify PrinTE library LTR-RT entries by REXdb protein-domain content (DANTE 0.2.11, Viridiplantae_v4.0).

Rule (the coding-unit rule of approach 4): coding = >= 2 distinct Class_I|LTR domain names of one
(strand, superfamily) inside the internal region [LTRlen, len - LTRlen); protein-less = no Class_I|LTR domain record
anywhere in the entry (any score); partial = everything else. Entries without an LTRlen tag (internal- or LTR-only
pieces) are classified on the whole entry and flagged has_ltr = 0.
Writes lib_class.tsv: family, header, len, ltr_len, int_len, n_dom, dom_names, cls.
"""
import os
import re
import sys
sys.path.insert(0, os.path.join(os.environ["BENCH_ROOT"], "scripts"))
from evidence import dante_domains, coding_unit  # noqa: E402

d = os.path.dirname(os.path.abspath(__file__)) + "/lib"
names = dict(l.rstrip("\n").split("\t") for l in open(f"{d}/names.tsv"))
lens, cur = {}, None
for l in open(f"{d}/ltr_entries.fa"):
    if l.startswith(">"):
        cur = l[1:].strip(); lens[cur] = 0
    else:
        lens[cur] += len(l.strip())
doms = dante_domains(f"{d}/dante.gff3")
with open(f"{d}/lib_class.tsv", "w") as f:
    f.write("family\theader\tlen\tltr_len\tint_len\thas_ltr\tn_dom\tdom_names\tcls\n")
    for sid, h in names.items():
        L = lens[sid]
        m = re.search(r"~LTRlen:(\d+)", h)
        ll = int(m.group(1)) if m else 0
        s, e = (ll, L - ll) if ll else (0, L)
        recs = doms.get(sid, ([], []))[1]
        nd = len(recs)
        cod = coding_unit(doms, sid, s, e)
        cls = "coding" if cod else ("proteinless" if nd == 0 else "partial")
        fam = h.split("#")[0]
        f.write(f"{fam}\t{h}\t{L}\t{ll}\t{e - s}\t{int(bool(ll))}\t{nd}\t{','.join(sorted({r[3] for r in recs}))}\t{cls}\n")
