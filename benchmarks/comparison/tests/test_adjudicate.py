"""Approach 4: per-call evidence (scripts/adjudicate.py, scripts/evidence.py) and summaries (scripts/adj_summary.py)."""
import csv
import os
import random
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts"))
import adj_summary as su  # noqa: E402
import adjudicate as adj  # noqa: E402
import evidence as ev  # noqa: E402

LENS = {"chrA": 400_000, "chrB": 300_000}


def rm_line(s1, e, cls, name="X", chrom="chrA"):
    return f"  300  10.0  0.0  0.0  {chrom}  {s1}  {e}  (0)  +  {name}  {cls}  1  {e - s1 + 1}  (0)  1\n"


@pytest.fixture
def toy(tmp_path):
    root, ds = tmp_path, "toy"
    rnd = random.Random(7)
    seqs = {c: "".join(rnd.choice("ACGT") for _ in range(n)) for c, n in LENS.items()}
    a = seqs["chrA"]
    # call1 10,000-20,000: TG..CA ends + TSD 'ACGTA'; chrB 100,000-260,000 is N (null must avoid it)
    a = a[:9995] + "ACGTA" + "TG" + a[10002:19998] + "CA" + "ACGTA" + a[20005:]
    seqs["chrA"] = a
    b = seqs["chrB"]
    seqs["chrB"] = b[:100_000] + "N" * 160_000 + b[260_000:]
    d = root / "data" / ds
    (d / "genes").mkdir(parents=True)
    off = 0
    with open(d / "genome.fa", "w") as f, open(d / "genome.fa.fai", "w") as fi:
        for c, s in seqs.items():
            f.write(f">{c}\n{s}\n")
            off += len(c) + 2
            fi.write(f"{c}\t{len(s)}\t{off}\t{len(s)}\t{len(s) + 1}\n")
            off += len(s) + 1
    with open(d / "names.tsv", "w") as f:
        f.write("new_id\toriginal_id\tlength\tn_count\tnonacgtn_converted\n")
        for c, n in LENS.items():
            f.write(f"{c}\t{c}\t{n}\t0\t0\n")
    # genes: real CDS inside call3; a retroviral CDS inside call1 (must be ignored)
    with open(d / "genes" / "genes.gff3", "w") as f:
        f.write("chrA\tx\tCDS\t50101\t50700\t.\t+\t0\tParent=g1.1;product=kinase\n")
        f.write("chrA\tx\tCDS\t12001\t12600\t.\t+\t0\tParent=g2.1;product=endogenous retrovirus K Pol\n")
    m = root / "mask" / ds / "curated"
    m.mkdir(parents=True)
    with open(m / f"{ds}.fa.out", "w") as f:
        f.write("   SW  perc perc perc  query  position\n\n")
        f.write(rm_line(10001, 20000, "LTR/Gypsy"))           # call1 supported
        f.write(rm_line(50001, 50500, "DNA/hAT"))             # call3 '5' LTR'
        f.write(rm_line(59501, 60000, "DNA/hAT"))             # call3 '3' LTR'
        f.write(rm_line(70001, 80000, "Satellite/cent"))      # call4
    r = root / "runs" / ds / "dante_ltr"
    r.mkdir(parents=True)
    with open(r / "dante.gff3", "w") as f:
        dom = "chrA\tdante\tprotein_domain\t{}\t{}\t99\t{}\t.\tName={};Final_Classification=Class_I|LTR|{}\n"
        f.write(dom.format(13001, 13600, "+", "RT", "Ty3/gypsy"))     # call1: RT + INT, one strand/superfamily
        f.write(dom.format(15001, 15600, "+", "INT", "Ty3/gypsy"))
        f.write(dom.format(31001, 31600, "+", "RT", "Ty3/gypsy"))     # call2: one domain only
        f.write(dom.format(91001, 91600, "+", "RT", "Ty1/copia"))     # call5: two strands
        f.write(dom.format(93001, 93600, "-", "INT", "Ty1/copia"))
    t = root / "eval" / "trf" / ds
    t.mkdir(parents=True)
    (t / "trf.bed").write_text("chrA\t70000\t80000\n")
    c = root / "calls" / ds
    c.mkdir(parents=True)
    calls = [("chrA", 10000, 20000), ("chrA", 30000, 40000), ("chrA", 50000, 60000), ("chrA", 70000, 80000),
             ("chrA", 90000, 95000), ("chrB", 90000, 99000)]
    with open(c / "t.ids.tsv", "w") as f:
        f.write("chrom\tstart\tend\tstrand\tid\tparent\tlltr_end\trltr_start\ttsd\tclass\n")
        for i, (ch, s, e) in enumerate(calls):
            f.write(f"{ch}\t{s}\t{e}\t+\tcall{i + 1}\t.\t.\t.\t.\tLTR\n")
    k = root / "eval" / "k2p" / ds / "t"
    k.mkdir(parents=True)
    with open(k / "annot.tsv", "w") as f:
        f.write("seq_id\tseq_len\tstatus\tltr5_start\tltr5_end\tltr3_start\tltr3_end\tk2p\n")
        f.write("chrA:10001-20000\t10000\tpass\t1\t500\t9501\t10000\t0.02\n")
        f.write("chrA:50001-60000\t10000\tpass\t1\t500\t9501\t10000\t0.15\n")
    return root, ds


def run(root, ds, tmp_path, seed=1):
    out = tmp_path / f"out{seed}.tsv"
    adj.main(["--ds", ds, "--tool", "t", "--root", str(root), "--out", str(out), "--seed", str(seed)])
    return {r["id"]: r for r in csv.DictReader(open(out), delimiter="\t")}


def test_flags(toy, tmp_path):
    root, ds = toy
    r = run(root, ds, tmp_path)
    c1 = r["call1"]
    assert (c1["supported"], c1["coding"], c1["gene"], c1["junk"]) == ("1", "1", "0", "0")   # retroviral CDS ignored
    assert (c1["age"], c1["inner"], c1["tgca"], c1["tsd"]) == ("0.01-0.03", "kmer2ltr", "1", "1")
    assert r["call2"]["coding"] == "0" and r["call2"]["supported"] == "0"                     # one domain only
    c3 = r["call3"]
    assert (c3["gene"], c3["nonltr"], c3["junk"], c3["age"]) == ("1", "1", "1", ">=0.1")
    c4 = r["call4"]
    assert (c4["sat"], c4["trf"], c4["junk"]) == ("1", "1", "1")
    assert r["call5"]["coding"] == "0"                                                         # strands differ
    assert r["call2"]["age"] == "NA" and r["call2"]["inner"] == "middle80"


def test_null(toy, tmp_path):
    root, ds = toy
    r1, r2 = run(root, ds, tmp_path, 1), run(root, ds, tmp_path, 1)
    assert [x["null_shift"] for x in r1.values()] == [x["null_shift"] for x in r2.values()]   # seeded
    for x in r1.values():
        if x["null_shift"] == "NA":
            continue
        d, s, e = int(x["null_shift"]), int(x["start"]), int(x["end"])
        assert 20_000 <= abs(d) <= 200_000
        assert s + d >= e or e + d <= s                                                        # no overlap
        assert 0 <= s + d and e + d <= LENS[x["chrom"]]
    b = r1["call6"]                                       # chrB 90-99 kb: shifts into 100-260 kb N are rejected
    if b["null_shift"] != "NA":
        d = int(b["null_shift"])
        assert 90_000 + d >= 259_100 or 99_000 + d <= 100_900


def test_coding_unit_superfamily(tmp_path):
    p = tmp_path / "d.gff3"
    dom = "c\tdante\tprotein_domain\t{}\t{}\t9\t+\t.\tName={};Final_Classification=Class_I|LTR|{}\n"
    p.write_text(dom.format(101, 200, "RT", "Ty1/copia") + dom.format(301, 400, "INT", "Ty3/gypsy") +
                 dom.format(501, 600, "RH", "Ty3/gypsy") + dom.format(701, 800, "GAG", "Ty3/gypsy"))
    d = ev.dante_domains(str(p))
    assert not ev.coding_unit(d, "c", 100, 450)       # RT copia + INT gypsy: different superfamilies
    assert ev.coding_unit(d, "c", 300, 600)           # INT + RH gypsy
    assert not ev.coding_unit(d, "c", 300, 550)       # RH not fully inside
    assert ev.coding_unit(d, "c", 0, 10_000)          # a span holding a whole nested element counts (null controls)


def test_cds_keywords(tmp_path):
    g = tmp_path / "g.gff3"
    g.write_text("s\tx\tCDS\t1\t400\t.\t+\t0\tParent=a.1;product=endogenous retrovirus group K Gag\n"
                 "s\tx\tCDS\t501\t900\t.\t+\t0\tParent=b.1;product=kinase\n")
    n = tmp_path / "n.tsv"
    n.write_text("new_id\toriginal_id\nc1\ts\n")
    m1 = ev.cds_mask(str(g), str(n), {"c1": 1000}, ev.TE_NOTE)
    m2 = ev.cds_mask(str(g), str(n), {"c1": 1000}, ev.TE_NOTE_V2)
    assert m1["c1"].sum() == 800 and m2["c1"].sum() == 400


def rows(spec):
    return [dict(chrom=c, start=s, end=e, supported=sup, coding=cod) for c, s, e, sup, cod in spec]


def test_pool_recall():
    a = rows([("c", 0, 10_000, "1", "0"), ("c", 50_000, 60_000, "1", "0"), ("c", 90_000, 91_000, "0", "0")])
    b = rows([("c", 200, 10_100, "1", "0"), ("c", 120_000, 130_000, "1", "0")])
    host = rows([("c", 40_000, 80_000, "1", "0")])        # host around a's 50-60 kb tenant: RO 0.25, separate locus
    res = su.element_library({"a": a, "b": b, "h": host}, ["a", "b", "h"], n_boot=20)
    assert res["a"]["n_loci_lib"] == 4
    assert res["a"]["recall_lib"] == pytest.approx(2 / 4) and res["b"]["recall_lib"] == pytest.approx(2 / 4)
    assert res["h"]["recall_lib"] == pytest.approx(1 / 4)
    assert res["a"]["precision_lib"] == pytest.approx(2 / 3)


def frows(n_sup, p_sup, n_uns, p_uns, p_null, seed=0):
    rng = random.Random(seed)
    out = []
    for _ in range(n_sup):
        out.append(dict(age="<0.01", supported="1", coding=str(int(rng.random() < p_sup)), null_coding="0"))
    for _ in range(n_uns):
        out.append(dict(age="<0.01", supported="0", coding=str(int(rng.random() < p_uns)),
                        null_coding=str(int(rng.random() < p_null))))
    return out


def test_f_real():
    f0, _, _ = su.f_real(frows(2000, 0.8, 2000, 0.1, 0.1), n_boot=50)
    f1, ci, _ = su.f_real(frows(2000, 0.8, 2000, 0.8, 0.1), n_boot=50)
    assert abs(f0) < 0.1 and abs(f1 - 1) < 0.1 and ci[0] < f1 < ci[1]


def test_wilson_logistic():
    lo, hi = su.wilson(50, 100)
    assert lo < 0.5 < hi
    x = np.linspace(-2, 0, 400)
    y = (np.random.default_rng(0).random(400) < 1 / (1 + np.exp(-(1 - 3 * (x + 1))))).astype(float)
    slope, ci = su.logistic(x, y)
    assert ci[0] < slope < ci[1] < 0
