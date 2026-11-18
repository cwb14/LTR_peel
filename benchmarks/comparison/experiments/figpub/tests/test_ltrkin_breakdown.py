"""Tiny synthetic fixtures for ltrkin_breakdown: nesting rule, GFF3 1-based closed -> 0-based half-open, protein class.
Run with the ltrkin_geno env (score.py imports pysam through pool.py)."""
import os
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import ltrkin_breakdown as K  # noqa: E402
import nonauto_sim  # noqa: E402


def call(tool, s, e, chrom="c1"):
    return dict(chrom=chrom, start=s, end=e, tool=tool)


def locus(i, s, e, chrom="c1"):
    return dict(locus_id=f"L{i}", chrom=chrom, start=s, end=e)


def test_consensus_needs_two_callers_and_takes_medians():
    calls = [call("a", 0, 10000), call("b", 0, 10010), call("a", 30000, 31000),       # 2-caller cluster; 1-caller
             call("a", 50000, 60000), call("b", 50000, 60000), call("c", 50002, 60004)]
    cons = K.consensus_elements(calls)
    assert [(e["start"], e["end"], e["callers"]) for e in cons] == [(0, 10005, "a,b"), (50000, 60000, "a,b,c")]


def test_nesting_rule():
    calls = [call("a", 0, 10000), call("b", 0, 10010),               # A: consensus [0, 10005)
             call("a", 20000, 25000), call("b", 20000, 25000),       # B: [20000, 25000), matches L2
             call("a", 30000, 31000),                                # one caller only: no consensus element
             call("a", 42000, 44000), call("b", 42000, 44000),       # D: inside L4, not matching it
             call("a", 49800, 55300), call("b", 49800, 55300),       # E: encloses L5, only 500 bp longer
             call("a", 60000, 80000), call("b", 60000, 80000),       # F: encloses L6 (>= 1 kb longer)
             call("a", 66000, 68000), call("b", 66000, 68000)]       # G: inside L6
    cons = K.consensus_elements(calls)
    loci = [locus(0, 0, 10000),        # matches A, contains pooled locus L1 -> host
            locus(1, 3000, 6000),      # inside A (10005 >= 3000 + 1000) -> tenant
            locus(2, 20000, 25000),    # only its own consensus element -> unnested
            locus(3, 29000, 36000),    # contains a one-caller cluster only -> unnested
            locus(4, 40000, 48000),    # contains D -> host
            locus(5, 50000, 55000),    # inside E, but E only 500 bp longer -> unnested
            locus(6, 64000, 70000)]    # inside F and contains G -> both -> host
    n = K.nesting(loci, cons)
    assert {k: v["nesting"] for k, v in n.items()} == {"L0": "host", "L1": "tenant", "L2": "unnested",
                                                        "L3": "unnested", "L4": "host", "L5": "unnested", "L6": "host"}
    assert n["L0"]["contained_loci"] == "L1" and n["L0"]["n_contained"] == 0 and n["L0"]["n_self_cons"] == 1
    assert n["L1"]["enclosing"] == "c1:0-10005" and n["L1"]["enclosing_len"] == "10005"
    assert n["L4"]["contained"] == "c1:42000-44000" and n["L4"]["contained_ro"] == "0.250"
    assert n["L6"]["is_tenant"] == 1 and n["L6"]["is_host"] == 1


def test_gff_1based_to_0based_and_renaming(tmp_path):
    gff = tmp_path / "dante.gff3"
    attrs = "Name={n};Final_Classification=Class_I|LTR|Ty1/copia|Ale"
    gff.write_text("##gff-version 3\n"
                   + f"c0001\tdante\tprotein_domain\t11\t20\t50\t+\t.\t{attrs.format(n='RT')}\n"
                   + f"c0001\tdante\tprotein_domain\t31\t40\t50\t+\t.\t{attrs.format(n='INT')}\n")
    names = tmp_path / "names.tsv"
    names.write_text("new_id\toriginal_id\tlength\tn_count\tnonacgtn_converted\nc0001\tCM1\t1000\t0\t0\n")
    fai = tmp_path / "ref.fa.fai"
    fai.write_text("CM1\t1000\t5\t60\t61\n")
    doms = K.load_domains(str(gff), str(names), str(fai))
    assert list(doms) == ["CM1"]
    assert [r[:2] for r in doms["CM1"][1]] == [(10, 20), (30, 40)]          # 1-based 11..20 -> [10, 20)
    # off-by-one sensitive: [10, 20) is inside a locus starting at 10, not one starting at 11
    assert K.domain_evidence(doms, "CM1", 10, 20)[0] == 1
    assert K.domain_evidence(doms, "CM1", 11, 30)[0] == 0
    # protein class (nonauto rule): 2 domain names of one strand/superfamily in the middle 80% -> coding
    assert nonauto_sim.call_class(doms, dict(chrom="CM1", start=5, end=45)) == "coding"
    assert nonauto_sim.call_class(doms, dict(chrom="CM1", start=5, end=25)) == "partial"
    assert nonauto_sim.call_class(doms, dict(chrom="CM1", start=100, end=500)) == "proteinless"


def test_fasta_length_mismatch_fails(tmp_path):
    gff = tmp_path / "dante.gff3"
    gff.write_text("##gff-version 3\n")
    names = tmp_path / "names.tsv"
    names.write_text("new_id\toriginal_id\tlength\tn_count\tnonacgtn_converted\nc0001\tCM1\t1000\t0\t0\n")
    fai = tmp_path / "ref.fa.fai"
    fai.write_text("CM1\t999\t5\t60\t61\n")
    try:
        K.load_domains(str(gff), str(names), str(fai))
    except SystemExit as e:
        assert "length" in str(e)
    else:
        raise AssertionError("length mismatch not caught")


def test_cat_of_matches_fig1_mapping():
    assert K.cat_of(("missed", None), ("missed", False)) == "missed"
    assert K.cat_of(("partial", None), ("partial", False)) == "partial"
    assert K.cat_of(("found", 3), ("F/E", True)) == "correct"
    for g in ("F/F", "F/N", "F/S", "no_call", "other", "not_genotypable"):
        assert K.cat_of(("found", 3), (g, False)) == "wrong"
