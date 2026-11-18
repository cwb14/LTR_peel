"""Tiny synthetic fixture for sim_class_prf_joint (0-based half-open coordinates throughout)."""
import math
import os
import sys

import pytest

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sim_class_prf_joint as J  # noqa: E402

LIB = {"famC": {"cls": "coding"}, "famP": {"cls": "proteinless"}, "famQ": {"cls": "partial"}}


def t(cls, s, e, id_, fam):
    return dict(id=id_, chrom="chr1", start=s, end=e, depth=0, tenant_in_ltr=0, std=int(cls != "host"), family=fam,
                **{"class": cls})


def c(s, e, id_):
    return dict(chrom="chr1", start=s, end=e, id=id_, parent=".")


def doms_of(recs, chrom="chr1"):
    """evidence.dante_domains layout: {chrom: (starts, records sorted)}, record = (start0, end, strand, name, sf)."""
    recs = sorted(recs)
    return {chrom: ([r[0] for r in recs], recs)}


def fixture():
    # host H 1000-11000 (coding family) split by tenant T 5000-7000 (protein-less family); un-nested U 20000-25000
    # (coding family); un-nested M 30000-34000 (partial family) that no call finds
    truth = [t("host", 1000, 11000, "H", "famC"), t("tenant", 5000, 7000, "T", "famP"),
             t("intact", 20000, 25000, "U", "famC"), t("intact", 30000, 34000, "M", "famQ")]
    feats = [dict(chrom="chr1", start=1000, end=5000, cat="LTR_piece"),
             dict(chrom="chr1", start=7000, end=11000, cat="LTR_piece"),
             dict(chrom="chr1", start=5000, end=7000, cat="LTR_intact"),
             dict(chrom="chr1", start=20000, end=25000, cat="LTR_intact"),
             dict(chrom="chr1", start=30000, end=34000, cat="LTR_intact")]
    # genome DANTE domains: RT + INT in the host's first piece, GAG + RT in U; none in the tenant or M
    doms = doms_of([(2500, 3000, "+", "RT", "Ty1/copia"), (3200, 3600, "+", "INT", "Ty1/copia"),
                    (21000, 21500, "+", "GAG", "Ty1/copia"), (22000, 22500, "+", "RT", "Ty1/copia")])
    calls = [c(1000, 11000, "c_host"),       # matches H; RT + INT in its middle 80% -> coding
             c(5000, 7000, "c_tenant"),      # matches T; no domain -> proteinless
             c(20000, 25000, "c_un"),        # matches U -> coding
             c(1000, 6000, "c_partial"),     # unmatched, partial_host -> host; RT + INT inside -> coding
             c(2400, 3100, "c_frag"),        # unmatched, partial_host -> host; RT only -> partial
             c(40000, 41000, "c_none")]      # unmatched, no feature -> unnested; no domain -> proteinless
    return truth, feats, doms, calls


def test_both_systems_on_one_matching():
    truth, feats, doms, calls = fixture()
    ecls = J.element_classes(truth, LIB)
    assert ecls == ["coding", "proteinless", "coding", "partial"]
    k, kk, removed, c2t, t2c = J.joint_counts(calls, truth, feats, ecls, doms)
    assert c2t == {0: 0, 1: 1, 2: 2} and kk == k and removed == 0
    n = k["nesting"]
    assert n["host"] == dict(n_truth=1, tp_truth=1, n_calls=3, tp_calls=1)
    assert n["tenant"] == dict(n_truth=1, tp_truth=1, n_calls=1, tp_calls=1)
    assert n["unnested"] == dict(n_truth=2, tp_truth=1, n_calls=2, tp_calls=1)
    p = k["protein"]
    assert p["coding"] == dict(n_truth=2, tp_truth=2, n_calls=3, tp_calls=2)
    assert p["partial"] == dict(n_truth=1, tp_truth=0, n_calls=1, tp_calls=0)
    assert p["proteinless"] == dict(n_truth=1, tp_truth=1, n_calls=2, tp_calls=1)
    # both systems sum to the same totals (the ext row)
    for kk in k.values():
        assert {x: sum(v[x] for v in kk.values()) for x in J.KEYS} == dict(n_truth=4, tp_truth=3, n_calls=6, tp_calls=3)
    pr, rc, f = J.SC.prf(p["coding"])
    assert (pr, rc) == (2 / 3, 1.0) and abs(f - 0.8) < 1e-12
    pr, rc, f = J.SC.prf(p["proteinless"])
    assert (pr, rc) == (0.5, 1.0) and abs(f - 2 / 3) < 1e-12
    pr, rc, f = J.SC.prf(p["partial"])
    assert (pr, rc) == (0.0, 0.0) and math.isnan(f)          # P = R = 0: F1 NA (sim_class_prf convention)


def sums(kk):
    return {s: {x: sum(v[x] for v in kc.values()) for x in J.KEYS} for s, kc in kk.items()}


def test_hasltr_removes_element_and_its_call():
    truth, feats, doms, calls = fixture()
    ecls = J.element_classes(truth, LIB)
    # remove U (matched by c_un): it leaves the truth, c_un leaves the calls
    k, kh, removed, _, _ = J.joint_counts(calls, truth, feats, ecls, doms, [True, True, False, True])
    assert removed == 1
    assert kh["nesting"]["unnested"] == dict(n_truth=1, tp_truth=0, n_calls=1, tp_calls=0)
    assert kh["nesting"]["host"] == k["nesting"]["host"] and kh["nesting"]["tenant"] == k["nesting"]["tenant"]
    assert kh["protein"]["coding"] == dict(n_truth=1, tp_truth=1, n_calls=2, tp_calls=1)
    assert kh["protein"]["proteinless"] == k["protein"]["proteinless"]
    s = sums(kh)
    assert s["nesting"] == s["protein"] == dict(n_truth=3, tp_truth=2, n_calls=5, tp_calls=2)


def test_hasltr_unmatched_calls_on_a_removed_element_stay():
    truth, feats, doms, calls = fixture()
    ecls = J.element_classes(truth, LIB)
    # remove host H: c_host leaves; the unmatched partial-host calls c_partial, c_frag stay as host false positives
    _, kh, removed, _, _ = J.joint_counts(calls, truth, feats, ecls, doms, [False, True, True, True])
    assert removed == 1
    assert kh["nesting"]["host"] == dict(n_truth=0, tp_truth=0, n_calls=2, tp_calls=0)
    p, r, _ = J.SC.prf(kh["nesting"]["host"])
    assert p == 0.0 and math.isnan(r)
    assert kh["protein"]["coding"] == dict(n_truth=1, tp_truth=1, n_calls=2, tp_calls=1)
    assert kh["protein"]["partial"] == dict(n_truth=1, tp_truth=0, n_calls=1, tp_calls=0)
    # removing a missed element (M) removes no call
    _, kh, removed, _, _ = J.joint_counts(calls, truth, feats, ecls, doms, [True, True, True, False])
    assert removed == 0
    s = sums(kh)
    assert s["nesting"] == s["protein"] == dict(n_truth=3, tp_truth=3, n_calls=6, tp_calls=3)


def test_crosstab():
    truth, _, _, _ = fixture()
    ecls = J.element_classes(truth, LIB)
    assert J.crosstab(truth, ecls) == {("host", "coding"): 1, ("tenant", "proteinless"): 1,
                                       ("unnested", "coding"): 1, ("unnested", "partial"): 1}
    assert J.crosstab(truth, ecls, [False, True, False, True]) == {("tenant", "proteinless"): 1,
                                                                   ("unnested", "partial"): 1}


def test_unknown_family_fails():
    truth, _, _, _ = fixture()
    truth.append(t("intact", 50000, 52000, "X", "famX"))
    with pytest.raises(SystemExit):
        J.element_classes(truth, LIB)


def test_summary_skips_undefined_genomes():
    rows = [dict(set="S2", seed=s, tool="x", system="protein", cls="partial", n_truth=10, tp_truth=5, n_calls=nc,
                 tp_calls=tc, precision=p, recall=0.5, f1=f, role="primary")
            for s, nc, tc, p, f in (("1", 2, 1, 0.5, 0.5), ("2", 0, 0, math.nan, math.nan), ("3", 10, 7, 0.7, 0.5833))]
    (d,) = J.summarize(rows)
    assert d["precision_n_genomes"] == 2 and abs(d["precision_mean"] - 0.6) < 1e-12
    assert abs(d["precision_sd"] - 0.1414213562) < 1e-9
    assert d["recall_n_genomes"] == 3 and d["recall_sd"] == 0.0
    assert (d["n_truth"], d["n_calls"]) == (30, 12)


def test_terminal_repeat_and_library_rule():
    L = 5000
    assert J.terminal_repeat([(1, 400, 4601, 5000)], L, 50) == 400
    assert J.terminal_repeat([(10, 400, 4611, 4990)], L, 50) == 400           # arms within 50 bp of the ends
    assert J.terminal_repeat([(100, 500, 4000, 4400)], L, 50) is None         # not terminal
    assert J.terminal_repeat([(1, 3000, 2001, 5000)], L, 50) is None          # overlapping arms
    # GAG inside the 5' repeat, chromodomain in the internal region: coding on the whole entry, partial between repeats
    doms = doms_of([(100, 500, "+", "GAG", "Ty3/gypsy"), (4300, 4450, "+", "CHD", "Ty3/gypsy")], "L0001")
    assert J.lib_class_of(doms, "L0001", 0, L, 2) == "coding"
    assert J.lib_class_of(doms, "L0001", 400, L - 400, 2) == "partial"
    assert J.lib_class_of({}, "L0002", 0, L, 0) == "proteinless"
