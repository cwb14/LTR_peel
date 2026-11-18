"""Tiny synthetic fixture for sim_class_prf.class_counts (0-based half-open coordinates throughout)."""
import os
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sim_class_prf as S  # noqa: E402


def t(cls, s, e, id_):
    return dict(id=id_, chrom="chr1", start=s, end=e, depth=0, tenant_in_ltr=0, std=int(cls != "host"), **{"class": cls})


def test_one_host_one_tenant_one_unnested():
    # host 1000-11000 split by tenant 5000-7000; un-nested element 20000-25000
    truth = [t("host", 1000, 11000, "H"), t("tenant", 5000, 7000, "T"), t("intact", 20000, 25000, "U")]
    feats = [dict(chrom="chr1", start=1000, end=5000, cat="LTR_piece"),
             dict(chrom="chr1", start=7000, end=11000, cat="LTR_piece"),
             dict(chrom="chr1", start=5000, end=7000, cat="LTR_intact"),
             dict(chrom="chr1", start=20000, end=25000, cat="LTR_intact")]
    calls = [dict(chrom="chr1", start=1000, end=11000, id="c_host", parent="."),       # matches the host
             dict(chrom="chr1", start=5000, end=7000, id="c_tenant", parent="."),      # matches the tenant
             dict(chrom="chr1", start=20000, end=25000, id="c_un", parent="."),        # matches the un-nested element
             dict(chrom="chr1", start=1000, end=6000, id="c_partial", parent="."),     # partial call on the host
             dict(chrom="chr1", start=40000, end=41000, id="c_none", parent=".")]      # on no element
    k, cats, c2t, t2c, unmatched = S.class_counts(calls, truth, feats)
    assert cats == {3: "partial_host", 4: "no_feature"}
    assert k["host"] == dict(n_truth=1, tp_truth=1, n_calls=2, tp_calls=1)
    assert k["tenant"] == dict(n_truth=1, tp_truth=1, n_calls=1, tp_calls=1)
    assert k["unnested"] == dict(n_truth=1, tp_truth=1, n_calls=2, tp_calls=1)
    p, r, f = S.prf(k["host"])
    assert (p, r) == (0.5, 1.0) and abs(f - 2 / 3) < 1e-12


def test_host_tenant_counts_as_host_and_fp_mapping():
    assert S.TRUTH_CLASS["host_tenant"] == "host"
    assert S.fp_class("partial_host_tenant") == "host"
    assert S.fp_class("partial_tenant") == "tenant"
    assert S.fp_class("partial_intact") == "unnested"
    for c in ("LTR_frag", "LTR_solo", "gene", "nonLTR_TE", "nonLTR_frag", "no_feature", "anything_else"):
        assert S.fp_class(c) == "unnested"
