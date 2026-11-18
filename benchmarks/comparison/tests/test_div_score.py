import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
import div_score as ds  # noqa: E402


def test_k2p_identical_and_known():
    assert ds.k2p("ACGTACGTAC", "ACGTACGTAC") == 0
    assert f"{ds.k2p('ACGT', 'ACGT'):.4f}" == "0.0000"   # not -0.0000
    # 1 transition (A->G) and 1 transversion (C->A) in 10 sites: P = Q = 0.1
    d = ds.k2p("ACGTACGTAC", "GAGTACGTAC")
    want = -0.5 * math.log(1 - 0.2 - 0.1) - 0.25 * math.log(1 - 0.2)
    assert abs(d - want) < 1e-12
    assert math.isnan(ds.k2p("NNNN", "ACGT"))  # no comparable sites


def test_snap_names():
    assert ds.snap(0) == "burnin"
    assert ds.snap(3) == "gen4620000_final"


def test_d50_interpolates_between_timepoints():
    pts = [(0, 0.0, 0.9), (1, 0.01, 0.8), (2, 0.05, 0.3)]
    # half of 0.9 = 0.45: between x 0.01 (0.8) and 0.05 (0.3) -> 0.01 + 0.04 * (0.8 - 0.45) / 0.5 = 0.038
    assert ds.d50(pts) == "0.0380"


def test_d50_not_reached_and_incomplete():
    full = [(k, 0.04 * k, 0.9) for k in range(11)]
    assert ds.d50(full) == ">0.4000"
    assert ds.d50(full[:5]) == "incomplete"                      # later timepoints not scored yet
    assert ds.d50([(0, 0.0, 0.9), (2, 0.05, 0.1)]) == "incomplete"  # gap
    assert ds.d50([(1, 0.01, 0.5)]) == "NA"                      # no burn-in point


def test_wilson_bounds():
    lo, hi = ds.wilson(5, 10)
    assert 0 < lo < 0.5 < hi < 1
    assert all(math.isnan(x) for x in ds.wilson(0, 0))


def _run(tmp_path, ds, tool, files):
    d = tmp_path / "runs" / ds / tool
    (d / "out").mkdir(parents=True)
    (d / ".failed").write_text("rc=1\n")
    for name, text in files.items():
        (d / name).write_text(text)


def test_empty_candidate_exits(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    # LTR_retriever: empty rawLTR.scn (comments only) + its "No candidate" error -> no_candidates
    _run(tmp_path, "SD1_09", "ltr_retriever", {"SD1_09.fa.rawLTR.scn": "#comment\n\n",
                                               "run.log": "ERROR: No candidate is found in the file(s)\n"})
    assert ds.status("SD1_09", "ltr_retriever") == "no_candidates"
    # same error text but a candidate line present -> a real failure
    _run(tmp_path, "SD1_08", "ltr_retriever", {"SD1_08.fa.rawLTR.scn": "#c\n1 2 3\n",
                                               "run.log": "ERROR: No candidate is found\n"})
    assert ds.status("SD1_08", "ltr_retriever") == "failed"
    # Look4LTRs: empty element database + abort in identity training -> no_candidates; non-empty db -> failed
    _run(tmp_path, "SD1_10", "look4ltrs", {"out/db.fasta": "", "run.log": "Training Identity...\nAborted\n"})
    assert ds.status("SD1_10", "look4ltrs") == "no_candidates"
    _run(tmp_path, "SD2_10", "look4ltrs", {"out/db.fasta": ">e1\nACGT\n", "run.log": "Training Identity...\n"})
    assert ds.status("SD2_10", "look4ltrs") == "failed"
    # other tools never get the empty-candidate reading
    _run(tmp_path, "SD1_10", "hite", {"run.log": "No candidate is found\n"})
    assert ds.status("SD1_10", "hite") == "failed"
    # calls present -> ok
    (tmp_path / "calls" / "SD1_10").mkdir(parents=True)
    (tmp_path / "calls" / "SD1_10" / "ltrharvest.tsv").write_text("chrom\tstart\tend\n")
    assert ds.status("SD1_10", "ltrharvest") == "ok"


def test_element_stats_ltr_k2p_and_tsd(tmp_path):
    # element [5, 25) with 4-bp LTRs and 3-bp TSDs at [2, 5) and [25, 28)
    seq = "NN" + "AAG" + "ACGT" + "C" * 12 + "ACGA" + "AAG" + "TTTT"
    fa = tmp_path / "g.fa"
    fa.write_text(">c1\n" + seq + "\n")
    out = ds.element_stats(str(fa), [("c1", 5, 25, 4, 3), ("c1", 5, 25, 4, 0), ("c1", 1, 25, 4, 3)])
    d, same, p = out[("c1", 5, 25)]   # last tuple for this key wins: tsd_len 0 -> None
    assert same is None
    out = ds.element_stats(str(fa), [("c1", 5, 25, 4, 3)])
    d, same, p = out[("c1", 5, 25)]
    assert same is True and d > 0 and p == 0.25   # one transversion (T->A) in 4 LTR sites
    seq2 = seq[:25] + "AAT" + seq[28:]
    fa.write_text(">c1\n" + seq2 + "\n")
    assert ds.element_stats(str(fa), [("c1", 5, 25, 4, 3)])[("c1", 5, 25)][1] is False
    assert ds.element_stats(str(fa), [("c1", 1, 25, 4, 3)])[("c1", 1, 25)][1] is None   # left TSD off the sequence


def test_tsd_adjusted_reweights_strata():
    # identical-TSD recall 0.8, differ 0.1; observed strata 60/40, expected identical share 0.9 -> 0.8*0.9 + 0.1*0.1
    assert abs(ds.tsd_adjusted([60, 48, 40, 4], 0.9) - 0.73) < 1e-12
    assert ds.tsd_adjusted([100, 90, 0, 0], 1.0) == 0.9            # burn-in: nothing to reweight
    assert math.isnan(ds.tsd_adjusted([100, 90, 0, 0], 0.8))       # no differ stratum to reweight from
    assert ds.fmt(float("nan")) == "NA" and ds.fmt(0.5) == "0.5000"
