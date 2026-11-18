"""Every tool parser on fixture records cut from real smoke-test outputs (A. thaliana chr2, normalized ID c0001).

The same element (1-based 16330..24007) is reported identically by LTRquest, LTR_Finder, LTR_retriever and DANTE_LTR,
so its 0-based span must be 16329..24007 from each; solo LTRs (Look4LTRs) and domain-only DANTE elements are dropped.
"""
import csv
import io
import os
import subprocess
import sys

H = os.path.dirname(__file__)
F = os.path.join(H, "fixtures", "calls")
S = os.path.join(H, "..", "scripts", "calls_to_bed.py")


def parse(tool, keep=False, rundir=None):
    args = [sys.executable, S, tool, os.path.join(F, rundir or tool), os.path.join(F, "data", "fx"), "/dev/stdout"]
    r = subprocess.run(args + (["--keep-ids"] if keep else []), capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return list(csv.DictReader(io.StringIO(r.stdout), delimiter="\t"))


def spans(rows):
    return [(r["chrom"], int(r["start"]), int(r["end"])) for r in rows]


def test_shared_element_identical_across_tools():
    for tool in ("ltrquest", "ltr_finder", "ltr_retriever", "dante_ltr"):
        assert ("Chr2", 16329, 24007) in spans(parse(tool)), tool


def test_ltrquest_nest_link_and_ltrs():
    rows = {r["id"]: r for r in parse("ltrquest")}
    t = rows["c0001:16330-24007"]
    assert t["parent"] == "c0001:13072-28021" and t["parent"] in rows
    assert (int(t["lltr_end"]), int(t["rltr_start"])) == (17647, 22689)


def test_ltrquest_hc_subset_drops_link_to_excluded_host():
    rows = parse("ltrquest_hc", rundir="ltrquest")      # tenant is high confidence, its host medium
    assert spans(rows) == [("Chr2", 16329, 24007)]
    assert rows[0]["parent"] == "."


def test_harvest_and_ltrpred_one_based():
    assert spans(parse("ltrharvest"))[0] == ("Chr2", 2079, 10521)
    assert spans(parse("ltrpred"))[0] == ("Chr2", 2079, 10521)       # LTRpred runs LTRharvest internally


def test_ltrpred_too_few_loci_is_zero_calls_but_missing_gff_fails():
    assert parse("ltrpred", rundir="ltrpred_fewloci") == []           # LTRpred's designed no-prediction outcome
    args = [sys.executable, S, "ltrpred", os.path.join(F, "ltrpred_nogff"), os.path.join(F, "data", "fx"), "/dev/null"]
    r = subprocess.run(args, capture_output=True, text=True)
    assert r.returncode != 0 and "no file matches" in r.stderr          # missing LTRdigest output stays a failure


def test_look4ltrs_zero_based_and_solo_dropped():
    rows = parse("look4ltrs")
    assert spans(rows) == [("Chr2", 16321, 24009), ("Chr2", 17630, 45372)]
    assert rows[1]["strand"] == "-"


def test_hirepeat_index_and_hite_rename():
    assert spans(parse("hirepeat"))[0] == ("Chr2", 13071, 28030)
    assert spans(parse("hite"))[0] == ("Chr2", 2902329, 2907117)


def test_dante_domain_only_dropped():
    rows = parse("dante_ltr")
    assert len(rows) == 1 and rows[0]["id"].startswith("TE_0")


def test_keep_ids():
    assert parse("ltr_finder", keep=True)[0]["chrom"] == "c0001"
