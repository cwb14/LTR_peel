"""printe_truth.py + score_elements.py on a hand-made PrinTE BED with one nested host."""
import csv
import os
import subprocess
import sys

S = os.path.join(os.path.dirname(__file__), "..", "scripts")

# host H (LTRlen 100) of 1000 bp cut at offset 400 by tenant T (500 bp, TSD 'ACGTA', 5 bp gap after it);
# T is itself intact; plus a free intact element F, a solo, a gene, and a DNA TE.
BED = [
    ("chr1", 1000, 1400, "H#LTR/Gypsy~LTRlen:100;CUT_BY:T#LTR/Copia~LTRlen:50", "AAAAA", "+"),
    ("chr1", 1400, 1900, "T#LTR/Copia~LTRlen:50;NESTED_IN:H#LTR/Gypsy~LTRlen:100", "ACGTA", "-"),
    ("chr1", 1905, 2505, "H#LTR/Gypsy~LTRlen:100;CUT_BY:T#LTR/Copia~LTRlen:50", "AAAAA", "+"),
    ("chr1", 5000, 7000, "F#LTR/Copia~LTRlen:200", "CCCCC", "+"),
    ("chr1", 9000, 9300, "S#LTR/Gypsy~LTRlen:300_SOLO", "GGGGG", "+"),
    ("chr1", 12000, 13000, "gene1", "NA", "+"),
    ("chr1", 15000, 15500, "D#DNA/DTA", "TTT", "+"),
    # a second copy of H far away with the same TSD must NOT be merged with the first
    ("chr1", 40000, 41000, "H#LTR/Gypsy~LTRlen:100", "AAAAA", "+"),
]


def run(script, args):
    r = subprocess.run([sys.executable, os.path.join(S, script)] + args, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return r


def truth(tmp_path):
    bed = tmp_path / "t.bed"
    bed.write_text("".join("\t".join(map(str, x)) + "\n" for x in BED))
    run("printe_truth.py", [str(bed), str(tmp_path / "T")])
    return {r["id"]: r for r in csv.DictReader(open(tmp_path / "T.truth.tsv"), delimiter="\t")}


def test_truth_classes(tmp_path):
    t = truth(tmp_path)
    assert set(t) == {"chr1:1000-2505", "chr1:1400-1900", "chr1:5000-7000", "chr1:40000-41000"}
    h = t["chr1:1000-2505"]
    assert h["class"] == "host" and h["n_pieces"] == "2" and h["seq_len"] == "1000" and h["depth"] == "1"
    assert h["tenant_in_ltr"] == "0" and h["std"] == "0"          # offset 400: internal region
    assert t["chr1:1400-1900"]["class"] == "tenant" and t["chr1:1400-1900"]["parent"] == "chr1:1000-2505"
    assert t["chr1:1400-1900"]["std"] == "1" and t["chr1:5000-7000"]["std"] == "1"
    pairs = list(csv.DictReader(open(tmp_path / "T.pairs.tsv"), delimiter="\t"))
    assert [(p["host_id"], p["tenant_id"]) for p in pairs] == [("chr1:1000-2505", "chr1:1400-1900")]
    std = (tmp_path / "T.std.bed").read_text().split()
    assert len(std) == 9                                           # 3 entries x 3 columns (T, F, far H copy)


def test_score(tmp_path):
    truth(tmp_path)
    calls = tmp_path / "c.tsv"
    calls.write_text("chrom\tstart\tend\tid\tparent\n"
                     "chr1\t1000\t2505\th\t.\n"          # host, exact
                     "chr1\t1401\t1899\tt\th\n"          # tenant, 1 bp off each end, linked to host
                     "chr1\t9000\t9300\ts\t.\n"          # solo LTR -> FP
                     "chr1\t5000\t5600\tp\t.\n")         # partial of F -> FP
    run("score_elements.py", ["--calls", str(calls), "--truth", str(tmp_path / "T"), "--out", str(tmp_path / "s")])
    m = {(r["scope"], r["ro"]): r for r in csv.DictReader(open(tmp_path / "s.metrics.tsv"), delimiter="\t")}
    std, ext = m[("std", "0.9")], m[("ext", "0.9")]
    assert std["tp"] == "1" and std["n_truth"] == "3"                 # only the tenant matches the standard set
    assert ext["tp"] == "2" and ext["n_truth"] == "4" and ext["precision"] == "0.5000"
    assert m[("ext:host", "0.9")]["recall"] == "1.0000"
    assert m[("nest_links", "0.9")]["tp"] == "1" and m[("nest_links", "0.9")]["recall"] == "1.0000"
    an = {r["category"]: int(r["n"]) for r in csv.DictReader(open(tmp_path / "s.anatomy.tsv"), delimiter="\t")}
    assert an == {"LTR_solo": 1, "partial_intact": 1}
    b = {r["end"]: r for r in csv.DictReader(open(tmp_path / "s.boundaries.tsv"), delimiter="\t")}
    assert b["left"]["exact"] == "0.5000" and b["left"]["le5"] == "1.0000"
