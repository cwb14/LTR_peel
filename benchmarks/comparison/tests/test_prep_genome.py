import gzip
import os
import subprocess
import sys

SCRIPT = os.path.join(os.path.dirname(__file__), "..", "scripts", "prep_genome.py")


def run(args):
    return subprocess.run([sys.executable, SCRIPT] + args, capture_output=True, text=True)


def test_normalize_and_index(tmp_path):
    src = tmp_path / "in.fa.gz"
    seq1 = "acgtRYacgt" * 13          # 130 bp, lower case + IUPAC
    seq2 = "NNNNACGT"
    with gzip.open(src, "wt") as f:
        f.write(f">chr_1 desc\n{seq1[:70]}\n{seq1[70:]}\n>scaf_2\n{seq2}\n")
    out = tmp_path / "out"
    r = run([str(src), str(out), "--split"])
    assert r.returncode == 0, r.stderr
    fa = (out / "genome.fa").read_text()
    assert fa.startswith(">c0001\n")
    body1 = "".join(fa.split(">c0002")[0].splitlines()[1:])
    assert body1 == ("ACGTNNACGT" * 13)
    assert all(len(l) <= 60 for l in fa.splitlines())
    # .fai offsets must let us seek to the sequence start
    fai = [l.split("\t") for l in (out / "genome.fa.fai").read_text().splitlines()]
    assert fai[0][:2] == ["c0001", "130"] and fai[1][:2] == ["c0002", "8"]
    with open(out / "genome.fa", "rb") as f:
        f.seek(int(fai[1][2]))
        assert f.read(8) == b"NNNNACGT"
    names = (out / "names.tsv").read_text().splitlines()
    assert names[1].split("\t") == ["c0001", "chr_1", "130", "26", "26"]
    assert (out / "split" / "c0002.fa").read_text() == ">c0002\nNNNNACGT\n"


def test_duplicate_ids_fail(tmp_path):
    src = tmp_path / "in.fa"
    src.write_text(">a\nACGT\n>a\nACGT\n")
    r = run([str(src), str(tmp_path / "o")])
    assert r.returncode != 0 and "duplicate" in r.stderr


def test_bad_prefix(tmp_path):
    src = tmp_path / "in.fa"
    src.write_text(">a\nACGT\n")
    r = run([str(src), str(tmp_path / "o"), "--prefix", "c_"])
    assert r.returncode != 0
