"""phase_inputs.py --max-len: long spans skip Kmer2LTR and are listed as too_long (fake Kmer2LTR on PATH)."""
import csv
import os
import stat
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts"))
import phase_inputs  # noqa: E402

FAKE = """#!/usr/bin/env python3
import sys
fa, out = sys.argv[1], sys.argv[sys.argv.index("-o") + 1]
names = [l[1:].split()[0] for l in open(fa) if l.startswith(">")]
open(fa + ".seen", "w").write("\\n".join(names))
with open(out, "w") as f:
    f.write("seq_id\\tseq_len\\tstatus\\tltr5_start\\tltr5_end\\tltr3_start\\tltr3_end\\tk2p\\n")
    for n in names:
        f.write(f"{n}\\t500\\tpass\\t1\\t100\\t401\\t500\\t0.01\\n")
"""


def test_max_len(tmp_path, monkeypatch):
    seq = "ACGT" * 1250                                   # 5,000 bp
    (tmp_path / "g.fa").write_text(f">c1\n{seq}\n")
    (tmp_path / "g.fa.fai").write_text(f"c1\t5000\t4\t5000\t5001\n")
    (tmp_path / "calls.tsv").write_text("chrom\tstart\tend\nc1\t100\t600\nc1\t1000\t4000\n")
    bin_ = tmp_path / "bin"
    bin_.mkdir()
    k = bin_ / "Kmer2LTR"
    k.write_text(FAKE)
    k.chmod(k.stat().st_mode | stat.S_IEXEC)
    monkeypatch.setenv("PATH", f"{bin_}:{os.environ['PATH']}")
    out = tmp_path / "o"
    phase_inputs.main(["--calls", str(tmp_path / "calls.tsv"), "--genome", str(tmp_path / "g.fa"),
                       "--outdir", str(out), "-t", "1", "--max-len", "1000"])
    rows = {r["seq_id"]: r for r in csv.DictReader(open(out / "annot.tsv"), delimiter="\t")}
    assert rows["c1:101-600"]["status"] == "pass"
    assert rows["c1:1001-4000"]["status"] == "too_long" and rows["c1:1001-4000"]["k2p"] == "NA"
    assert (out / "elements.fa.seen").read_text().split() == ["c1:101-600"]      # Kmer2LTR saw only the short one
    assert [l.strip() for l in open(out / "elements.k.fa") if l.startswith(">")] == [">c1:101-600"]
