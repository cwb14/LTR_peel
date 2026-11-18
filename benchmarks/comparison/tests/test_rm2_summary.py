"""BM_RM2 family counts (scripts/rm2_summary.py) from the RepeatModeler2 script's family lists."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts"))
import rm2_summary  # noqa: E402


def test_counts(tmp_path):
    (tmp_path / "cur.fa").write_text(">A_LTR#LTR/Copia\nAC\n>A_I#LTR/Copia\nAC\n>B#DNA/hAT\nAC\n>C#LTR/Gypsy\nAC\n")
    w = tmp_path / "w"
    w.mkdir()
    (w / "perfect.families").write_text("A_LTR\n")
    (w / "good.families").write_text("A_I\nB\n")
    (w / "present.all.families").write_text("A_LTR\nA_I\nB\nC\n")
    out = tmp_path / "o.tsv"
    rm2_summary.main(["--workdir", str(w), "--curated", str(tmp_path / "cur.fa"), "--out", str(out)])
    rows = [l.split("\t") for l in out.read_text().splitlines()[1:]]
    assert rows[0] == ["LTR", "3", "1", "1", "2", "3"]
    assert rows[1] == ["all", "4", "1", "2", "3", "4"]
