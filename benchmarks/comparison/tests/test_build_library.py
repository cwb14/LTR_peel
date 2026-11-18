"""build_library.py: intact mode (cd-hit-est removes only near-identical copies) and native header normalization."""
import os
import random
import shutil
import subprocess
import sys

import pytest

H = os.path.dirname(__file__)
S = os.path.join(H, "..", "scripts", "build_library.py")
sys.path.insert(0, os.path.join(H, "..", "scripts"))
import build_library as bl  # noqa: E402


def read_fa(path):
    return dict(bl.read_fa(path))


def test_cdhit_word_sizes_follow_user_guide():
    assert [bl.cdhit_word(c) for c in (1.0, 0.95, 0.9, 0.88, 0.85, 0.8)] == [10, 10, 8, 7, 6, 5]


@pytest.mark.skipif(shutil.which("cd-hit-est") is None, reason="cd-hit-est not on PATH")
def test_intact_keeps_one_of_identical_copies(tmp_path):
    rnd = random.Random(7)
    ele_a = "".join(rnd.choice("ACGT") for _ in range(1500))
    ele_b = "".join(rnd.choice("ACGT") for _ in range(1500))
    pad = lambda n: "".join(rnd.choice("ACGT") for _ in range(n))  # noqa: E731
    g = pad(500) + ele_a + pad(500) + ele_a + pad(500) + ele_b + pad(500)
    gf = tmp_path / "g.fa"
    gf.write_text(">c0001\n" + g + "\n")
    (tmp_path / "g.fa.fai").write_text(f"c0001\t{len(g)}\t7\t{len(g)}\t{len(g) + 1}\n")
    s1, s2, s3 = 500, 2500, 4500                       # 0-based starts of the three copies
    calls = tmp_path / "calls.tsv"
    calls.write_text("chrom\tstart\tend\tclass\n" + "".join(
        f"c0001\t{s}\t{s + 1500}\t{c}\n" for s, c in ((s1, "Gypsy"), (s2, "Gypsy"), (s3, "Copia"))))
    out = tmp_path / "lib"
    r = subprocess.run([sys.executable, S, "--calls", str(calls), "--genome", str(gf), "--out", str(out),
                        "--prefix", "tool", "-t", "1", "--intact", "0.95"], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    lib = read_fa(str(out) + ".fa")
    assert len(lib) == 2
    assert sorted(lib.values()) == sorted([ele_a, ele_b])
    assert all(n.startswith("tool_") and n.split("#")[1] in ("LTR/Gypsy", "LTR/Copia") for n in lib)


@pytest.mark.skipif(shutil.which("cd-hit-est") is None, reason="cd-hit-est not on PATH")
def test_intact_passes_overlong_calls_through(tmp_path, monkeypatch):
    monkeypatch.setattr(bl, "CDHIT_MAX_SEQ", 2000)    # stand-in for cd-hit-est's 1 Mb limit
    rnd = random.Random(3)
    ele = "".join(rnd.choice("ACGT") for _ in range(1500))
    big = "".join(rnd.choice("ACGT") for _ in range(3000))
    g = ele + "A" * 100 + ele + "A" * 100 + big + "A" * 100 + big
    gf = tmp_path / "g.fa"
    gf.write_text(">c0001\n" + g + "\n")
    (tmp_path / "g.fa.fai").write_text(f"c0001\t{len(g)}\t7\t{len(g)}\t{len(g) + 1}\n")
    spans = [(0, 1500), (1600, 3100), (3200, 6200), (6300, 9300)]
    calls = tmp_path / "calls.tsv"
    calls.write_text("chrom\tstart\tend\tclass\n" + "".join(f"c0001\t{s}\t{e}\tGypsy\n" for s, e in spans))
    out = tmp_path / "lib"
    bl.main(["--calls", str(calls), "--genome", str(gf), "--out", str(out), "--prefix", "tool", "-t", "1",
             "--intact", "0.95"])
    lib = read_fa(str(out) + ".fa")
    assert sorted(lib.values()) == sorted([ele, big, big])   # short copies reduced, overlong ones kept unreduced
    assert {"tool_3#LTR/Gypsy", "tool_4#LTR/Gypsy"} <= set(lib)


def test_native_keeps_ltr_entries_only(tmp_path):
    nat = tmp_path / "nat.fa"
    nat.write_text(">fam1#LTR/Gypsy\nACGTACGT\n>fam2#DNA/hAT\nACGT\n>fam3#LINE/L1\nACGT\n>fam4#LTR/Copia x\nacgt\n")
    out = tmp_path / "n"
    r = subprocess.run([sys.executable, S, "--native", str(nat), "--out", str(out), "--prefix", "t"],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    lib = read_fa(str(out) + ".fa")
    assert lib == {"t_fam1#LTR/Gypsy": "ACGTACGT", "t_fam4#LTR/Copia": "ACGT"}
