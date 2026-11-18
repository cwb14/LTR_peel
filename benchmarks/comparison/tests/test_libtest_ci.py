"""libtest_ci.py on a tiny synthetic benchmark tree: exact lib-test LTR counts, the lib-test class regex, the totals
check, and paired differences (a library against itself has dF1 = 0 with a zero-width CI)."""
import csv
import os
import subprocess
import sys

S = os.path.join(os.path.dirname(__file__), "..", "scripts", "libtest_ci.py")
HDR = "   SW  perc perc perc  query      position in query           matching       repeat              position\n\n"


def out_line(chrom, s1, e, cls, name="x"):
    return f"  100 1.0 0.0 0.0 {chrom} {s1} {e} (0) + {name} {cls} 1 10 (0) 1\n"


def write(p, text):
    os.makedirs(os.path.dirname(p), exist_ok=True)
    with open(p, "w") as f:
        f.write(text)


def libtest_row(tp, fp, fn):
    return ("label\t#Metrics\tsens\tspec\taccu\tprec\tFDR\tF1\tTP\tTN\tFP\tFN\n"
            f"d/x\tr\t0\t0\t0\t0\t0\t0\t{tp}\t0\t{fp}\t{fn}\n")


def test_counts_regex_check_and_paired_delta(tmp_path):
    b = str(tmp_path)
    write(f"{b}/data/d/genome.fa.fai", "c1\t100\t4\t60\t61\nc2\t50\t110\t60\t61\n")
    # reference: c1 11-40 LTR, 61-70 DNA (not LTR), c2 1-10 'nonLTR' (lib-test's regex matches 'LTR' inside it)
    write(f"{b}/mask/d/curated/d.fa.out", HDR + out_line("c1", 11, 40, "LTR/Gypsy") + out_line("c1", 61, 70, "DNA/hAT")
          + out_line("c2", 1, 10, "nonLTR/x"))
    # test A: c1 21-50 (TP 20, FP 10), c2 nothing -> FN 10 (c1 11-20) + 10 (c2) = 20
    write(f"{b}/mask/d/A/d.fa.out", HDR + out_line("c1", 21, 50, "LTR/unknown"))
    write(f"{b}/eval/lib/d/A.tsv", libtest_row(20, 10, 20))
    # test B: identical to A but lib-test totals recorded differently -> MISMATCH
    write(f"{b}/mask/d/B/d.fa.out", HDR + out_line("c1", 21, 50, "LTR/unknown"))
    write(f"{b}/eval/lib/d/B.tsv", libtest_row(1, 1, 1))
    o = f"{b}/ci.tsv"
    r = subprocess.run([sys.executable, S, "d", "--base", b, "--block", "25", "--reps", "200", "--ref-lib", "A",
                        "-o", o], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    rows = {x["lib"]: x for x in csv.DictReader(open(o), delimiter="\t")}
    a = rows["A"]
    assert (int(a["TP"]), int(a["FP"]), int(a["FN"])) == (20, 10, 20)
    assert abs(float(a["F1"]) - 40 / 70) < 1e-4
    assert a["check"] == "ok" and rows["B"]["check"].startswith("MISMATCH")
    assert float(a["F1_lo"]) <= float(a["F1"]) <= float(a["F1_hi"])
    assert float(a["dF1"]) == 0 and float(a["dF1_lo"]) == 0 and float(a["dF1_hi"]) == 0
    assert int(a["n_blocks"]) == 4 + 2                    # 100 bp -> 4 blocks of 25, 50 bp -> 2
