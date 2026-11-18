import json
import os
import subprocess
import sys

SCRIPT = os.path.join(os.path.dirname(__file__), "..", "scripts", "bench_run.py")

# two child processes, each holding ~200 MB and burning ~1 s CPU, plus a grandchild
WORK = r"""
import os, sys, time, subprocess
code = "b = bytearray(200*1024*1024); [b.__setitem__(i, 1) for i in range(0, len(b), 4096)]; t=time.process_time()\nwhile time.process_time()-t < 1.0: pass\ntime.sleep(3)"
ps = [subprocess.Popen([sys.executable, '-c', 'import time\n' + code]) for _ in range(2)]
[p.wait() for p in ps]
"""


def test_tree_memory_cpu(tmp_path):
    (tmp_path / "w.py").write_text(WORK)
    r = subprocess.run([sys.executable, SCRIPT, str(tmp_path), "--interval", "0.5", "--",
                        sys.executable, "w.py"], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    res = json.loads((tmp_path / "resources.json").read_text())
    assert res["status"] == "ok" and (tmp_path / ".done").exists()
    assert 0.35 < res["peak_pss_gb"] < 1.0, res          # two ~0.21 GB children summed
    assert 0.18 < res["peak_single_rss_gb"] < 0.4, res   # largest single process
    assert res["cpu_s"] >= 1.8, res                       # 2 x 1 s busy loop
    assert res["cpu_proc_s"] >= 1.5, res


def test_timeout_kills(tmp_path):
    r = subprocess.run([sys.executable, SCRIPT, str(tmp_path), "--cap-h", str(2 / 3600), "--interval", "0.5", "--",
                        "bash", "-c", "sleep 300 & sleep 300; wait"], capture_output=True, text=True, timeout=120)
    res = json.loads((tmp_path / "resources.json").read_text())
    assert res["status"] == "timeout" and (tmp_path / ".failed").exists()
    assert res["wall_s"] < 90


def test_failure_marked(tmp_path):
    (tmp_path / "cmd.sh").write_text("exit 3\n")
    r = subprocess.run([sys.executable, SCRIPT, str(tmp_path)], capture_output=True, text=True)
    res = json.loads((tmp_path / "resources.json").read_text())
    assert r.returncode == 1 and res["status"] == "failed" and res["exit_code"] == 3
