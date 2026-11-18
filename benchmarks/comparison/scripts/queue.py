#!/usr/bin/env python3
"""Run jobs.tsv (make_jobs.py) under a thread budget, each job through bench_run.py (resource accounting + limits).

Resumable: skips jobs whose outdir has .done (and .failed unless --rerun-failed). Order: smallest genome first, so
early results arrive early and the largest (gorilla) jobs run last. A job starts when its threads fit the budget;
a job wider than the whole budget runs alone.
Launch detached:  setsid -f nohup python queue.py jobs.tsv --budget 80 > queue.log 2>&1 < /dev/null
"""
import argparse
import csv
import os
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent


def log(msg):
    print(time.strftime("%Y-%m-%d %H:%M:%S"), msg, flush=True)


def live_bench_runs():
    """{resolved outdir: pid} of every live bench_run.py process (first positional argument = outdir)."""
    out = {}
    for d in os.listdir("/proc"):
        if not d.isdigit():
            continue
        try:
            with open(f"/proc/{d}/cmdline", "rb") as f:
                argv = f.read().split(b"\0")
        except OSError:
            continue
        for i, x in enumerate(argv):
            if x.endswith(b"bench_run.py") and i + 1 < len(argv):
                out[str(Path(argv[i + 1].decode()).resolve())] = int(d)
                break
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description="Run jobs.tsv under a thread budget via bench_run.py.")
    ap.add_argument("jobs", help="jobs.tsv from make_jobs.py")
    ap.add_argument("--budget", type=int, default=80, help="total threads (default 80)")
    ap.add_argument("--cap-h", type=float, default=36.0, help="per-job wall limit, hours (default 36)")
    ap.add_argument("--mem-gb", type=float, default=200.0, help="per-job memory guard, GB PSS (default 200)")
    ap.add_argument("--rerun-failed", action="store_true", help="rerun jobs marked .failed")
    ap.add_argument("--only", nargs="*", default=None, help="restrict to these job_ids")
    ap.add_argument("--poll", type=float, default=20.0, help="seconds between scheduler checks (default 20)")
    a = ap.parse_args(argv)
    started_ids = set()
    adopted = {}          # job_id -> (pid, threads, outdir): bench_run processes left by an earlier queue

    def load_todo():
        """Re-read jobs.tsv every cycle, so jobs can be appended while the queue runs (one shared budget).
        A job whose bench_run.py is still alive (e.g. the previous queue died) is adopted, never started twice."""
        jobs = list(csv.DictReader(open(a.jobs), delimiter="\t"))
        live = live_bench_runs()
        out = []
        for j in jobs:
            od = Path(j["outdir"])
            if j["job_id"] in started_ids or j["job_id"] in adopted or \
                    (a.only is not None and j["job_id"] not in a.only):
                continue
            if str(od.resolve()) in live:
                adopted[j["job_id"]] = (live[str(od.resolve())], int(j["threads"]), od)
                log(f"adopted running {j['job_id']} (pid {adopted[j['job_id']][0]})")
                continue
            if (od / ".done").exists() or ((od / ".failed").exists() and not a.rerun_failed):
                continue
            out.append(j)
        out.sort(key=lambda j: (int(j["size"]), j["tool"]))
        return out, len(jobs)

    todo, n_jobs = load_todo()
    log(f"{len(todo)} of {n_jobs} jobs to run; {len(adopted)} adopted; budget {a.budget} threads; cap {a.cap_h} h; "
        f"mem {a.mem_gb} GB")
    running = {}
    while todo or running or adopted:
        for jid, (pid, th, od) in list(adopted.items()):
            if not os.path.exists(f"/proc/{pid}"):
                ok = (od / ".done").exists()
                log(f"{'done' if ok else 'FAILED'} {jid} (adopted)")
                del adopted[jid]
        new, _ = load_todo()
        known = {j["job_id"] for j in todo}
        for j in new:
            if j["job_id"] not in known:
                todo.append(j)
                log(f"queued new job {j['job_id']}")
        todo.sort(key=lambda j: (int(j["size"]), j["tool"]))
        for jid, (p, th) in list(running.items()):
            rc = p.poll()
            if rc is not None:
                log(f"{'done' if rc == 0 else 'FAILED'} {jid} rc={rc}")
                del running[jid]
        used = sum(th for _, th in running.values()) + sum(th for _, th, _ in adopted.values())
        started = True
        while todo and started:
            started = False
            for i, j in enumerate(todo):
                th = int(j["threads"])
                if used + th <= a.budget or not (running or adopted):
                    od = j["outdir"]
                    with open(os.path.join(od, "bench_run.log"), "ab") as lf:
                        p = subprocess.Popen([sys.executable, str(HERE / "bench_run.py"), od, "--cap-h", str(a.cap_h),
                                              "--mem-gb", str(a.mem_gb), "--threads", str(th), "--job-id", j["job_id"]],
                                             stdout=lf, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL)
                    running[j["job_id"]] = (p, th)
                    started_ids.add(j["job_id"])
                    used += th
                    log(f"start {j['job_id']} ({th} threads; {used}/{a.budget} in use)")
                    del todo[i]
                    started = True
                    break
        time.sleep(a.poll)
    log("all jobs finished")


if __name__ == "__main__":
    main()
