#!/usr/bin/env python3
"""Run one benchmark job and measure it.

Runs `bash cmd.sh` (or the command after --) in OUTDIR as a new session, output to OUTDIR/run.log. Every --interval
seconds it walks the job's process tree in /proc (descendants by parent PID, plus any process in the job's session)
and records the summed RSS and PSS (PSS splits shared pages, so forked workers are not double-counted).

Writes OUTDIR/resources.json:
  wall_s, cpu_s (user+sys of all waited-for descendants: getrusage RUSAGE_CHILDREN), cpu_proc_s (sum over every PID
  seen of its last utime+stime; also counts orphaned children), peak_rss_gb / peak_pss_gb (process tree, summed),
  peak_single_rss_gb (largest single process, = /usr/bin/time 'Maximum resident set size'), max_procs, exit_code,
  status (ok | failed | timeout | mem_guard), load averages at start and end.
and OUTDIR/.done (status ok) or OUTDIR/.failed.

Limits: wall time > --cap-h hours -> status timeout; tree PSS > --mem-gb -> status mem_guard. Either kills the job
(SIGTERM to every process of the tree, SIGKILL after 60 s).
"""
import argparse
import json
import os
import resource
import signal
import socket
import subprocess
import sys
import time

PAGE = os.sysconf("SC_PAGE_SIZE")
TICK = os.sysconf("SC_CLK_TCK")


def proc_table():
    """{pid: (ppid, sid, utime+stime seconds)} for every readable process."""
    out = {}
    for d in os.listdir("/proc"):
        if not d.isdigit():
            continue
        try:
            with open(f"/proc/{d}/stat") as f:
                s = f.read()
        except OSError:
            continue
        r = s.rfind(")")                       # comm may contain spaces and ')'
        f2 = s[r + 2:].split()
        # fields after comm: state(0) ppid(1) pgrp(2) session(3) ... utime(11) stime(12)
        if f2 and f2[0] == "Z":                # zombies hold no memory and are not alive
            continue
        try:
            out[int(d)] = (int(f2[1]), int(f2[3]), (int(f2[11]) + int(f2[12])) / TICK)
        except (IndexError, ValueError):
            continue
    return out


def tree(root, sid, table):
    kids = {}
    for pid, (ppid, _, _) in table.items():
        kids.setdefault(ppid, []).append(pid)
    seen, stack = set(), [root]
    while stack:
        p = stack.pop()
        if p in seen or p not in table:
            continue
        seen.add(p)
        stack.extend(kids.get(p, []))
    seen.update(pid for pid, (_, s, _) in table.items() if s == sid)
    return seen


def mem_of(pid):
    """(rss_bytes, pss_bytes) of one process; (0, 0) if it vanished."""
    rss = pss = 0
    try:
        with open(f"/proc/{pid}/smaps_rollup") as f:
            for line in f:
                if line.startswith("Rss:"):
                    rss = int(line.split()[1]) * 1024
                elif line.startswith("Pss:"):
                    pss = int(line.split()[1]) * 1024
    except OSError:
        try:
            with open(f"/proc/{pid}/statm") as f:
                rss = pss = int(f.read().split()[1]) * PAGE
        except OSError:
            pass
    return rss, pss


def kill_tree(root, sid, grace=60):
    for sig in (signal.SIGTERM, signal.SIGKILL):
        pids = tree(root, sid, proc_table())
        for p in pids:
            try:
                os.kill(p, sig)
            except OSError:
                pass
        if sig == signal.SIGTERM:
            t0 = time.time()
            while time.time() - t0 < grace and tree(root, sid, proc_table()):
                time.sleep(2)


def main(argv=None):
    ap = argparse.ArgumentParser(description="Run one job with wall/CPU/memory accounting and limits.")
    ap.add_argument("outdir", help="job directory (cwd of the job; holds cmd.sh unless a command follows --)")
    ap.add_argument("--cap-h", type=float, default=36.0, help="wall-clock limit in hours (default 36)")
    ap.add_argument("--mem-gb", type=float, default=200.0, help="kill if the tree's summed PSS exceeds this (default 200)")
    ap.add_argument("--interval", type=float, default=2.0, help="sampling interval in seconds (default 2)")
    ap.add_argument("--job-id", default=None, help="label stored in resources.json (default: outdir name)")
    ap.add_argument("--threads", type=int, default=None, help="threads the job was given (stored only)")
    ap.add_argument("-v", "--verbose", action="store_true", help="log every sample to stderr")
    ap.epilog = "Command: everything after '--' (default: bash cmd.sh in OUTDIR)."
    argv = list(sys.argv[1:] if argv is None else argv)
    cmd = []
    if "--" in argv:
        i = argv.index("--")
        argv, cmd = argv[:i], argv[i + 1:]
    a = ap.parse_args(argv)
    od = os.path.abspath(a.outdir)
    if not cmd:
        if not os.path.exists(os.path.join(od, "cmd.sh")):
            sys.exit(f"error: no command given and {od}/cmd.sh missing")
        cmd = ["bash", "cmd.sh"]
    for m in (".done", ".failed"):
        if os.path.exists(os.path.join(od, m)):
            os.remove(os.path.join(od, m))
    job = a.job_id or os.path.basename(od)
    load0 = os.getloadavg()
    t0 = time.time()
    print(f"bench_run: start {job}: {' '.join(cmd)}", file=sys.stderr, flush=True)
    with open(os.path.join(od, "run.log"), "ab") as log:
        p = subprocess.Popen(cmd, cwd=od, stdout=log, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                             start_new_session=True)
    sid = p.pid
    peak_rss = peak_pss = 0
    max_procs = 0
    cpu_seen = {}
    status = None
    while True:
        rc = p.poll()
        table = proc_table()
        pids = tree(p.pid, sid, table)
        for q in pids:
            cpu_seen[q] = table[q][2] if q in table else cpu_seen.get(q, 0.0)
        if rc is not None and not pids:
            break
        rss = pss = 0
        for q in pids:
            r1, p1 = mem_of(q)
            rss += r1
            pss += p1
        peak_rss, peak_pss = max(peak_rss, rss), max(peak_pss, pss)
        max_procs = max(max_procs, len(pids))
        if a.verbose:
            print(f"  t={time.time() - t0:.0f}s procs={len(pids)} rss={rss / 1e9:.2f}G pss={pss / 1e9:.2f}G",
                  file=sys.stderr, flush=True)
        if status is None and time.time() - t0 > a.cap_h * 3600:
            status = "timeout"
            print(f"bench_run: {job} exceeded {a.cap_h} h; killing", file=sys.stderr, flush=True)
            kill_tree(p.pid, sid)
        elif status is None and pss > a.mem_gb * 1e9:
            status = "mem_guard"
            print(f"bench_run: {job} PSS {pss / 1e9:.1f} GB > {a.mem_gb} GB; killing", file=sys.stderr, flush=True)
            kill_tree(p.pid, sid)
        if rc is not None and status is not None:
            break
        time.sleep(a.interval)
    rc = p.wait()
    wall = time.time() - t0
    ru = resource.getrusage(resource.RUSAGE_CHILDREN)
    if status is None:
        status = "ok" if rc == 0 else "failed"
    res = dict(job=job, cmd=" ".join(cmd), host=socket.gethostname(), threads=a.threads,
               start=time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(t0)),
               end=time.strftime("%Y-%m-%dT%H:%M:%S"), wall_s=round(wall, 1),
               cpu_s=round(ru.ru_utime + ru.ru_stime, 1), cpu_user_s=round(ru.ru_utime, 1),
               cpu_sys_s=round(ru.ru_stime, 1), cpu_proc_s=round(sum(cpu_seen.values()), 1),
               peak_rss_gb=round(peak_rss / 1e9, 3), peak_pss_gb=round(peak_pss / 1e9, 3),
               peak_single_rss_gb=round(ru.ru_maxrss * 1024 / 1e9, 3), max_procs=max_procs,
               exit_code=rc, status=status, cap_h=a.cap_h, mem_gb=a.mem_gb,
               load_start=[round(x, 1) for x in load0], load_end=[round(x, 1) for x in os.getloadavg()])
    with open(os.path.join(od, "resources.json"), "w") as f:
        json.dump(res, f, indent=1)
    with open(os.path.join(od, ".done" if status == "ok" else ".failed"), "w") as f:
        f.write(f"{status} {rc}\n")
    print(f"bench_run: end {job} status={status} rc={rc} wall={wall / 3600:.2f} h cpu={res['cpu_s'] / 3600:.2f} h "
          f"pss={res['peak_pss_gb']} GB", file=sys.stderr, flush=True)
    return 0 if status == "ok" else 1


if __name__ == "__main__":
    sys.exit(main())
