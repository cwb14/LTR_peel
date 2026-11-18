#!/usr/bin/env bash
# Pre-registered 36 h stop for one evaluation job, enforced on the server so it does not depend on a live session.
# Until DEADLINE: exit quietly as soon as ROOT (the job's snakemake) ends. At DEADLINE: SIGTERM every descendant of ROOT
# whose command line matches PATTERN (grep -E), never ROOT itself, so snakemake records the job as failed and exits
# cleanly (locks released); SIGKILL any that survive 60 s. The Snakefile's DEFERRED / NOT_RUN sets keep the main
# evaluation loop from restarting the stopped job.
# usage: setsid -f nohup bash scripts/deadline_stop.sh "YYYY-MM-DD HH:MM" ROOT_PID 'PATTERN' > logs/deadline_TAG.log 2>&1 < /dev/null
set -uo pipefail
DL=$(date -d "$1" +%s) || exit 1
ROOT=$2
PAT=$3
log() { echo "$(date '+%F %T') deadline_stop: $*"; }
descendants() {  # PIDs below $1, matching PAT
    ps -eo pid=,ppid= | awk -v r="$1" '{p[$1]=$2} END {for (x in p) {y=x; while (y in p && p[y] != 0) {y=p[y]; if (y==r) {print x; break}}}}' |
        while read -r d; do
            tr '\0' ' ' < "/proc/$d/cmdline" 2>/dev/null | grep -qE -- "$PAT" && echo "$d"
        done
}
log "armed: descendants of $ROOT matching '$PAT' at $1"
while [ "$(date +%s)" -lt "$DL" ]; do
    kill -0 "$ROOT" 2>/dev/null || { log "$ROOT ended before the deadline; nothing to do"; exit 0; }
    sleep 60
done
T=$(descendants "$ROOT" | tr '\n' ' ')
if [ -z "$T" ]; then
    log "deadline reached: no descendant of $ROOT matches '$PAT'; nothing to do"
    exit 0
fi
for d in $T; do
    log "deadline reached: SIGTERM $d: $(tr '\0' ' ' < /proc/$d/cmdline 2>/dev/null | cut -c1-200)"
    kill -TERM "$d" 2>/dev/null
done
sleep 60
for d in $T; do
    kill -0 "$d" 2>/dev/null && { kill -KILL "$d" 2>/dev/null; log "SIGKILL $d (alive 60 s after SIGTERM)"; }
done
log "done"
