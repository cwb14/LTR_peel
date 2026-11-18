#!/usr/bin/env bash
# Re-run the evaluation workflow every 30 min while tool runs finish; final pass after the tool queue ends.
# usage: setsid -f nohup bash scripts/eval_loop.sh [cores] > logs/eval_loop.log 2>&1 < /dev/null
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate bench_tools
B="${BENCH_ROOT:?set BENCH_ROOT to the benchmark working directory}"; cd $B
CORES=${1:-32}; export TMPDIR=$B/tmp/smk; mkdir -p $TMPDIR
while true; do
  finished=0; grep -q "all jobs finished" logs/queue_main.log && finished=1
  echo "$(date '+%F %T') snakemake pass (queue finished=$finished)"
  snakemake -s workflow/Snakefile --cores $CORES --keep-going --rerun-incomplete --default-resources tmpdir=$TMPDIR \
    > logs/snakemake_last.log 2>&1
  rc=$?                                   # capture before any $(...) resets $?
  nerr=$(grep -cE 'Error in rule|^[A-Za-z]*(Exception|Error):' logs/snakemake_last.log)
  echo "$(date '+%F %T') pass exit $rc : $(grep -c 'Finished jobid' logs/snakemake_last.log) jobs; $nerr errors"
  grep -A3 "Error in rule" logs/snakemake_last.log | grep -E "rule|jobid|output" | head -20
  grep -E -A2 '^[A-Za-z]*(Exception|Error):' logs/snakemake_last.log | head -8
  [ $finished = 1 ] && break
  sleep 1800
done
echo "$(date '+%F %T') eval loop done"
