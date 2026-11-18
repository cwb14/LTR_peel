#!/usr/bin/env bash
# Approach 1b: one dataset per seed x timepoint: data/SD<rep>_<kk> (prep_genome.py --split) + sims/truth/SD<rep>_<kk>.
# kk 00 = burn-in; kk = k -> gen (k x 1,540,000). usage: bash sims/div/prep_div.sh  (from bench/)
set -euo pipefail
B="${BENCH_ROOT:?set BENCH_ROOT to the benchmark working directory}"; PY="$(conda info --base)/envs/bench_tools/bin/python -I"
one() {
  rep=$1; k=$2; ds=$(printf "SD%d_%02d" $rep $k)
  if [ $k -eq 0 ]; then snap=burnin; else snap=gen$((k * 1540000))_final; fi
  $PY $B/scripts/prep_genome.py $B/sims/div/D$rep/$snap.fasta $B/data/$ds --split > $B/logs/prep_$ds.log 2>&1
  $PY $B/scripts/printe_truth.py -v $B/sims/div/D$rep/$snap.bed $B/sims/truth/$ds >> $B/logs/prep_$ds.log 2>&1
  echo "$ds ok"
}
export -f one; export B PY
for rep in 1 2 3; do for k in $(seq 0 10); do echo "$rep $k"; done; done | xargs -P 8 -n 2 bash -c 'one "$@"' _
