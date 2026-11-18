#!/usr/bin/env bash
# Approach 1b: generate the 3 divergence-series genomes from params.tsv (run dirs sims/div/D<rep>/), in parallel.
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate printe_bench
S="${BENCH_ROOT:?set BENCH_ROOT to the benchmark working directory}"/sims; export TMPDIR=$S/div/tmp; mkdir -p $TMPDIR
rep=0
grep -v '^#' $S/div/params.tsv | tail -n +2 | while IFS=$'\t' read set seed args; do
  rep=$((rep+1)); d=$S/div/D${rep}; mkdir -p $d
  ( cd $d && [ -s gen15400000_final.bed ] || { timeout 4h printe $args -s $seed -cl $S/lib_clean.fa > printe.out 2>&1;
      echo "$(date '+%F %T') D${rep} seed $seed exit $?"; } ) &
done
wait
