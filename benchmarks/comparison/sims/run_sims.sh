#!/usr/bin/env bash
# Generate the 6 benchmark genomes from params.tsv (run dirs sims/<set>_<rep>/), in parallel.
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate printe_bench
S="${BENCH_ROOT:?set BENCH_ROOT to the benchmark working directory}"/sims; export TMPDIR=$S/tmp; mkdir -p $TMPDIR
rep=0; prev=""
grep -v '^#' $S/params.tsv | tail -n +2 | while IFS=$'\t' read set seed args; do
  [ "$set" != "$prev" ] && rep=0; prev=$set; rep=$((rep+1)); d=$S/${set}_${rep}; mkdir -p $d
  ( cd $d && [ -s gen1000000_final.bed ] || { timeout 8h printe $args -m 1.3e-8 -ge 1000000 -st 250000 -t 8 -s $seed \
      -cl $S/lib_clean.fa > printe.out 2>&1; echo "$(date '+%F %T') ${set}_${rep} seed $seed exit $?"; } ) &
done
wait
