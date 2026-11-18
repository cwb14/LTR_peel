#!/usr/bin/env bash
# Does the library-reduction step change library-masking scores?
# The reduction exists only to make masking tractable; the ideal is masking with all intact calls. Per tool, every
# variant starts from the same tandem-cleaned calls (EDTA protocol step 2 output, libs/<ds>/<tool>.uniform.work/raw.fa.cln)
# and differs only in the reduction:
#   all          no reduction (the reference the reductions should reproduce)
#   edta         cleanup_nested.pl (the benchmark library; its existing mask is reused)
#   edta_rep     cleanup_nested.pl again on the same input (run-to-run noise: Perl hash order + thread race)
#   cd95         cd-hit-est -c 0.95 only
#   cd95_nested  cd-hit-est -c 0.95, then cleanup_nested.pl (the maize/gorilla fallback)
# Each -> rmfast (EDTA params, as workflow rule mask) -> EDTA lib-test.py LTR row vs mask/<ds>/curated.
# usage: bash run.sh <ds> <threads> <rep:0|1> <tool> [<tool> ...]     -> experiments/libreduce/<ds>/results.tsv
set -euo pipefail
B="${BENCH_ROOT:?set BENCH_ROOT to the benchmark working directory}"
ROOT=$(dirname $B)
E=$B/experiments/libreduce
ds=$1; T=$2; REP=$3; shift 3
source "$(conda info --base)/etc/profile.d/conda.sh"
export TMPDIR=$B/tmp/libreduce_$ds; mkdir -p $TMPDIR
RES=$E/$ds/results.tsv; mkdir -p $E/$ds
[ -s $RES ] || printf "ds\ttool\tvariant\tn_seq\tlib_mb\treduce_s\tmask_cpu_s\tsens\tprec\tF1\tTP\tFP\tFN\n" > $RES

nested() {  # $1 input fasta, $2 output fasta, $3 workdir; prints seconds
  mkdir -p $3; cp $1 $3/in.fa
  local t0=$(date +%s)
  (cd $3 && perl $ROOT/EDTA/bin/cleanup_nested.pl -in in.fa -threads $T -minlen 80 -cov 0.95 > log 2>&1)
  awk '/^>/{print $1; next}{print}' $3/in.fa.cln > $2
  echo $(( $(date +%s) - t0 ))
}

for tool in "$@"; do
  W=$E/$ds/$tool; mkdir -p $W
  src=$B/libs/$ds/$tool.uniform.work/raw.fa.cln
  [ -s $src ] || { echo "skip $tool: no $src"; continue; }
  conda activate ltrr_bench
  declare -A RS=()
  [ -s $W/all.fa ] || awk '/^>/{print $1; next}{print}' $src > $W/all.fa; RS[all]=0
  [ -s $W/edta.fa ] || cp $B/libs/$ds/$tool.uniform.fa $W/edta.fa; RS[edta]=NA
  if [ ! -s $W/cd95.fa ]; then
    t0=$(date +%s); cd-hit-est -i $W/all.fa -o $W/cd95.fa -c 0.95 -T $T -M 0 -d 0 > $W/cd95.log 2>&1
    echo $(( $(date +%s) - t0 )) > $W/cd95.s
  fi; RS[cd95]=$(cat $W/cd95.s)
  if [ ! -s $W/cd95_nested.fa ]; then nested $W/cd95.fa $W/cd95_nested.fa $W/nest_cd95 > $W/cd95_nested.s; fi
  RS[cd95_nested]=$(( $(cat $W/cd95.s) + $(cat $W/cd95_nested.s) ))
  VARS="all edta cd95 cd95_nested"
  if [ "$REP" = 1 ]; then
    [ -s $W/edta_rep.fa ] || nested $W/all.fa $W/edta_rep.fa $W/nest_rep > $W/edta_rep.s
    RS[edta_rep]=$(cat $W/edta_rep.s); VARS="$VARS edta_rep"
  fi
  for v in $VARS; do
    M=$W/mask_$v; mkdir -p $M
    if [ ! -s $M/$ds.fa.out ]; then
      if [ $v = edta ] && [ -s $B/mask/$ds/$tool.uniform/$ds.fa.out ]; then
        cp $B/mask/$ds/$tool.uniform/{$ds.fa.out,rmfast.log} $M/
      else
        ( export RAYON_NUM_THREADS=$T OMP_NUM_THREADS=$T MASKERAIDE_HOME=$B/tmp/maskeraide_home
          cd $M && ln -sf $B/data/$ds/genome.fa $ds.fa && \
          $B/build/MaskerAide/target/release/rmfast -pa $T -q -no_is -norna -nolow -div 40 -cutoff 225 \
            -lib $W/$v.fa -dir . $ds.fa > rmfast.log 2>&1 )
      fi
    fi
    conda activate bench_tools
    [ -s $M/$ds.fa.out.repeatmasker.lib.report ] || \
      python $ROOT/EDTA/lib-test.py --genome $B/mask/$ds/curated/$ds.fa --reference $B/mask/$ds/curated/$ds.fa.out \
        --test $M/$ds.fa.out --extended_report > $M/libtest.log 2>&1
    conda activate ltrr_bench
    n=$(grep -c '>' $W/$v.fa); mb=$(awk '!/^>/{n+=length($0)}END{printf "%.2f", n/1e6}' $W/$v.fa)
    cpu=$(grep -oP 'engines total CPU \K[0-9.]+' $M/rmfast.log | tail -1)
    awk -v p="$ds\t$tool\t$v\t$n\t$mb\t${RS[$v]}\t${cpu:-NA}" 'BEGIN{OFS="\t"} /[.]ltr[.]lib[.]report/{print p,$2,$5,$7,$8,$10,$11}' \
      $M/$ds.fa.out.repeatmasker.lib.report >> $RES
  done
  unset RS
done
echo "libreduce $ds done $(date '+%F %T')"
