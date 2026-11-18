#!/usr/bin/env bash
# Run one masking job outside the workflow, with the exact command of Snakefile rule mask (for DEFERRED targets).
# usage: setsid -f nohup nice -n 10 bash scripts/mask_deferred.sh DS LIB > mask/DS/LIB/deferred.log 2>&1 < /dev/null
set -euo pipefail
B="${BENCH_ROOT:?set BENCH_ROOT to the benchmark working directory}"; ds=$1; lib=$2; threads=16
[ -s $B/libs/$ds/$lib.fa ] || { echo "error: no library $B/libs/$ds/$lib.fa" >&2; exit 1; }
export RAYON_NUM_THREADS=$threads OMP_NUM_THREADS=$threads MASKERAIDE_HOME=$B/tmp/maskeraide_home
export TMPDIR=$B/tmp/mask_${ds}_${lib}; mkdir -p $TMPDIR $B/mask/$ds/$lib; cd $B/mask/$ds/$lib
rm -rf RM_*                                     # scratch of an earlier, stopped run
ln -sf $B/data/$ds/genome.fa $ds.fa
echo "$(date '+%F %T') start rmfast $ds $lib"
$B/build/MaskerAide/target/release/rmfast -pa $threads -q -no_is -norna -nolow -div 40 -cutoff 225 \
  -lib $B/libs/$ds/$lib.fa -dir . $ds.fa > rmfast.log 2>&1
echo "$(date '+%F %T') done rmfast $ds $lib rc=$?"
