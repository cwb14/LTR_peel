#!/usr/bin/env bash
# Build compiled tools from COPIES of the clones (clones stay pristine) and pull containers.
# usage: setsid -f nohup bash build_tools.sh > ../logs/build_tools.log 2>&1 < /dev/null
set -uo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
BENCH="${BENCH_ROOT:?set BENCH_ROOT to the benchmark working directory}"; ROOT=$(dirname "$BENCH")   # tool clones (LTR_peel, PrinTE, HiRepEat, HiTE, ...) sit next to the working directory
B=$BENCH/build
C=$BENCH/containers
export TMPDIR=$BENCH/tmp/build; mkdir -p "$TMPDIR" "$B" "$C"
log() { echo "$(date '+%F %T') $*"; }

copy() { [ -d "$B/$1" ] || cp -a "$ROOT/$1" "$B/$1"; }

# MaskerAide (Rust; HDF5 from the env)
( copy MaskerAide; conda activate maskeraide
  cd $B/MaskerAide && export CARGO_HOME=$B/cargo_home && \
  HDF5_DIR=$CONDA_PREFIX RUSTFLAGS="-C link-args=-Wl,-rpath,$CONDA_PREFIX/lib" cargo build --release > $B/maskeraide_build.log 2>&1 \
    && log "MaskerAide built: $(ls target/release/ | tr '\n' ' ')" || log "MaskerAide FAILED (see build/maskeraide_build.log)" ) &

# Look4LTRs (cmake from bench_tools, system g++)
( copy Look4LTRs; conda activate bench_tools
  cd $B/Look4LTRs && mkdir -p bin && cd bin && CC=/usr/bin/gcc CXX=/usr/bin/g++ cmake .. > $B/look4ltrs_build.log 2>&1 && \
  make -j 8 look4ltrs >> $B/look4ltrs_build.log 2>&1 && log "Look4LTRs built: $(ls $B/Look4LTRs/bin/look4ltrs)" \
    || log "Look4LTRs FAILED (see build/look4ltrs_build.log)" ) &

# HiRepEat (make inside the copy, in its env)
( copy HiRepEat; conda activate hirepeat_bench
  cd $B/HiRepEat && make > $B/hirepeat_build.log 2>&1 && log "HiRepEat built" || log "HiRepEat FAILED (see build/hirepeat_build.log)" ) &

# Containers (cache and temp under BENCH_ROOT, not the system disk)
( conda activate apptainer
  export APPTAINER_CACHEDIR=$BENCH/tmp/apptainer_cache APPTAINER_TMPDIR=$BENCH/tmp/apptainer_tmp
  mkdir -p $APPTAINER_CACHEDIR $APPTAINER_TMPDIR
  [ -s $C/hite_3.3.3.sif ] || apptainer pull $C/hite_3.3.3.sif docker://kanghu/hite:3.3.3 > $B/hite_pull.log 2>&1 \
     && log "HiTE image ok" || log "HiTE pull FAILED"
  [ -s $C/ltrpred.sif ] || apptainer pull $C/ltrpred.sif docker://drostlab/ltrpred > $B/ltrpred_pull.log 2>&1 \
     && log "LTRpred image ok" || log "LTRpred pull FAILED" ) &

wait
log "build_tools done"
