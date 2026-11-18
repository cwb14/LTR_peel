#!/usr/bin/env bash
# HiTE via its README's recommended conda route (the 3.3.3 image lacks the LtrDetector binary -> zero LTR candidates).
# Replicates configure.py without editing ~/.bashrc: chmod tools, recompile LtrDetector; PATH is set in each cmd.sh.
set -uo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
B="${BENCH_ROOT:?set BENCH_ROOT to the benchmark working directory}"; ROOT=$(dirname "$B")
export TMPDIR=$B/tmp/envs; mkdir -p $TMPDIR
[ -d $B/build/HiTE ] || cp -a $ROOT/HiTE $B/build/HiTE
conda env list | awk '{print $1}' | grep -qx hite_bench || mamba env create -y -n hite_bench -f $ROOT/HiTE/environment.yml
conda env export -n hite_bench > $B/envs/hite_bench.yml
conda activate hite_bench
cd $B/build/HiTE
chmod +x tools/* module/* bin/NeuralTE/tools/* bin/LTR_FINDER_parallel-master/bin/LTR_FINDER.x86_64-1.0.7/ltr_finder bin/FiLTR-main/tools/* 2>/dev/null
cd bin/FiLTR-main/bin/LtrDetector && rm -rf bin && cd src && make bin && make tr -j 8
ls -la $B/build/HiTE/bin/FiLTR-main/bin/LtrDetector/bin/LtrDetector && echo "HiTE build ok"
