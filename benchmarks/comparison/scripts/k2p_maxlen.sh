#!/usr/bin/env bash
# Uniform K2P (Kmer2LTR via phase_inputs.py) for one call set with spans > 100 kb left out (status too_long, K2P NA).
# For HiRepEat call sets with Mb-long spans (finch failed twice; gorilla up to 26.8 Mb). usage: k2p_maxlen.sh DS TOOL
set -euo pipefail
cd "${BENCH_ROOT:?set BENCH_ROOT to the benchmark working directory}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate kmer2ltr
export PATH=$(conda info --base)/envs/bench_tools/bin:$PATH
ds=$1; t=$2; o=eval/k2p/$ds/$t
$(conda info --base)/envs/bench_tools/bin/python -I scripts/phase_inputs.py --calls calls/$ds/$t.ids.tsv \
    --genome data/$ds/genome.fa --outdir $o -t 16 --max-len 100000
rm -f $o/elements.fa $o/elements.k.fa
echo "k2p_maxlen $ds $t done"
