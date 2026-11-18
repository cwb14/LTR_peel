#!/usr/bin/env bash
# finch HiRepEat: Kmer2LTR workers die on its Mb-long spans (failed twice); spans > 100 kb get K2P NA (too_long)
set -euo pipefail
cd "${BENCH_ROOT:?set BENCH_ROOT to the benchmark working directory}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate kmer2ltr
export PATH=$(conda info --base)/envs/bench_tools/bin:$PATH
ds=${1:-finch}; t=hirepeat; o=eval/k2p/$ds/$t
$(conda info --base)/envs/bench_tools/bin/python -I scripts/phase_inputs.py --calls calls/$ds/$t.ids.tsv --genome data/$ds/genome.fa --outdir $o -t 16 --max-len 100000
rm -f $o/elements.fa $o/elements.k.fa
$(conda info --base)/envs/bench_tools/bin/python -I scripts/adjudicate.py --ds $ds --tool $t --out eval/adj/$ds/$t.calls.tsv
echo "k2p+adjudicate $ds hirepeat done"
