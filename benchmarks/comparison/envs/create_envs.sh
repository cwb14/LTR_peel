#!/usr/bin/env bash
# Create the benchmark's conda environments, one after another (sequential: one mamba package cache).
# Each env is exported to envs/<name>.yml afterwards. Re-running skips envs that already exist.
# usage: setsid -f nohup bash create_envs.sh > ../logs/create_envs.log 2>&1 < /dev/null
set -uo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
BENCH="${BENCH_ROOT:?set BENCH_ROOT to the benchmark working directory}"; ROOT=$(dirname "$BENCH")   # tool clones (LTR_peel, PrinTE, HiRepEat, HiTE, ...) sit next to the working directory
ENVDIR=$BENCH/envs
export TMPDIR=$BENCH/tmp/envs; mkdir -p "$TMPDIR"
export PIP_NO_CACHE_DIR=1

have() { conda env list | awk '{print $1}' | grep -qx "$1"; }
done_msg() { echo "$(date '+%F %T') [$1] $2"; }
export_env() { conda env export -n "$1" > "$ENVDIR/$1.yml" 2>/dev/null; done_msg "$1" "exported"; }

mk() {   # mk NAME -- mamba create args...
    local n=$1; shift
    if have "$n"; then done_msg "$n" "exists, skipped"; return 0; fi
    done_msg "$n" "creating"
    if mamba create -y -n "$n" "$@"; then done_msg "$n" "created"; else done_msg "$n" "FAILED"; return 1; fi
}

# 1. LTR_peel from the clone (pinned commit listed on the benchmark page)
if ! have ltr_peel_bench; then
    done_msg ltr_peel_bench creating
    mamba env create -y -n ltr_peel_bench -f $ROOT/LTR_peel/environment.yml && \
    conda run -n ltr_peel_bench pip install $ROOT/LTR_peel && done_msg ltr_peel_bench created || done_msg ltr_peel_bench FAILED
fi
export_env ltr_peel_bench

# 2. PrinTE from the clone
if ! have printe_bench; then
    done_msg printe_bench creating
    mamba env create -y -n printe_bench -f $ROOT/PrinTE/environment.yml && \
    conda run -n printe_bench pip install $ROOT/PrinTE && done_msg printe_bench created || done_msg printe_bench FAILED
fi
export_env printe_bench

# 3. LTR_retriever deps per its README (step-by-step route) + LTRharvest + LTR_Finder
mk ltrr_bench -c conda-forge -c bioconda perl perl-text-soundex libstdcxx-ng cd-hit 'repeatmasker<4.1.5' rmblast \
    tesorter trf blast genometools-genometools=1.6.6 ltr_finder=1.07
export_env ltrr_bench

# 4. MaskerAide build env (Rust >= 1.85, HDF5)
mk maskeraide -c conda-forge 'rust>=1.85' hdf5 pkg-config c-compiler cxx-compiler make
export_env maskeraide

# 5. general bench tools: scoring, lib-test.py (pyfastx, bitarray), EDTA cleanup scripts (perl, blast, trf), cmake
mk bench_tools -c conda-forge -c bioconda python=3.11 pyfastx bitarray numpy pandas scipy pytest cmake samtools \
    seqkit bedtools trf blast cd-hit perl edlib matplotlib-base
export_env bench_tools

# 6. HiRepEat: python 3.11 + its requirements, TEsorter 1.5.1 / hmmer 3.4 / rmblast as in its container env
if ! have hirepeat_bench; then
    mk hirepeat_bench -c conda-forge -c bioconda python=3.11 pip tesorter=1.5.1 hmmer=3.4 rmblast c-compiler make zlib && \
    conda run -n hirepeat_bench pip install -r $ROOT/HiRepEat/requirements.txt && done_msg hirepeat_bench "pip ok"
fi
export_env hirepeat_bench

# 7. DANTE + DANTE_LTR (petrnovak channel; dante_ltr pinned to the clone's version)
mk dante_bench -c bioconda -c conda-forge -c petrnovak dante_ltr=0.6.5.0 dante=0.2.11
export_env dante_bench

done_msg all "finished"
