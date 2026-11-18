#!/usr/bin/env bash
# Normalize all real genomes (prep_genome.py) in parallel; one log per dataset in logs/prep_<ds>.log
set -uo pipefail
B="${BENCH_ROOT:?set BENCH_ROOT to the benchmark working directory}"
PY=$(conda info --base)/envs/bench_tools/bin/python
declare -A IN=(
  [rice]=$B/data/rice/raw/MSU7_all.con.fa
  [arab]=$B/data/arab/raw/GCF_000001735.4_TAIR10.1_genomic.fna.gz
  [maize]=$B/data/maize/raw/Zm-B73-REFERENCE-NAM-5.0.fa.gz
  [finch]=$B/data/finch/raw/GCA_048771995.1_bTaeGut7.mat_genomic.fna.gz
  [gorilla]=$B/data/gorilla/raw/GCA_029281585.3_NHGRI_mGorGor1-v2.1_pri_genomic.fna.gz
  [asue]=$B/data/asue/raw/Arabidopsis_suecica.4x.GCA_905175345.1.fa.gz
  [pann]=$B/data/pann/raw/Poa_annua.4x.GCA_026929855.1.fa.gz
)
for ds in "${!IN[@]}"; do
  [ -s $B/data/$ds/genome.fa.fai ] && { echo "skip $ds"; continue; }
  $PY -I $B/scripts/prep_genome.py "${IN[$ds]}" $B/data/$ds --split > $B/logs/prep_$ds.log 2>&1 &
done
wait
echo "prep_all done $(date '+%F %T')"
