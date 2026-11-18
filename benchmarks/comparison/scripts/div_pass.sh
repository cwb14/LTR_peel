#!/usr/bin/env bash
# Approach 1b (post hoc) evaluation loop: every 20 min, score finished divergence-series runs (runs/SD*/<tool>/.done).
#   1. calls -> score_elements.py through the Snakefile's generic rules (calls, calls_derived, sim_score), targets only
#      (outputs disjoint from every other pass, so --nolock);
#   2. PrinTE's own evaluator (bedtools.py --printe-intact -r 0.9) on the same calls, run directly: rule sim_printe
#      hard-codes sims/<ds>/gen1000000_final.bed (S1/S2 layout).
# Ends after one last pass once the div queue (scripts/queue.py runs/jobs_div.tsv) has exited.
# usage: setsid -f nohup bash scripts/div_pass.sh > logs/div_pass.log 2>&1 < /dev/null
set -uo pipefail
cd "${BENCH_ROOT:?set BENCH_ROOT to the benchmark working directory}"
SM="nice -n 10 $(conda info --base)/envs/bench_tools/bin/snakemake -s workflow/Snakefile"
OPT="--cores 8 --nolock --keep-going --rerun-incomplete --default-resources tmpdir=$PWD/tmp/smk"
log() { echo "$(date '+%F %T') div_pass: $*"; }
snap() { local k=$((10#${1#*_})); if [ $k -eq 0 ]; then echo burnin; else echo gen$((k * 1540000))_final; fi; }
final=0
while :; do
    pgrep -f "scripts/queue[.]py runs/jobs_div" > /dev/null || final=1
    T=()
    for d in runs/SD*/*/; do
        d=${d%/}; ds=$(basename "$(dirname "$d")"); tool=$(basename "$d")
        [ -e "$d/.done" ] || continue
        [ -s "eval/sim/$ds/$tool.metrics.tsv" ] || T+=("eval/sim/$ds/$tool.metrics.tsv")
        if [ "$tool" = ltrquest ] && [ ! -s "eval/sim/$ds/ltrquest_hc.metrics.tsv" ]; then
            T+=("eval/sim/$ds/ltrquest_hc.metrics.tsv")
        fi
    done
    if [ ${#T[@]} -gt 0 ]; then
        log "${#T[@]} targets"
        $SM "${T[@]}" $OPT >> logs/div_pass_smk.log 2>&1
        log "snakemake exit $?"
    fi
    n=0
    for c in calls/SD*/*.tsv; do
        [ -e "$c" ] || continue
        case $c in *.ids.tsv) continue ;; esac
        ds=$(basename "$(dirname "$c")"); tool=$(basename "$c" .tsv); o=eval/sim/$ds/$tool.printe.txt
        [ -s "$o" ] && continue
        bed=sims/div/D${ds:2:1}/$(snap "$ds").bed
        (
            source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate printe_bench
            awk -F'\t' 'NR>1{print $1":"$2+1"-"$3}' "$c" > "$o.pass"
            if [ -s "$o.pass" ]; then
                python -m printe.util.bedtools -pass_scn "$o.pass" -bed "$bed" --printe-intact -r 0.9 > "$o"
            else
                echo 'no calls' > "$o"
            fi
            rm -f "$o.pass"
        )
        n=$((n + 1))
    done
    [ $n -gt 0 ] && log "printe cross-check: $n"
    [ $final -eq 1 ] && { log "done"; break; }
    sleep 1200
done
