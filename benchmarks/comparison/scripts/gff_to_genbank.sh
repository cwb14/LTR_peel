#!/usr/bin/env bash
# Rename a RefSeq GFF's seqids (NC_/NW_) to the paired GenBank accessions (the names data/<ds>/names.tsv knows), using
# the assembly report's GenBank-Accn (col 5) / RefSeq-Accn (col 7) columns; sequences without a pair (mito) are dropped.
# usage: gff_to_genbank.sh ASSEMBLY_REPORT GFF.gz OUT.gff3.gz
set -euo pipefail
[ $# -eq 3 ] || { echo "usage: $0 ASSEMBLY_REPORT GFF.gz OUT.gff3.gz" >&2; exit 1; }
awk -F'\t' -v OFS='\t' 'NR==FNR { if ($0 !~ /^#/ && $7 != "na" && $5 != "na") m[$7] = $5; next }
    /^#/ { print; next } ($1 in m) { $1 = m[$1]; print }' "$1" <(zcat "$2") | gzip > "$3.part"
mv "$3.part" "$3"
echo "gff_to_genbank: $3 ($(zcat "$3" | grep -vc '^#') features)" >&2
