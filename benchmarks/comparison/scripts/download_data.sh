#!/usr/bin/env bash
# Download benchmark genomes, curated libraries and related-species proteomes into bench/data/<ds>/raw/.
# NCBI files are checked against the assembly's md5checksums.txt. Re-running skips files already verified.
# usage: setsid -f nohup bash download_data.sh > ../logs/download_data.log 2>&1 < /dev/null
set -uo pipefail
D="${BENCH_ROOT:?set BENCH_ROOT to the benchmark working directory}"/data
NCBI=https://ftp.ncbi.nlm.nih.gov/genomes/all
log() { echo "$(date '+%F %T') $*"; }

fetch() {  # fetch URL OUT
    local url=$1 out=$2
    [ -s "$out" ] && { log "have $out"; return 0; }
    mkdir -p "$(dirname "$out")"
    log "get $url"
    curl -fsSL --retry 5 --retry-delay 20 --max-time 7200 -o "$out.part" "$url" && mv "$out.part" "$out" \
        || { log "FAILED $url"; rm -f "$out.part"; return 1; }
}

ncbi() {  # ncbi ACC_NAME SUFFIX OUTDIR  (e.g. GCF_000231095.2_ObraRS2 protein.faa.gz)
    local an=$1 suf=$2 od=$3 acc=${1%%_[A-Za-z]*}
    local p1=${an:0:3} n=${an:4:9}
    local base=$NCBI/$p1/${n:0:3}/${n:3:3}/${n:6:3}/$an
    fetch "$base/${an}_$suf" "$od/${an}_$suf" || return 1
    fetch "$base/md5checksums.txt" "$od/${an}.md5checksums.txt" || return 1
    local want; want=$(awk -v f="./${an}_$suf" '$2==f{print $1}' "$od/${an}.md5checksums.txt")
    local got; got=$(md5sum "$od/${an}_$suf" | cut -d' ' -f1)
    if [ "$want" = "$got" ]; then log "md5 ok ${an}_$suf"; else log "MD5 MISMATCH ${an}_$suf want=$want got=$got"; return 1; fi
}

gdrive() { fetch "https://drive.usercontent.google.com/download?id=$1&export=download" "$2"; }

# genomes
fetch https://rice.uga.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_7.0/all.dir/all.con \
      $D/rice/raw/MSU7_all.con.fa
ncbi GCF_000001735.4_TAIR10.1 genomic.fna.gz $D/arab/raw
ncbi GCA_048771995.1_bTaeGut7.mat genomic.fna.gz $D/finch/raw
ncbi GCA_029281585.3_NHGRI_mGorGor1-v2.1_pri genomic.fna.gz $D/gorilla/raw

# curated libraries (HiRepEat preprint, tinyurl.com/HiRepEat-curated -> Google Drive)
gdrive 1_1gQ9Z79fXofMekEDRNEP8jMlxDzMtdc $D/finch/raw/taeniopygia_guttata_curated.fa
gdrive 1kYOXNWE67JmYRrrHYUkOoJ3j3IqKNmgo $D/gorilla/raw/gorilla_gorilla_curated_lib.fa

# related-species proteomes for LTRquest --proteins (closest RefSeq relative, never the target)
ncbi GCF_000231095.2_ObraRS2 protein.faa.gz $D/proteins          # rice
ncbi GCF_000003195.3_Sorghum_bicolor_NCBIv3 protein.faa.gz $D/proteins   # maize
ncbi GCF_022539315.1_serCan2020 protein.faa.gz $D/proteins       # zebra finch
ncbi GCF_000001405.40_GRCh38.p14 protein.faa.gz $D/proteins      # gorilla
ncbi GCF_019359855.2_Kyuss_2.0 protein.faa.gz $D/proteins        # Poa annua
ncbi GCF_000004255.2_v.1.0 protein.faa.gz $D/proteins            # Arabidopsis thaliana, A. suecica
log "all done"

# gene annotations for approach 4: RefSeq annotation of the same assemblies (seqids converted to GenBank
# names by scripts/gff_to_genbank.sh); maize Zm00001eb.1 from a local copy of the MaizeGDB release (set MAIZE_GENES_GFF)
ncbi GCF_048771995.1_bTaeGut7.mat genomic.gff.gz $D/finch/genes/raw
ncbi GCF_048771995.1_bTaeGut7.mat assembly_report.txt $D/finch/genes/raw
ncbi GCF_029281585.2_NHGRI_mGorGor1-v2.1_pri genomic.gff.gz $D/gorilla/genes/raw
ncbi GCF_029281585.2_NHGRI_mGorGor1-v2.1_pri assembly_report.txt $D/gorilla/genes/raw
[ -s $D/maize/genes/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz ] || { mkdir -p $D/maize/genes && \
    cp "${MAIZE_GENES_GFF:?set MAIZE_GENES_GFF to a local copy of Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz (MaizeGDB)}" $D/maize/genes/; }
