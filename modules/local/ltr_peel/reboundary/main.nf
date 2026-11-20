// Template re-boundarying. One pooled task: an element's templates may sit in any
// genome, so the samples cannot be split here the way the per-genome stages are. It rewrites the _clean_ tables and FASTAs and
// writes one <prefix>_reboundary.tsv per sample, which LTR_PEEL_GFF3 reads to
// follow the elements it renamed.

process LTR_PEEL_REBOUNDARY {
    tag "pooled"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container 'ghcr.io/cwb14/ltr_peel:1.0.1'

    input:
    val(ids)
    path(tables, stageAs: 'in/*')
    path(fastas, stageAs: 'in/*')
    path(genomes, stageAs: 'genomes/g??/*')

    output:
    // No tuple to wrap them, so no parentheses: a bare `path(...), emit:` parses
    // only under the strict syntax parser, and Nextflow 24.x does not have it.
    path "*_depth*_clean_ltr.tsv", arity: '1..*', emit: tsv
    path "*_depth*_clean_ltr.fa" , arity: '1..*', emit: fasta
    path "*_reboundary.tsv"      , arity: '1..*', emit: sidecar
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def genome_list = (genomes instanceof List ? genomes : [genomes]).join(' ')
    // Rewritten in place, so copied out of the staging directory first: editing
    // through Nextflow's input symlinks would corrupt the upstream task's outputs.
    """
    cp in/*_clean_ltr.tsv in/*_clean_ltr.fa .

    ltr_peel-reboundary \\
        --indir . \\
        --prefix ${ids.join(' ')} \\
        --genome ${genome_list} \\
        --threads ${task.cpus} \\
        --mutation-rate ${params.mutation_rate} \\
        --tools-dir \${LTR_PEEL_TOOLS_DIR:-/opt/ltr_peel/tools} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ltr_peel: \$(python -c 'import ltr_peel; print(ltr_peel.__version__)')
        blastn: \$(blastn -version 2>&1 | head -1 | sed 's/^blastn: //')
    END_VERSIONS
    """

    stub:
    """
    cp in/*_clean_ltr.tsv in/*_clean_ltr.fa .
    for id in ${ids.join(' ')}; do
      printf '#old_seq_id\\tnew_seq_id\\tdecision\\n' > \${id}_reboundary.tsv
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ltr_peel: 1.0.1
    END_VERSIONS
    """
}
