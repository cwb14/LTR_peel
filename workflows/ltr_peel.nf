/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LTR_peel: iterative nested LTR-RT detection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { LTR_PEEL_DETECT_ROUNDS } from '../subworkflows/local/ltr_peel_detect_rounds/main'
include { LTR_PEEL_RECONCILE     } from '../modules/local/ltr_peel/reconcile/main'
include { LTR_PEEL_CLUSTER       } from '../modules/local/ltr_peel/cluster/main'
include { LTR_PEEL_FLAGFP        } from '../modules/local/ltr_peel/flagfp/main'
include { LTR_PEEL_RECOVERSTRAND_ALIGN } from '../modules/local/ltr_peel/recoverstrand/main'
include { LTR_PEEL_RECOVERSTRAND_APPLY } from '../modules/local/ltr_peel/recoverstrand/main'
include { LTR_PEEL_ANNOTATE      } from '../modules/local/ltr_peel/annotate/main'
include { LTR_PEEL_GFF3          } from '../modules/local/ltr_peel/gff3/main'
include { LTR_PEEL_REBOUNDARY   } from '../modules/local/ltr_peel/reboundary/main'
include { LTR_PEEL_PLOTS         } from '../modules/local/ltr_peel/plots/main'

workflow LTR_PEEL {

    take:
    ch_samples  // channel: [ meta, genome, proteins ]

    main:
    // nextflow_schema.json is documentation, not a gate: nothing calls
    // validateParameters(). Check here so a typo'd preset fails now rather than
    // at the recovery stage, which sits after every expensive one. The CLI does
    // the same at ltr_peel.sh's argument validation.
    // null, false and 'off' all mean no strand recovery. The PPT fallback is on
    // by default, so with recovery off it is simply unused rather than an error.
    // Switches are read through toString().toBoolean(): Nextflow 24.x turns
    // `--reboundary false` into the Boolean false, but newer releases pass the String
    // "false" (seen on 26.10), which Groovy counts as true, so the off switch
    // would be ignored.
    def strand_on = params.strand_recovery
        && !(params.strand_recovery.toString() in ['off', 'false'])
    def reboundary_on = params.reboundary.toString().toBoolean()
    if (strand_on
            && !(params.strand_recovery in ['conservative', 'balanced', 'sensitive'])) {
        error "--strand_recovery must be conservative, balanced, sensitive or off " +
              "(got '${params.strand_recovery}')"
    }
    if (!(params.fastltr in ['all', 'rescue', 'off'])) {
        error "--fastltr must be all, rescue or off (got '${params.fastltr}')"
    }
    // Matched on the text: a value given on the command line may arrive as
    // a Long or a String rather than an Integer, and must still be accepted.
    if (!(params.ltr_timeout.toString() ==~ /\d+/)) {
        error "--ltr_timeout must be a whole number of seconds, 0 for no cap " +
              "(got '${params.ltr_timeout}')"
    }

    ch_versions = Channel.empty()

    //
    // Phase 1 - per genome: detect, mask, detect again, until the well runs dry.
    //
    LTR_PEEL_DETECT_ROUNDS(ch_samples)
    ch_versions = ch_versions.mix(LTR_PEEL_DETECT_ROUNDS.out.versions)

    //
    // Phase 2 - per genome: pool the rounds and resolve containment.
    //
    // Round N cannot know whether what it just found sits inside something round
    // N-1 found, so containment is settled once, here, over the union. The two
    // lists are sorted by round so their indices line up.
    //
    ch_reconcile_in = LTR_PEEL_DETECT_ROUNDS.out.results
        .map { meta, round, tsv, fasta, workdir -> [ meta, [ round, tsv, fasta, workdir ] ] }
        .groupTuple()
        .map { meta, per_round ->
            def ordered = per_round.sort { a, b -> a[0] <=> b[0] }
            [ meta.id,
              meta,
              ordered.collect { it[1] },
              ordered.collect { it[2] } ]
        }
        // The genome too: a host a tenant in one of its LTRs split into two
        // calls is rebuilt from the bases beyond that tenant.
        .join(ch_samples.map { meta, genome, _proteins -> [ meta.id, genome ] })
        .map { _id, meta, tsvs, fastas, genome -> [ meta, tsvs, fastas, genome ] }

    LTR_PEEL_RECONCILE(ch_reconcile_in)
    ch_versions = ch_versions.mix(LTR_PEEL_RECONCILE.out.versions)

    //
    // Phase 3 - pooled across every genome: one clustering pass, one
    // false-positive call.
    //
    // Pooling is the point. Families computed per species are not comparable
    // between species, and a repeat that looks convincing in one genome but
    // wrong across several only gets caught when they are judged together.
    //
    def pool = [ id: params.family_prefix ]

    ch_pooled_fasta = LTR_PEEL_RECONCILE.out.fasta.map { _meta, files -> files }.collect()
    ch_pooled_tsv   = LTR_PEEL_RECONCILE.out.tsv.map   { _meta, files -> files }.collect()

    LTR_PEEL_CLUSTER(ch_pooled_fasta.map { files -> [ pool, files ] })
    ch_versions = ch_versions.mix(LTR_PEEL_CLUSTER.out.versions)

    // ltr_peel.flag_fp needs a genome to write its (here unused) masked FASTA.
    // The CLI would re-run the whole pipeline on that mask when false positives
    // are pervasive; see docs/nextflow.md for why this pipeline reports instead.
    ch_first_genome = ch_samples.map { _meta, genome, _proteins -> genome }.first()

    // The tandem-array purge judges each call against its own genome's flanks,
    // so it needs every sample's genome, not just the first. Sorted by name so
    // the task hash, and with it -resume, does not depend on arrival order.
    // Off, an empty list flows in its place and stages nothing.
    ch_tandem_genomes = params.skip_tandem_filter.toString().toBoolean()
        ? Channel.value([])
        : ch_samples.map { _meta, genome, _proteins -> genome }
            .toSortedList { a, b -> a.name <=> b.name }

    LTR_PEEL_FLAGFP(
        LTR_PEEL_CLUSTER.out.consensus_cluster
            .join(LTR_PEEL_CLUSTER.out.internal_cluster)
            .join(LTR_PEEL_CLUSTER.out.consensus_fasta)
            .combine(ch_pooled_tsv.map { files -> [files] })
            .combine(ch_pooled_fasta.map { files -> [files] })
            .combine(ch_first_genome)
            .combine(ch_tandem_genomes.map { files -> [files] })
    )
    ch_versions = ch_versions.mix(LTR_PEEL_FLAGFP.out.versions)

    //
    // Phase 4 - back to per genome, judged against the one pooled cluster table,
    // so `family` denotes the same family in every genome's output.
    //
    // The FP stage returns one flat pile of cleaned files for every genome at
    // once; they are re-attributed by the prefix in their own names, which is
    // the same prefix the per-genome tables were written with.
    //
    ch_clean_by_id = LTR_PEEL_FLAGFP.out.tsv.map { _meta, files -> files }.flatten()
        .mix(LTR_PEEL_FLAGFP.out.fasta.map { _meta, files -> files }.flatten())
        .map { f -> [ f.name.replaceFirst(/_depth\d+_clean_ltr\.(tsv|fa)$/, ''), f ] }
        .groupTuple()

    // The rejoin map rides with the raw tables, so every stage that stages them
    // (annotate, recover-strand) finds it beside them by name.
    ch_raw_by_id = LTR_PEEL_RECONCILE.out.tsv.map { meta, files -> [ meta.id, files ] }
        .join(LTR_PEEL_RECONCILE.out.fasta.map { meta, files -> [ meta.id, files ] })
        .join(LTR_PEEL_RECONCILE.out.rejoined.map { meta, file -> [ meta.id, file ] })
        .map { id, tsvs, fastas, rejoined -> [ id, tsvs + fastas + [ rejoined ] ] }
    ch_rejoined_by_id = LTR_PEEL_RECONCILE.out.rejoined.map { meta, file -> [ meta.id, file ] }

    ch_work_by_id = LTR_PEEL_DETECT_ROUNDS.out.results
        .map { meta, _round, _tsv, _fasta, workdir -> [ meta.id, workdir ] }
        .groupTuple()

    // One bundle per genome: everything the annotation, GFF3 and plotting stages
    // read, keyed by sample so each stage can take just the slice it needs.
    ch_bundle = ch_samples
        .map { meta, genome, _proteins -> [ meta.id, meta, genome ] }
        .join(ch_raw_by_id)
        .join(ch_clean_by_id)
        .join(ch_work_by_id)
        .combine(LTR_PEEL_CLUSTER.out.consensus_cluster.map { _meta, tsv -> tsv })
        .map { _id, meta, genome, raw, clean, workdirs, cluster ->
            def files = raw + clean
            [ meta,
              files.findAll { it.name.endsWith('.tsv') },
              files.findAll { it.name.endsWith('.fa')  },
              workdirs, genome, cluster ]
        }

    // Strand recovery, first half. On unless --strand_recovery off. It calls strand for the
    // elements the annotator's cascade cannot reach and writes the sidecar the
    // annotator reads back as a fourth tier. With recovery off, an empty list
    // flows in its place: it stages nothing and the modules render no flag for
    // it, so the run is exactly what it was before.
    if (strand_on) {
        LTR_PEEL_RECOVERSTRAND_ALIGN(
            ch_bundle.map { meta, tsvs, _fastas, workdirs, genome, _cluster ->
                [ meta, tsvs, workdirs, genome ]
            }
        )
        ch_versions = ch_versions.mix(LTR_PEEL_RECOVERSTRAND_ALIGN.out.versions)
        ch_recovery = LTR_PEEL_RECOVERSTRAND_ALIGN.out.recovery
    } else {
        ch_recovery = ch_bundle.map { meta, _t, _f, _w, _g, _c -> [ meta, [] ] }
    }

    LTR_PEEL_ANNOTATE(
        ch_bundle.map { meta, tsvs, _fastas, workdirs, _genome, cluster ->
            [ meta, tsvs, workdirs, cluster ]
        }.join(ch_recovery)
    )
    ch_versions = ch_versions.mix(LTR_PEEL_ANNOTATE.out.versions)

    // Second half: bring the _clean_ depth FASTAs into line with the strand
    // column the annotator just wrote, so a recovered minus element is stored in
    // coding sense exactly like a TEsorter2-called one. The unfiltered depth
    // FASTAs are left as detected. The GFF3 carries no sequence,
    // so it is indifferent to this and could equally have run first.
    if (strand_on) {
        LTR_PEEL_RECOVERSTRAND_APPLY(
            LTR_PEEL_ANNOTATE.out.tsv
                .join(ch_bundle.map { meta, _t, fastas, _w, _g, _c -> [ meta, fastas ] })
                .join(ch_recovery)
        )
        ch_versions = ch_versions.mix(LTR_PEEL_RECOVERSTRAND_APPLY.out.versions)
        ch_tables = LTR_PEEL_RECOVERSTRAND_APPLY.out.tsv
        ch_fastas = LTR_PEEL_RECOVERSTRAND_APPLY.out.fasta
    } else {
        ch_tables = LTR_PEEL_ANNOTATE.out.tsv
        ch_fastas = ch_bundle.map { meta, _t, fastas, _w, _g, _c -> [ meta, fastas ] }
    }

    // Re-boundarying, when asked for: one pooled task over every sample, between
    // the annotated tables and the GFF3. Off, an empty map flows in its place.
    ch_meta_by_id = ch_bundle.map { meta, _t, _f, _w, _g, _c -> [ meta.id, meta ] }
    if (reboundary_on) {
        // Kept so a sample the pooled process returns nothing for can fall
        // through unchanged below: every one of its _clean_ tables was purged
        // as a false positive, so it had nothing to re-boundary. This mirrors
        // ltr_peel.sh's second GFF3 loop, which writes that genome's GFF3
        // straight from these same pre-reboundary annotated tables when
        // ${p}_reboundary.tsv does not exist for it.
        ch_tables_prerb = ch_tables
        ch_fastas_prerb = ch_fastas

        ch_rb_in = ch_tables
            .join(ch_fastas)
            .join(ch_bundle.map { meta, _t, _f, _w, genome, _c -> [ meta, genome ] })
            .toSortedList { a, b -> a[0].id <=> b[0].id }
        LTR_PEEL_REBOUNDARY(
            ch_rb_in.map { rows -> rows.collect { it[0].id } },
            ch_rb_in.map { rows -> rows.collect { it[1] }.flatten() },
            ch_rb_in.map { rows -> rows.collect { it[2] }.flatten() },
            ch_rb_in.map { rows -> rows.collect { it[3] } }
        )
        ch_versions = ch_versions.mix(LTR_PEEL_REBOUNDARY.out.versions)

        // Ids REBOUNDARY actually wrote files for: those with a surviving
        // _clean_ table to begin with. Any other id that entered this stage is
        // simply absent from these outputs -- there was nothing to pool for it
        // -- and is filled back in below from ch_*_prerb, unchanged.
        ch_rb_tables = LTR_PEEL_REBOUNDARY.out.tsv.flatten()
            .map { f -> [ f.name.replaceFirst(/_depth\d+_clean_ltr\.tsv$/, ''), f ] }
            .groupTuple()
        ch_rb_fastas = LTR_PEEL_REBOUNDARY.out.fasta.flatten()
            .map { f -> [ f.name.replaceFirst(/_depth\d+_clean_ltr\.fa$/, ''), f ] }
            .groupTuple()
        ch_rb_sidecar = LTR_PEEL_REBOUNDARY.out.sidecar.flatten()
            .map { f -> [ f.name.replaceFirst(/_reboundary\.tsv$/, ''), f ] }

        // remainder: true so an id with no reboundary output still passes
        // through the join instead of being dropped; ch_meta_by_id covers
        // every id that entered this stage, so it is always the left side.
        ch_tables = ch_meta_by_id.join(ch_rb_tables, remainder: true)
            .join(ch_tables_prerb.map { meta, files -> [ meta.id, files ] })
            .map { _id, meta, rb, orig -> [ meta, rb ?: orig ] }
        ch_fastas = ch_meta_by_id.join(ch_rb_fastas, remainder: true)
            .join(ch_fastas_prerb.map { meta, files -> [ meta.id, files ] })
            .map { _id, meta, rb, orig -> [ meta, rb ?: orig ] }
        ch_rbmap = ch_meta_by_id.join(ch_rb_sidecar, remainder: true)
            .map { _id, meta, f -> [ meta, f ?: [] ] }
    } else {
        ch_rbmap = ch_bundle.map { meta, _t, _f, _w, _g, _c -> [ meta, [] ] }
    }

    // From here on the tables are the ANNOTATED ones, not the reconciler's.
    ch_annotated = ch_tables
        .join(ch_fastas)
        .join(ch_bundle.map { meta, _tsvs, _fastas, workdirs, genome, cluster ->
            [ meta, workdirs, genome, cluster ]
        })

    // The GFF3 stage reads the annotated tables, which do not carry the rejoin
    // map, so it is added back: strand_source must name the same tier the
    // annotator used.
    LTR_PEEL_GFF3(
        ch_annotated.map { meta, tables, _fastas, workdirs, genome, cluster ->
            [ meta.id, meta, tables, workdirs, genome, cluster ]
        }.join(ch_rejoined_by_id).map { _id, meta, tables, workdirs, genome, cluster, rejoined ->
            [ meta, (tables instanceof List ? tables : [ tables ]) + [ rejoined ],
              workdirs, genome, cluster ]
        }.join(ch_recovery).join(ch_rbmap)
    )
    ch_versions = ch_versions.mix(LTR_PEEL_GFF3.out.versions)

    //
    // Phase 5 - plots.
    //
    if (!params.skip_plots.toString().toBoolean()) {
        LTR_PEEL_PLOTS(
            ch_annotated.map { meta, tables, fastas, workdirs, genome, cluster ->
                [ meta, tables, fastas, workdirs, genome, cluster ]
            }
        )
        ch_versions = ch_versions.mix(LTR_PEEL_PLOTS.out.versions)
    }

    emit:
    depth_tables = ch_tables
    gff3         = LTR_PEEL_GFF3.out.gff3
    versions     = ch_versions
}
