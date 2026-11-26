include { COLLECT_POOL } from '../modules/assessment/collect_pool'
include { ALIGN_POOL as SELF_ALIGN; ALIGN_POOL as DATABASE_ALIGN } from '../modules/assessment/align_pool'
include { FILTER_POOL } from '../modules/assessment/filter_pool'
include { ANALYSE_POOL } from '../modules/assessment/analyse_pool'
include { CLUSTER_MISSING } from '../modules/assessment/cluster_missing'

workflow ASSESSMENT {
    take:
    database
    marker_map

    main:
    analysis_script = file("${projectDir}/scripts/assessment/analyse_pool.R", checkIfExists: true)
    COLLECT_POOL(marker_map, file(params.assessment_rules, checkIfExists: true), params.nr_database,
                 file("${projectDir}/scripts/assessment/collect_pool.bash", checkIfExists: true))

    pools = COLLECT_POOL.out.pool
        .map { dir -> file("${dir}/*.fasta") }
        .flatten()
        .map { fasta -> tuple(fasta.baseName, fasta, fasta, 'self') }
    SELF_ALIGN(pools)
    FILTER_POOL(SELF_ALIGN.out.alignments, analysis_script)

    references = database
        .map { dir -> file("${dir}/markers/*/*/*_protein.fasta") }
        .flatten()
        .map { fasta -> tuple(fasta.name.replaceFirst(/_protein\.fasta$/, ''), fasta) }
    DATABASE_ALIGN(FILTER_POOL.out.proteins.join(references)
        .map { id, pool, ref -> tuple(id, pool, ref, 'database') })
    ANALYSE_POOL(DATABASE_ALIGN.out.alignments, database, analysis_script)
    CLUSTER_MISSING(ANALYSE_POOL.out.missing)

    emit:
    results = ANALYSE_POOL.out.results
    clusters = CLUSTER_MISSING.out.clusters
}
