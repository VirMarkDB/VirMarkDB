include { SUMMARY_STATISTICS } from '../modules/statistics/summary_statistics'

workflow STATISTICS {
    take:
    database
    manifest
    marker_map
    genome_markers

    main:
    SUMMARY_STATISTICS(database, manifest, marker_map, genome_markers,
                       file("${projectDir}/scripts/statistics/summary_statistics.R", checkIfExists: true))

    emit:
    statistics = SUMMARY_STATISTICS.out.statistics
}
