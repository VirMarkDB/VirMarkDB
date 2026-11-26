#!/usr/bin/env nextflow

include { DATABASE_GENERATION } from './subworkflows/database_generation'
include { ASSESSMENT } from './subworkflows/assessment'
include { STATISTICS } from './subworkflows/statistics'

workflow {
    DATABASE_GENERATION()

    if (params.run_assessment) {
        if (!params.nr_database) error 'Assessment requires --nr_database /path/to/blast/nr'
        ASSESSMENT(DATABASE_GENERATION.out.database, DATABASE_GENERATION.out.marker_map)
    }

    if (params.run_statistics) {
        STATISTICS(DATABASE_GENERATION.out.database, DATABASE_GENERATION.out.manifest,
                   DATABASE_GENERATION.out.marker_map, DATABASE_GENERATION.out.genome_markers)
    }
}
