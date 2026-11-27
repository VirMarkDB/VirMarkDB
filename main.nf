#!/usr/bin/env nextflow

include { DATABASE_GENERATION } from './subworkflows/database_generation'
include { ASSESSMENT } from './subworkflows/assessment'
include { STATISTICS } from './subworkflows/statistics'

workflow {
    main:
    DATABASE_GENERATION()

    if (params.run_assessment) {
        if (!params.nr_database) error 'Assessment requires --nr_database /path/to/blast/nr'
        ASSESSMENT(DATABASE_GENERATION.out.database, DATABASE_GENERATION.out.marker_map)
    }

    if (params.run_statistics) {
        STATISTICS(DATABASE_GENERATION.out.database, DATABASE_GENERATION.out.manifest,
                   DATABASE_GENERATION.out.marker_map, DATABASE_GENERATION.out.genome_markers)
    }

    // Render after Nextflow has finished writing its native DAG and reports.
    onComplete:
    if (params.metro_map && workflow.success) {
        def command = ['python3', "${projectDir}/scripts/utils/render_metro.py".toString(),
                       file(params.dag_file).toString()]
        Runtime.getRuntime().addShutdownHook(new Thread({ ->
            try {
                def renderer = new ProcessBuilder(command).redirectErrorStream(true).start()
                renderer.inputStream.eachLine { line -> System.err.println(line) }
                if (renderer.waitFor() != 0) {
                    System.err.println('Metro rendering failed; workflow results are unaffected.')
                }
            }
            catch (Exception error) {
                System.err.println("Metro rendering unavailable: ${error.message}")
            }
        }))
    }
}
