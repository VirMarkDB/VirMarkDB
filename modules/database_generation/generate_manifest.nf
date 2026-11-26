#!/usr/bin/env nextflow


process GENERATE_MANIFEST {
    publishDir "${params.outdir}/config", mode: 'copy'

    input:
    path inputs, stageAs: 'input'
    path script_r
    val ictv_version

    output:
    path 'manifest_genomes.tsv', emit: manifest
    path 'marker_taxo_map.tsv', emit: marker_map
    path 'genome_marker_map.tsv', emit: genome_markers
    path 'dup_ids_check.tsv', optional: true, emit: duplicates

    script:
    """
    Rscript "${script_r}" "${ictv_version}"
    """
}
