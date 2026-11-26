process SUMMARY_STATISTICS {
    publishDir "${params.outdir}/statistics", mode: 'copy'
    input:
    path database
    path manifest, stageAs: 'manifest_genomes.tsv'
    path marker_map, stageAs: 'marker_taxo_map.tsv'
    path genome_markers, stageAs: 'genome_marker_map.tsv'
    path script_r
    output:
    path 'analysis', emit: statistics
    script:
    """
    Rscript "${script_r}" "${database}"
    """
}
