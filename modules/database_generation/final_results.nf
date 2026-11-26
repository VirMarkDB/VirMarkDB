process FINAL_RESULTS {

    publishDir "${params.outdir_generation}", mode: 'copy'

    input:
    path manifest, stageAs: 'manifest_genomes.tsv'
    path marker_map, stageAs: 'marker_taxo_map.tsv'
    path genome_markers, stageAs: 'genome_marker_map.tsv'
    path hits, stageAs: 'hits/*'
    path proteins, stageAs: 'proteins/*'
    path cds, stageAs: 'cds/*'
    path remove_sequences, stageAs: 'remove_specific_ORF_sequences.txt'
    path script_r

    output:
    path 'VirMarkDB', emit: database

    script:
    """
    Rscript "${script_r}"
    """
}
