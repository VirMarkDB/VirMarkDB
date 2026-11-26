process FILTER_POOL {
    tag "${marker_id}"
    publishDir "${params.outdir}/assessment/filtered_pool", mode: 'copy', pattern: 'filtered/*'
    input:
    tuple val(marker_id), path(pool), path(alignment)
    path script_r
    output:
    tuple val(marker_id), path("filtered/${marker_id}.fasta"), emit: proteins
    path 'filtered/*_summary.tsv', emit: summary
    path 'filtered/*.pdf', optional: true, emit: plots
    script:
    """
    Rscript "${script_r}" filter "${marker_id}" "${pool}" "${alignment}"
    """
}
