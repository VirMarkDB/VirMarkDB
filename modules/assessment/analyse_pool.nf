process ANALYSE_POOL {
    tag "${marker_id}"
    publishDir "${params.outdir}/assessment/results", mode: 'copy'
    input:
    tuple val(marker_id), path(pool), path(alignment)
    path database
    path script_r
    output:
    path 'results', emit: results
    tuple val(marker_id), path("missing/${marker_id}.fasta"), emit: missing
    path 'hits', emit: hits
    script:
    """
    Rscript "${script_r}" analyse "${marker_id}" "${pool}" "${alignment}" "${database}"
    """
}
