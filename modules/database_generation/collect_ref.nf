process COLLECT_REF {

    // Seuls les profils RdRp-scan sont publiés ; les FASTA restent dans work/.
    publishDir "${params.outdir_generation}", mode: 'copy', pattern: 'hmms/rdrpscan'

    input:
    path marker_map
    path manual_references
    path reference_data
    path script_r

    output:
    path 'references_protein/*.fasta', optional: true, emit: proteins
    path 'reference_sources', emit: sources
    path 'hmms/rdrpscan', optional: true, emit: rdrpscan

    script:
    """
    Rscript "${script_r}" "${marker_map}" "${manual_references}" "${reference_data}"
    """
}
