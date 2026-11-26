process COLLECT_POOL {
    publishDir "${params.outdir}/assessment/pool_protein", mode: 'copy'
    input:
    path marker_map
    path rules
    val nr_database
    path script_bash
    output:
    path 'pool', emit: pool
    script:
    """
    bash "${script_bash}" "${nr_database}" "${marker_map}" "${rules}"
    """
}
