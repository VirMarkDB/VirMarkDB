process HMM_SEARCH {

    tag "${ref_id} (${tool})"
    publishDir "${params.outdir_generation}/search", mode: 'copy'

    input:
    tuple val(ref_id), path(proteins), val(tool), path(profile)

    output:
    tuple val(ref_id), path("${ref_id}.tbl"), emit: hits
    tuple val(ref_id), path("${ref_id}.log"), emit: logs

    script:
    // RdRp-scan uses a pressed database; other references use individual HMM files.
    def hmm_file = tool == 'hmmscan' ? "${profile}/RdRp_HMM_profile_CLUSTALO.db" : profile
    """
    cat *.faa > pool.fasta
    ${tool} --cpu ${task.cpus} --noali -E ${params.hmm_evalue} --domtblout "${ref_id}.tbl" "${hmm_file}" pool.fasta > "${ref_id}.log"
    """
}
