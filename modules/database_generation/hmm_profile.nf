process HMM_PROFILE {

    tag "${ref_id}"
    publishDir "${params.outdir_generation}/hmms", mode: 'copy'

    input:
    tuple val(ref_id), path(proteins)

    output:
    tuple val(ref_id), path("${ref_id}.hmm"), emit: profiles

    script:
    """
    # Aligner les protéines de référence, puis construire leur profil HMM.
    mafft --amino --thread ${task.cpus} --ep 0 --genafpair --maxiterate 1000 "${proteins}" > "${ref_id}.aln"
    hmmbuild --amino --cpu ${task.cpus} -n "${ref_id}" "${ref_id}.hmm" "${ref_id}.aln"
    """
}
