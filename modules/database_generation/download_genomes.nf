process DOWNLOAD_GENOMES {

    tag "batch_${batch_id}"
    maxForks params.genome_download_parallelism

    publishDir "${params.outdir_generation}/viruses", mode: 'copy'

    input:
    tuple val(batch_id), val(genome_ids)
    val references_ready

    output:
    tuple val(genome_ids), path('*/*.fna', arity: '1..*'), emit: genomes

    script:
    """
    for virus_id in ${genome_ids.join(' ')}; do
        mkdir -p "\$virus_id"
        ncbi_id=\${virus_id%_partial}
        efetch -db nuccore -id "\$ncbi_id" -format fasta > "\$virus_id/\$virus_id.fna"
        test -s "\$virus_id/\$virus_id.fna"
        sleep ${params.genome_download_parallelism.toInteger() * 0.4}
    done
    """
}
