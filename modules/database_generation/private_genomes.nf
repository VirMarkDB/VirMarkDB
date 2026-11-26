process PRIVATE_GENOMES {

    tag "${virus_id}"

    publishDir "${params.outdir_generation}/viruses", mode: 'copy'

    input:
    tuple val(virus_id), path(genome)

    output:
    tuple val(virus_id), path("${virus_id}/${virus_id}.fna"), emit: genomes

    script:
    """
    mkdir -p "${virus_id}"
    cp "${genome}" "${virus_id}/${virus_id}.fna"
    """

}
