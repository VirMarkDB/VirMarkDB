process PRODIGAL {

    tag "batch_${batch_id}"

    publishDir "${params.outdir_generation}/viruses", mode: 'copy'

    input:
    tuple val(batch_id), val(genome_ids), path(genomes)

    output:
    tuple val(batch_id), val(genome_ids), path('*/*.faa', arity: '1..*'), emit: proteins
    tuple val(batch_id), val(genome_ids), path('*/*.ffn'), emit: cds
    tuple val(batch_id), val(genome_ids), path('*/*.gff'), emit: annotations

    script:
    """
    for virus_id in ${genome_ids.join(' ')}; do
        mkdir -p "\$virus_id"
        prodigal -i "\$virus_id.fna" -a "\$virus_id/\$virus_id.faa" -d "\$virus_id/\$virus_id.ffn" -f gff -o "\$virus_id/\$virus_id.gff" -p meta
    done
    """
}
