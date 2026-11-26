process ALIGN_POOL {
    tag "${marker_id} (${mode})"
    label 'process_medium'
    publishDir { "${params.outdir}/assessment/alignments/${mode}" }, mode: 'copy', pattern: '*.m8'
    input:
    tuple val(marker_id), path(query), path(target), val(mode)
    output:
    tuple val(marker_id), path(query), path("${marker_id}.m8"), emit: alignments
    script:
    """
    if [ -s "${query}" ] && [ -s "${target}" ]; then
        mmseqs easy-search "${query}" "${target}" "${marker_id}.m8" tmp \
            --threads ${task.cpus} --format-output 'query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits'
    else
        touch "${marker_id}.m8"
    fi
    """
}
