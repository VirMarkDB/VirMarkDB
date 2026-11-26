process CLUSTER_MISSING {
    tag "${marker_id}"
    label 'process_medium'
    publishDir "${params.outdir}/assessment/clustering", mode: 'copy'
    input:
    tuple val(marker_id), path(proteins)
    output:
    path "${marker_id}", emit: clusters
    script:
    """
    mkdir -p "${marker_id}"
    if [ -s "${proteins}" ]; then
        cd-hit -i "${proteins}" -o "${marker_id}/${marker_id}_cdhit95.fasta" \
            -c 0.95 -n 5 -d 0 -T ${task.cpus} -M ${task.memory.toMega()} > "${marker_id}/clustering.log" 2>&1
    else
        touch "${marker_id}/${marker_id}_cdhit95.fasta" "${marker_id}/${marker_id}_cdhit95.fasta.clstr"
    fi
    """
}
