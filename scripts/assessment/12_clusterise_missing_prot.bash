#!/usr/bin/env bash

module load cd-hit/4.8.1

in_dir="output/assessment/alignment_VirMarkDB_pool/missing_prot"
out_dir="output/assessment/alignment_VirMarkDB_pool/missing_prot/clustering"

for fasta in "${in_dir}"/*.fasta; do

    marker=$(basename "${fasta}" .fasta)

    mkdir -p "${out_dir}/${marker}"

    echo "Input : ${fasta}"
    echo "Marker: ${marker}"

    cd-hit \
        -i "${fasta}" \
        -o "${out_dir}/${marker}/${marker}_cdhit95.fasta" \
        -c 0.95 \
        -n 5 \
        -d 0 \
        -M 0 \
        > "${out_dir}/${marker}/${marker}_cdhit95.log" 2>&1

done

