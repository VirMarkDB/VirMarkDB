#!/usr/bin/env bash
#SBATCH --job-name=align_VirMarkDB_pool
#SBATCH --partition=fast
#SBATCH --cpus-per-task=50
#SBATCH --output=align_VirMarkDB_pool.out

module load mmseqs2/15.6f452

mkdir -p output/assessment/alignment_VirMarkDB_pool/alignment/tmp

VirMarkDB_db="output/VirMarkDB/markers"
pool_dir="output/assessment/pool_protein/pool_protein_filt"
out_dir="output/assessment/alignment_VirMarkDB_pool/alignment"

for ref in "$VirMarkDB_db"/*/*/*_protein.fasta; do
    full_name=$(basename "$ref" _protein.fasta)
    query="${pool_dir}/${full_name}.fasta"

    echo "Running marker_group_id: $full_name"
    echo "Using pool marker: $full_name"
    echo "$full_name"
    echo "$ref"

    mmseqs easy-search \
        "$query" \
        "$ref" \
        "${out_dir}/${full_name}.m8" \
        output/assessment/alignment_VirMarkDB_pool/alignment/tmp \
        --threads 50
done

echo "All alignments done"