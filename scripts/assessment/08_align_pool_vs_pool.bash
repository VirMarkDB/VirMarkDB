#!/usr/bin/env bash
#SBATCH --job-name=align_VirMarkDB_VirMarkDB
#SBATCH --partition=fast
#SBATCH --cpus-per-task=50
#SBATCH --output=align_VirMarkDB_VirMarkDB.out

module load mmseqs2/15.6f452

mkdir -p output/assessment/pool_protein/auto_alignment/tmp/tmp_aln
mkdir -p output/assessment/pool_protein/auto_alignment/blast

pool="output/assessment/pool_protein/pool_protein_raw/"

for marker in "$pool"*.fasta; do
    name=$(basename "$marker" .fasta)

    echo "Running: $name"

    mmseqs easy-search \
        "$marker" \
        "$marker" \
        "output/assessment/pool_protein/auto_alignment/blast/${name}.m8" \
        output/assessment/pool_protein/auto_alignment/tmp/tmp_aln \
        --threads 50 
done