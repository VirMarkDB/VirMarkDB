#!/usr/bin/env bash

module load mafft/7.525

IN="output/references_protein"
OUT1="output/hmm/aln"

mkdir -p "$OUT1"

for fasta in "$IN"/*.fasta; do
  stem="$(basename "$fasta" .fasta)"

  mafft --ep 0 --genafpair --maxiterate 1000 "$fasta" > "${OUT1}/${stem}.aln"

done

echo "OK"
echo "MAFFT alignments : $OUT1"
