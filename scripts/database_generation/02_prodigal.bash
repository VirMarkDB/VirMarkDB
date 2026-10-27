#!/usr/bin/env bash

module load prodigal/2.6.3

GENOME_DIR="output/genomes"
OUT="output/orfs"

mkdir -p "$OUT"

for genome in "$GENOME_DIR"/*.fna; do
  id="$(basename "$genome" .fna)"

  faa="$OUT/${id}.faa"
  fna="$OUT/${id}.fna"
  gff="$OUT/${id}.gff"

  if [[ -s "$faa" && -s "$fna" && -s "$gff" ]]; then
    echo "[SKIP] $id already processed"
    continue
  fi

  prodigal -i "$genome" -a "$OUT/${id}.faa" -d "$OUT/${id}.fna" -f gff -o "$OUT/${id}.gff" -p meta

done
