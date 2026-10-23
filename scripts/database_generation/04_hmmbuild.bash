#!/usr/bin/env bash

module load hmmer/3.3.2

IN="output/hmm/aln"
OUT="output/hmm/hmms"
mkdir -p "$OUT"

for aln in "$IN"/*.aln; do
  base="$(basename "$aln" .aln)"
  hmmbuild "$OUT/${base}.hmm" "$aln"
done
