#!/usr/bin/env bash

module load hmmer/3.3.2

ORF_DIR="output/orfs"
HMM_DIR="output/hmm/hmms"

MAP_MARKERS="input/mapfile.tsv"
MAP_GENOMES="output/config/genome_marker_map.tsv"

OUT_DIR="output/hmm/search"
POOL_DIR="output/hmm/pools"

# External RdRp-scan database
RDRPSCAN_DB="${HMM_DIR}/rdrpscan/RdRp_HMM_profile_CLUSTALO.db"

mkdir -p "$OUT_DIR" "$POOL_DIR"


# For each marker_group_id defined in input/mapfile.tsv
tail -n +2 "$MAP_MARKERS" |
  cut -f5 |
  sort -u |
  while read -r marker; do

    [[ -z "$marker" ]] && continue
    [[ "$marker" == "NA" ]] && continue

    echo "[MARKER] $marker"

    pool="$POOL_DIR/${marker}.faa"
    tbl="$OUT_DIR/${marker}.tbl"
    log="$OUT_DIR/${marker}.log"

    # ----------------------------------------------------------
    # Build the protein pool for genomes associated with marker
    # ----------------------------------------------------------

    > "$pool"

    awk -F'\t' -v marker="$marker" '
      NR == 1 {
        for (i = 1; i <= NF; i++) {
          if ($i == "virus_id" || $i == "genome_id") genome_col = i
          if ($i == "marker_group_id") marker_col = i
        }
        next
      }

      $marker_col == marker {
        print $genome_col
      }
    ' "$MAP_GENOMES" |
      sort -u |
      while read -r genome; do

        faa="$ORF_DIR/${genome}.faa"

        if [[ -s "$faa" ]]; then
          cat "$faa" >> "$pool"
        else
          echo "[WARN] missing ORF file: $faa"
        fi

      done


    # ----------------------------------------------------------
    # Search marker
    # ----------------------------------------------------------

    if [[ "$marker" == "rdrp_orthornavirae" ]]; then

      echo "[SEARCH] hmmscan RdRp-scan -> $marker"

      if [[ ! -s "${RDRPSCAN_DB}.h3m" ]]; then
        echo "[ERROR] RdRp-scan database not found:"
        echo "        $RDRPSCAN_DB"
        exit 1
      fi

      hmmscan \
        --cpu "${SLURM_CPUS_PER_TASK:-1}" \
        --noali \
        -E 1e-6 \
        --domtblout "$tbl" \
        "$RDRPSCAN_DB" \
        "$pool" \
        > "$log"

    else

      hmm="$HMM_DIR/${marker}.hmm"

      echo "[SEARCH] hmmsearch -> $marker"

      if [[ ! -s "$hmm" ]]; then
        echo "[ERROR] HMM not found: $hmm"
        exit 1
      fi

      hmmsearch \
        --cpu "${SLURM_CPUS_PER_TASK:-1}" \
        --noali \
        -E 1e-5 \
        --domtblout "$tbl" \
        "$hmm" \
        "$pool" \
        > "$log"

    fi

    echo "[OK] $tbl"

done
