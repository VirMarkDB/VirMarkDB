#!/usr/bin/env bash
#SBATCH --job-name=make_bench_pool
#SBATCH --partition=fast
#SBATCH --cpus-per-task=50
#SBATCH --output=make_bench_pool.out

module load blast/2.16.0

OUTDIR="output/assessment/pool_protein/pool_protein_raw"
METADIR="${OUTDIR}/metadata"
TAXDB_DIR="output/assessment/ncbi_taxonomy_config"
SHARED_NR_DIR="/shared/bank/nr/nr_2025-07-20/blast"

mkdir -p "$OUTDIR" "$METADIR" "$TAXDB_DIR"

if [[ ! -s "$TAXDB_DIR/taxonomy4blast.sqlite3" ]]; then
  echo "Downloading current NCBI BLAST taxonomy database..."
  wget -q -O "$TAXDB_DIR/taxdb.tar.gz" https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
  tar -xzf "$TAXDB_DIR/taxdb.tar.gz" -C "$TAXDB_DIR"
  rm "$TAXDB_DIR/taxdb.tar.gz"
fi

export BLASTDB="${TAXDB_DIR}:${SHARED_NR_DIR}"

DB="nr"

declare -A TAXIDS=(
  [nucleocytoviricota]="2732007"
  [orthornavirae]="2732396"
  [pantevenvirales]="3420545"
)

dump_taxon_titles() {
  local taxon="$1"
  local taxid="${TAXIDS[$taxon]}"
  local outfile="${METADIR}/${taxon}_titles.tsv"

  if [[ -s "$outfile" ]]; then
    echo "Using cached title file: $outfile"
    return 0
  fi

  echo "Extracting titles for ${taxon} / taxid ${taxid}"

  blastdbcmd \
    -db "$DB" \
    -dbtype prot \
    -taxids "$taxid" \
    -target_only \
    -outfmt "%a---%T---%S---%l---%t" \
    > "$outfile"
}

fetch_marker() {
  local marker="$1"
  local taxon="$2"
  local pattern="$3"
  local min_len="$4"
  local max_len="$5"

  local titles="${METADIR}/${taxon}_titles.tsv"
  local ids="${METADIR}/${marker}.accessions.txt"
  local audit="${METADIR}/${marker}.selected_titles.tsv"
  local outfile="${OUTDIR}/${marker}.fasta"

  dump_taxon_titles "$taxon"

  # Conservative exclusion terms
  local exclude='partial|fragment|fragments|hypothetical|low quality|frameshift|pseudogene|(^|[^A-Za-z])MAG([^A-Za-z]|$)'

  echo "Filtering marker: $marker"

  awk -F '---' \
    -v IGNORECASE=1 \
    -v pat="$pattern" \
    -v excl="$exclude" \
    -v min_len="$min_len" \
    -v max_len="$max_len" \
    -v audit="$audit" '
      ($5 ~ pat) &&
      ($5 !~ excl) &&
      ($4 >= min_len) &&
      ($4 <= max_len) {
        print $0 > audit
        print $1
      }
    ' "$titles" \
    | sort -u \
    > "$ids"

  if [[ -s "$ids" ]]; then
    blastdbcmd \
      -db "$DB" \
      -dbtype prot \
      -entry_batch "$ids" \
      -out "$outfile"
  else
    : > "$outfile"
  fi

  local n
  n=$(grep -c '^>' "$outfile" || true)

  echo -e "${marker}\t${taxon}\t${TAXIDS[$taxon]}\t${n}" \
    >> "${METADIR}/marker_counts.tsv"
}

rm -f "${METADIR}/marker_counts.tsv"
echo -e "marker\ttaxon\ttaxid\tn_sequences" > "${METADIR}/marker_counts.tsv"

# Nucleocytoviricota markers
fetch_marker "majorcapsid_nucleocytoviricota" "nucleocytoviricota" \
  'major capsid protein|major capsid' \
  250 1000

fetch_marker "dnapol_nucleocytoviricota" "nucleocytoviricota" \
  'DNA polymerase B|B-family DNA polymerase|family B DNA polymerase|(^|[^A-Za-z])PolB([^A-Za-z]|$)|(^|[^A-Za-z])pPolB([^A-Za-z]|$)' \
  300 1600

fetch_marker "atpase_nucleocytoviricota" "nucleocytoviricota" \
  'packaging ATPase|DNA packaging ATPase|virion packaging ATPase|A32-like ATPase|A32 ATPase' \
  200 900

fetch_marker "primase_nucleocytoviricota" "nucleocytoviricota" \
  '(^|[^A-Za-z])primase([^A-Za-z]|$)|DNA primase|primase-helicase' \
  150 1400

fetch_marker "vltf3_nucleocytoviricota" "nucleocytoviricota" \
  'late transcription factor 3|(^|[^A-Za-z])VLTF3([^A-Za-z]|$)|viral late transcription factor 3' \
  100 700

fetch_marker "tf2s_nucleocytoviricota" "nucleocytoviricota" \
  'transcription elongation factor SII|(^|[^A-Za-z])TFIIS([^A-Za-z]|$)|transcription factor S-II|transcription factor SII' \
  100 700

fetch_marker "rnapol1_nucleocytoviricota" "nucleocytoviricota" \
  'DNA-directed RNA polymerase subunit 1|RNA polymerase subunit 1|RNA polymerase large subunit|RPO1|RPB1' \
  500 2500

fetch_marker "rnapol2_nucleocytoviricota" "nucleocytoviricota" \
  'DNA-directed RNA polymerase subunit 2|RNA polymerase subunit 2|RNA polymerase second largest subunit|RPO2|RPB2' \
  300 1800

# Orthornavirae marker
fetch_marker "rdrp_orthornavirae" "orthornavirae" \
  'RNA-dependent RNA polymerase|RNA-directed RNA polymerase|RNA-dependent RNA-directed RNA polymerase|RdRp|RDRP|large polymerase protein|polymerase protein L|RNA polymerase L' \
  250 3500

# Pantevenvirales marker
fetch_marker "g23_pantevenvirales" "pantevenvirales" \
  '(^|[^A-Za-z0-9])(gp23|g23)([^A-Za-z0-9]|$)|gene 23 protein|major capsid protein|capsid protein gp23|T4-like major capsid protein|T-even major capsid protein' \
  250 900

echo
echo "Final marker counts:"
column -t -s $'\t' "${METADIR}/marker_counts.tsv"

