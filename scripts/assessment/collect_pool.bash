#!/usr/bin/env bash
set -euo pipefail

OUTDIR="pool"
METADIR="${OUTDIR}/metadata"
TAXDB_DIR="taxonomy"
if [[ $# -ne 3 ]]; then
  echo "Usage: collect_pool.bash <nr_database> <marker_map.tsv> <assessment_rules.tsv>" >&2
  exit 1
fi
DB="$1"
MARKER_MAP="$2"
RULES="$3"
SHARED_NR_DIR="$(dirname "$DB")"

# Validate and select rules before accessing NR or downloading taxonomy data.
# Join by marker_group_id; unused rules are allowed for subset mapfiles.
awk -F '\t' -v OFS='\t' '
  FNR==1 {
    delete cols
    for(i=1;i<=NF;i++) cols[$i]=i
    if (!("marker_group_id" in cols) || !("group_id" in cols)) {
      print "Missing marker_group_id or group_id column in " FILENAME > "/dev/stderr"; exit 1
    }
    if (NR==FNR) {mi=cols["marker_group_id"]; gi=cols["group_id"]}
    else {
      split("marker_group_id group_id taxid include_pattern exclude_pattern min_length max_length", required, " ")
      for(i=1;i<=7;i++) if (!(required[i] in cols)) {
        print "Missing rule column: " required[i] > "/dev/stderr"; exit 1
      }
    }
    next
  }
  NR==FNR {
    if (($mi in wanted) && wanted[$mi]!=$gi) {
      print "Conflicting marker groups: " $mi > "/dev/stderr"; exit 1
    }
    wanted[$mi]=$gi; next
  }
  {
    id=$(cols["marker_group_id"])
    if (!(id in wanted)) next
    group=$(cols["group_id"]); taxid=$(cols["taxid"])
    low=$(cols["min_length"]); high=$(cols["max_length"])
    pattern=$(cols["include_pattern"]); exclude=$(cols["exclude_pattern"])
    if (seen[id]++ || group!=wanted[id] || id!~/^[A-Za-z0-9_.-]+$/ ||
        group!~/^[A-Za-z0-9_.-]+$/ || taxid!~/^[0-9]+$/ || taxid+0<1 ||
        low!~/^[0-9]+$/ || high!~/^[0-9]+$/ || low+0>high+0 || pattern=="" || exclude=="") {
      print "Invalid or duplicate assessment rule: " id > "/dev/stderr"; failed=1; next
    }
    print id,group,taxid,pattern,exclude,low,high
  }
  END {
    for(id in wanted) if (!(id in seen)) {
      print "Missing assessment rule: " id > "/dev/stderr"; failed=1
    }
    if (failed) exit 1
  }
' "$MARKER_MAP" "$RULES" > selected_rules.tsv

mkdir -p "$OUTDIR" "$METADIR" "$TAXDB_DIR"

if [[ ! -s "$TAXDB_DIR/taxonomy4blast.sqlite3" ]]; then
  echo "Downloading current NCBI BLAST taxonomy database..."
  wget -q -O "$TAXDB_DIR/taxdb.tar.gz" https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
  tar -xzf "$TAXDB_DIR/taxdb.tar.gz" -C "$TAXDB_DIR"
  rm "$TAXDB_DIR/taxdb.tar.gz"
fi

export BLASTDB="${TAXDB_DIR}:${SHARED_NR_DIR}"

dump_taxon_titles() {
  local taxon="$1"
  local taxid="$2"
  local outfile="${METADIR}/${taxon}_${taxid}_titles.tsv"

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
  local taxid="$3"
  local pattern="$4"
  local exclude="$5"
  local min_len="$6"
  local max_len="$7"

  local titles="${METADIR}/${taxon}_${taxid}_titles.tsv"
  local ids="${METADIR}/${marker}.accessions.txt"
  local audit="${METADIR}/${marker}.selected_titles.tsv"
  local outfile="${OUTDIR}/${marker}.fasta"

  dump_taxon_titles "$taxon" "$taxid"

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

  echo -e "${marker}\t${taxon}\t${taxid}\t${n}" \
    >> "${METADIR}/marker_counts.tsv"
}

rm -f "${METADIR}/marker_counts.tsv"
echo -e "marker\ttaxon\ttaxid\tn_sequences" > "${METADIR}/marker_counts.tsv"

# The mapfile controls marker selection; the rule table supplies extraction criteria.
cp selected_rules.tsv "${METADIR}/selected_rules.tsv"
while IFS=$'\t' read -r marker group taxid pattern exclude min_len max_len; do
  fetch_marker "$marker" "$group" "$taxid" "$pattern" "$exclude" "$min_len" "$max_len"
done < selected_rules.tsv

echo
echo "Final marker counts:"
cat "${METADIR}/marker_counts.tsv"

