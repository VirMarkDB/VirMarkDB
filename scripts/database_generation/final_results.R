library(tidyverse)

# Nextflow stages inputs in separate directories to avoid filename collisions.
HMM_DIR      <- "hits"
ORF_DIR      <- "proteins"
CDS_DIR      <- "cds"
MANIFEST_TSV <- "manifest_genomes.tsv"
MAP_TSV      <- "marker_taxo_map.tsv"
GENOME_MAP   <- "genome_marker_map.tsv"
OUT_DIR      <- "VirMarkDB"
OUT_LOGS     <- file.path(OUT_DIR, "logs")
REMOVE_SEQS  <- "remove_specific_ORF_sequences.txt"

OUT_TABLES  <- file.path(OUT_DIR, "virus_informations")
OUT_MARKERS <- file.path(OUT_DIR, "markers")

dir.create(OUT_TABLES, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_MARKERS, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_LOGS, recursive = TRUE, showWarnings = FALSE)

taxonomy_all <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")


# Inputs
manifest <- read_tsv(MANIFEST_TSV, show_col_types = FALSE)

# These files are already cleaned by GENERATE_MANIFEST.
marker_map <- read_tsv(MAP_TSV, show_col_types = FALSE) %>%
  distinct(marker_group_id, group_id, marker)
genome_markers <- read_tsv(GENOME_MAP, show_col_types = FALSE)

# Identify prebuilt RdRp profiles from the cleaned map, not their filenames.
rdrp_ids <- marker_map %>%
  filter(group_id == "orthornavirae", marker == "rdrp") %>%
  pull(marker_group_id)

# Read HMMER domtblout tables
tbl_files <- list.files(HMM_DIR, pattern = "\\.tbl$", full.names = TRUE)

domtblout_cols <- c(
  "orf_name", "target_accession", "target_length_aa",
  "hmm_profile_id", "hmm_accession", "hmm_length_aa",
  "evalue", "score", "bias", "domain_number", "domain_total",
  "c_evalue", "i_evalue", "domain_score", "domain_bias",
  "hmm_from", "hmm_to", "ali_from", "ali_to",
  "env_from", "env_to", "acc"
)

read_domtblout <- function(tbl) {
  marker_group_id_from_file <- basename(tbl) %>% str_remove("\\.tbl$")

  lines <- readLines(tbl, warn = FALSE)
  lines <- lines[!startsWith(lines, "#") & nzchar(trimws(lines))]

  if (length(lines) == 0) return(NULL)

  # Only the first 22 fields are fixed; descriptions may contain spaces.
  fields <- strsplit(trimws(lines[nzchar(trimws(lines))]), "\\s+")
  if (any(lengths(fields) < 22L)) stop("Invalid domtblout table: ", tbl)
  hmmer <- as_tibble(do.call(rbind, lapply(fields, head, 22L)),
                      .name_repair = "minimal")
  names(hmmer) <- domtblout_cols

  # hmmscan reverses the protein and profile columns compared with hmmsearch.
  if (marker_group_id_from_file %in% rdrp_ids) {
    columns <- c("orf_name", "target_length_aa", "hmm_profile_id", "hmm_length_aa")
    hmmer[columns] <- hmmer[columns[c(3, 4, 1, 2)]]
  }

  hmmer %>%
    transmute(
      marker_group_id = marker_group_id_from_file,
      orf_name, hmm_profile_id,
      across(c(target_length_aa, hmm_length_aa, ali_from, ali_to, env_from, env_to), as.integer),
      across(c(evalue, score, i_evalue, domain_score), as.numeric)
    ) %>%
    distinct()
}

# Keep column types even when every search has no hits.
empty_hits <- tibble(
  marker_group_id = character(), orf_name = character(), hmm_profile_id = character(),
  target_length_aa = integer(), hmm_length_aa = integer(),
  ali_from = integer(), ali_to = integer(), env_from = integer(), env_to = integer(),
  evalue = double(), score = double(), i_evalue = double(), domain_score = double()
)
tblout <- bind_rows(empty_hits, map_dfr(tbl_files, read_domtblout))

# Track sequence counts throughout the selection steps
filter_log <- tibble()
log_step <- function(data, step) {
  counts <- data %>%
    group_by(marker_group_id) %>%
    summarise(n_hits = n(), n_orfs = n_distinct(orf_name), .groups = "drop")
  # Retain zero counts when the last ORF for a marker is removed.
  filter_log <<- bind_rows(filter_log, marker_map %>%
    distinct(marker_group_id) %>%
    left_join(counts, by = "marker_group_id") %>%
    mutate(across(c(n_hits, n_orfs), ~replace_na(.x, 0L)), step = step))
}
log_step(tblout, "01_Initial HMMER hits")


# Recover PRODIGAL descriptions from original ORF fasta
# This is required because hmmscan reports the HMM description instead of the ORF description.
read_orf_headers <- function(faa) {
  headers <- readLines(faa, warn = FALSE)
  headers <- headers[str_starts(headers, ">")]

  if (length(headers) == 0) return(tibble(orf_name = character(), description = character(), virus_id = character()))

  tibble(
    orf_name = str_extract(headers, "^>\\S+") %>% str_remove("^>"),
    virus_id = sub("\\.faa$", "", basename(faa)),
    description = str_remove(headers, "^>\\S+\\s*")
  )
}

orf_headers <- list.files(ORF_DIR, pattern = "\\.faa$", full.names = TRUE) %>%
  map_dfr(read_orf_headers)
if (nrow(orf_headers) == 0L) {
  orf_headers <- tibble(orf_name = character(), description = character(), virus_id = character())
}
if (anyDuplicated(orf_headers$orf_name)) {
  stop("Prodigal ORF identifiers must be unique across genomes.")
}

# Contig names may differ from manifest IDs, especially for private genomes.
tblout <- tblout %>%
  left_join(orf_headers, by = "orf_name")
if (anyNA(tblout$virus_id)) stop("HMMER hits contain ORFs absent from the staged proteins.")
if (any(!tblout$virus_id %in% manifest$virus_id)) {
  stop("HMMER hits contain genomes absent from manifest_genomes.tsv.")
}
# Reuse the generated genome/marker associations instead of rebuilding taxonomic rules.
if (nrow(anti_join(tblout, genome_markers, by = c("virus_id", "marker_group_id")))) {
  stop("HMMER hits contain genome/marker pairs absent from genome_marker_map.tsv.")
}


# Remove specified unwanted ORF sequences
rem_seqs <- readLines(REMOVE_SEQS, warn = FALSE)

tblout <- tblout %>%
  filter(!orf_name %in% rem_seqs)

log_step(tblout, "02_Manual removal")


# Keep best hit and credible additional copies
score_threshold_multicopy <- 0.8
length_threshold_multicopy <- 0.8

# First keep the best HMM hit for each ORF.
# This is important for RdRp-scan because one ORF may match several HMM profiles.
# Then rank distinct ORFs within each virus x marker group.
tblout_ranked <- tblout %>%
  group_by(virus_id, marker_group_id, orf_name) %>%
  arrange(desc(score), evalue, desc(domain_score), .by_group = TRUE) %>%  slice(1) %>%
  ungroup() %>%
  group_by(virus_id, marker_group_id) %>%
  arrange(desc(score), evalue, .by_group = TRUE) %>%
  mutate(
    copy_rank = row_number(),
    score_ratio = score / first(score),
    length_ratio = target_length_aa / first(target_length_aa)
  ) %>%
  ungroup()

log_step(tblout_ranked, "03_Best HMM hit per ORF")

# Keep best hit + probable extra copies (except for RdRp)
selected_orfs <- tblout_ranked %>%
  left_join(marker_map, by = "marker_group_id") %>%
  filter(
    copy_rank == 1 |
      (!marker_group_id %in% rdrp_ids &
         score_ratio >= score_threshold_multicopy & length_ratio >= length_threshold_multicopy)
  ) %>%
  left_join(manifest, by = "virus_id")

log_step(selected_orfs, "04_ORF selection")

# Read each selected genome once; both FASTAs use the same Prodigal ORF IDs.
aa_sequences <- vector("list", nrow(selected_orfs))
nt_sequences <- vector("list", nrow(selected_orfs))
selected_orfs$marker_length_AA <- integer(nrow(selected_orfs))
selected_orfs$marker_length_NT <- integer(nrow(selected_orfs))

for (virus in unique(selected_orfs$virus_id)) {
  aa <- Biostrings::readAAStringSet(file.path(ORF_DIR, paste0(virus, ".faa")))
  nt <- Biostrings::readDNAStringSet(file.path(CDS_DIR, paste0(virus, ".ffn")))
  names(aa) <- sub(" .*", "", names(aa))
  names(nt) <- sub(" .*", "", names(nt))

  for (i in which(selected_orfs$virus_id == virus)) {
    r <- selected_orfs[i, ]
    if (!r$orf_name %in% names(aa) || !r$orf_name %in% names(nt)) {
      stop("Missing protein or coding sequence for ORF: ", r$orf_name)
    }
    aa_seq <- aa[r$orf_name]
    nt_seq <- nt[r$orf_name]

    # RdRp exports the domain envelope; nucleotide coordinates are relative to the CDS.
    if (r$marker_group_id %in% rdrp_ids) {
      start <- max(1L, r$env_from)
      end <- min(Biostrings::width(aa_seq), r$env_to)
      aa_seq <- Biostrings::subseq(aa_seq, start = start, end = end)
      nt_seq <- Biostrings::subseq(nt_seq, start = (start - 1L) * 3L + 1L,
                                 end = min(Biostrings::width(nt_seq), end * 3L))
    }

    aa_seq <- Biostrings::AAStringSet(sub("\\*$", "", as.character(aa_seq)))
    aa_sequences[[i]] <- aa_seq
    nt_sequences[[i]] <- nt_seq
    selected_orfs$marker_length_AA[i] <- Biostrings::width(aa_seq)
    selected_orfs$marker_length_NT[i] <- Biostrings::width(nt_seq)
  }
}

# Filter markers shorter than 25% of their own median
selected_orfs <- selected_orfs %>%
  group_by(marker_group_id) %>%
  mutate(min_length_AA = 0.25 * median(marker_length_AA)) %>%
  ungroup()

keep <- selected_orfs$marker_length_AA >= selected_orfs$min_length_AA

selected_orfs <- selected_orfs %>%
  filter(keep) %>%
  select(-min_length_AA)

aa_sequences <- aa_sequences[keep]
nt_sequences <- nt_sequences[keep]

log_step(selected_orfs, "05_Length filter")

filter_log <- filter_log %>%
  arrange(marker_group_id, step) %>%
  group_by(marker_group_id) %>%
  mutate(n_removed_hits = lag(n_hits, default = first(n_hits)) - n_hits,
         n_removed_orfs = lag(n_orfs, default = first(n_orfs)) - n_orfs) %>%
  ungroup()

write_tsv(filter_log, file.path(OUT_LOGS, "filter_summary.tsv"))
print(filter_log, n = Inf)


# Mapping from each filtered table row to its extracted sequence.
selected_orfs$sequence_index <- seq_len(nrow(selected_orfs))

# EXPORT TABLES
marker_group_ids <- sort(unique(selected_orfs$marker_group_id))


# Write table with virus composition across markers + taxonomy
virus_compo_taxo <- selected_orfs %>%
  group_by(virus_id, marker_group_id) %>%
  summarise(orf_names = str_c(sort(unique(orf_name)), collapse = ";"), .groups = "drop") %>%
  pivot_wider(
    names_from = marker_group_id,
    values_from = orf_names,
    names_glue = "{marker_group_id}_{.value}",
    values_fill = NA_character_
  ) %>%
  left_join(
    manifest %>% select(virus_id, ICTV_ID, Virus_names, all_of(taxonomy_all)),
    by = "virus_id"
  )

orf_cols <- str_c(marker_group_ids, "_orf_names")

virus_compo_taxo <- virus_compo_taxo %>%
  select(virus_id, ICTV_ID, Virus_names, all_of(orf_cols), all_of(taxonomy_all))

write_tsv(virus_compo_taxo, file.path(OUT_TABLES, "virus_compo_taxo.tsv"))


# Write table with virus metadata
virus_metadata <- manifest %>%
  select(
    virus_id, ICTV_ID, Virus_names, Virus_names_abrv, Host_source, Source,
    Kingdom, Phylum, Class, Order, Family, Genus, Species
  ) %>%
  rename(Origin_source = Source) %>%
  filter(virus_id %in% virus_compo_taxo$virus_id)

write_tsv(virus_metadata, file.path(OUT_TABLES, "virus_metadata.tsv"))


# Summaries use the selected rows directly; marker names come from the mapfile.
shared_orfs <- selected_orfs %>%
  group_by(orf_name) %>%
  summarise(
    n_marker_groups = n_distinct(marker_group_id),
    marker_group_ids = str_c(sort(unique(marker_group_id)), collapse = ";"),
    .groups = "drop"
  ) %>%
  filter(n_marker_groups > 1)
write_tsv(shared_orfs, file.path(OUT_LOGS, "log_orf_shared_across_markers.tsv"))

family_counts <- selected_orfs %>%
  distinct(Family, marker, virus_id) %>%
  count(Family, marker, name = "n_virus_genomes")
marker_summary <- family_counts %>%
  pivot_wider(names_from = marker, values_from = n_virus_genomes, values_fill = 0) %>%
  arrange(Family)
marker_summary <- bind_rows(marker_summary,
  summarise(marker_summary, across(where(is.numeric), sum), Family = "Total genomes"))
write_lines(knitr::kable(marker_summary, format = "pipe"),
            file.path(OUT_LOGS, "family_marker_summary_wide.md"))

# Facets require at least one selected marker.
if (nrow(selected_orfs) > 0L) {
  p_length <- ggplot(selected_orfs, aes(marker_group_id, marker_length_AA)) +
    geom_boxplot(outlier.size = 0.5) + coord_flip() +
    facet_wrap(~marker, scales = "free_y") +
    labs(x = "Marker group", y = "ORF length", title = "ORF length distributions") +
    theme_bw(base_size = 11)
  ggsave(file.path(OUT_LOGS, "plot_orf_length_boxplot.png"), p_length, width = 12, height = 8)

  p_family <- family_counts %>%
    mutate(Family = fct_reorder(Family, n_virus_genomes, .fun = sum)) %>%
    ggplot(aes(Family, n_virus_genomes, fill = marker)) +
    geom_col(position = "dodge") + coord_flip() +
    facet_wrap(~marker, scales = "free_y") +
    labs(x = "Family", y = "Virus genomes", fill = "Marker",
         title = "Number of virus genomes by family and marker") +
    theme_bw(base_size = 11)
  ggsave(file.path(OUT_LOGS, "plot_family_marker_summary_barplot.png"), p_family, width = 12, height = 8)
}

# EXPORT FINAL FASTA AND PER-MARKER DESCRIPTION TABLES
targets <- selected_orfs %>% distinct(group_id, marker, marker_group_id) %>% arrange(group_id, marker)

for (i in seq_len(nrow(targets))) {
  g <- targets$group_id[i]
  m <- targets$marker[i]
  mg <- targets$marker_group_id[i]

  out_m <- file.path(OUT_MARKERS, g, m)
  dir.create(out_m, recursive = TRUE, showWarnings = FALSE)

  df_m <- selected_orfs %>%
    filter(marker_group_id == mg) %>%
    arrange(virus_id, copy_rank, evalue, desc(score)) %>%
    mutate(description = str_remove(description, "^# ")) %>%
    separate(description, into = c("position_start", "position_end", "prodigal_strand", "prodigal_info"),
             sep = " # ", extra = "merge", fill = "right", convert = TRUE) %>%
    rename(ORF_length_aa = target_length_aa)

  df_m2 <- df_m %>%
    select(orf_name, virus_id, group_id, marker,
           marker_length_AA, marker_length_NT,
           position_start, position_end, ORF_length_aa, ali_from, ali_to, env_from, env_to,
           evalue, score, Species, Virus_names)
  write_tsv(df_m2, file.path(out_m, str_c(mg, "_virus_orf_description.tsv")))

  headers <- df_m %>%
    unite("taxonomy", all_of(taxonomy_all), sep = ";") %>%
    transmute(header = str_c(orf_name, " ", marker, " ", taxonomy)) %>%
    pull(header)
  aa_m <- do.call(c, aa_sequences[df_m$sequence_index])
  nt_m <- do.call(c, nt_sequences[df_m$sequence_index])
  names(aa_m) <- headers
  names(nt_m) <- headers

  Biostrings::writeXStringSet(aa_m, file.path(out_m, str_c(mg, "_protein.fasta")), format = "fasta", width = 60)
  Biostrings::writeXStringSet(nt_m, file.path(out_m, str_c(mg, "_nucleotide.fasta")), format = "fasta", width = 60)
}
