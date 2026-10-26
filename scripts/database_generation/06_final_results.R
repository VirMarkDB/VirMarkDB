library(tidyverse)

# CONFIG
HMM_DIR         <- "output/hmm/search"
ORF_DIR         <- "output/orfs"
REF_PROTEIN_DIR <- "output/references_protein"
MANIFEST_TSV    <- "output/config/manifest_genomes.tsv"
MAP_TSV         <- "output/config/marker_taxo_map.tsv"
OUT_DIR         <- "output/VirMarkDB"
OUT_LOGS        <- "output/hmm/logs"
REMOVE_SEQS     <- "input/remove_specific_ORF_sequences.txt"

OUT_TABLES  <- file.path(OUT_DIR, "virus_informations")
OUT_MARKERS <- file.path(OUT_DIR, "markers")

dir.create(OUT_TABLES, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_MARKERS, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_LOGS, recursive = TRUE, showWarnings = FALSE)

taxonomy_all <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")


# Inputs
manifest <- read_tsv(MANIFEST_TSV, show_col_types = FALSE)

marker_map <- read_tsv(MAP_TSV, show_col_types = FALSE) %>%
  mutate(marker_group_id = if_else(is.na(marker_group_id) | marker_group_id == "",
                                   str_c(group_id, "__", marker), marker_group_id)) %>%
  distinct(marker_group_id, group_id, marker)


# Read hmmer DOMBLOUT tables
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
  lines <- lines[!startsWith(lines, "#")]

  if (length(lines) == 0) return(tibble())

  hmmer_raw <- readr::read_table(
    paste0(paste(lines, collapse = "\n"), "\n"),
    col_names = FALSE,
    show_col_types = FALSE,
    progress = FALSE
  )

  hmmer <- hmmer_raw %>% select(1:22)
  names(hmmer) <- domtblout_cols

  # RdRp-scan is run with hmmscan:
  # target = HMM and query = ORF.
  # Swap target/query fields to recover the standard VirMarkDB structure.
  if (marker_group_id_from_file == "rdrp_orthornavirae") {
    target_name <- hmmer$orf_name
    target_accession <- hmmer$target_accession
    target_length <- hmmer$target_length_aa

    query_name <- hmmer$hmm_profile_id
    query_accession <- hmmer$hmm_accession
    query_length <- hmmer$hmm_length_aa

    hmmer$orf_name <- query_name
    hmmer$target_accession <- query_accession
    hmmer$target_length_aa <- query_length

    hmmer$hmm_profile_id <- target_name
    hmmer$hmm_accession <- target_accession
    hmmer$hmm_length_aa <- target_length
  }

  hmmer <- hmmer %>%
    mutate(
      marker_group_id = marker_group_id_from_file,
      virus_id = str_remove(orf_name, "_[0-9]+$"),
      orf_index = str_extract(orf_name, "[0-9]+$"),
      .before = 1
    ) %>%
    mutate(
      across(c(target_length_aa, hmm_length_aa, domain_number, domain_total,
               hmm_from, hmm_to, ali_from, ali_to, env_from, env_to, orf_index), as.integer),
      across(c(evalue, score, bias, c_evalue, i_evalue, domain_score, domain_bias, acc), as.numeric),
      virus_id_cut = str_remove(virus_id, "\\.[0-9]+$")
    ) %>%
    select(
      virus_id, virus_id_cut, orf_name, target_length_aa,
      marker_group_id, hmm_profile_id, hmm_length_aa,
      evalue, score, i_evalue, domain_score,
      ali_from, ali_to, env_from, env_to
    ) %>%
    distinct()

  return(hmmer)
}

tblout <- map_dfr(tbl_files, read_domtblout)

# Track sequence counts throughout the selection steps
filter_log <- tibble()
log_step <- function(data, step) {
  filter_log <<- bind_rows(filter_log, data %>%
    group_by(marker_group_id) %>%
    summarise(n_hits = n(), n_orfs = n_distinct(orf_name), .groups = "drop") %>%
    mutate(step = step))
}
log_step(tblout, "01_Initial HMMER hits")


# Recover PRODIGAL descriptions from original ORF fasta
# This is required because hmmscan reports the HMM description instead of the ORF description.
read_orf_headers <- function(faa) {
  headers <- readLines(faa, warn = FALSE)
  headers <- headers[str_starts(headers, ">")]

  if (length(headers) == 0) return(tibble())

  tibble(
    orf_name = str_extract(headers, "^>\\S+") %>% str_remove("^>"),
    description = str_remove(headers, "^>\\S+\\s*")
  )
}

orf_headers <- list.files(ORF_DIR, pattern = "\\.faa$", full.names = TRUE) %>%
  map_dfr(read_orf_headers) %>%
  distinct(orf_name, .keep_all = TRUE)

tblout <- tblout %>% left_join(orf_headers, by = "orf_name")


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
    best_score = first(score),
    best_orf_length = first(target_length_aa),
    score_ratio = score / best_score,
    length_ratio = target_length_aa / best_orf_length
  ) %>%
  rename(virus_id_tblout = virus_id) %>%
  ungroup()

log_step(tblout_ranked, "03_Best HMM hit per ORF")

# Manage virus id in manifest
manifest <- manifest %>%
  mutate(virus_id_cut = str_remove(virus_id, "_partial"))


# Keep best hit + probable extra copies (except for RdRp)
selected_orfs <- tblout_ranked %>%
  left_join(marker_map, by = "marker_group_id") %>%
  filter(
    (marker_group_id == "rdrp_orthornavirae" & copy_rank == 1) |
    (marker_group_id != "rdrp_orthornavirae" &
       (copy_rank == 1 | (score_ratio >= score_threshold_multicopy & length_ratio >= length_threshold_multicopy)))
  ) %>%
  left_join(manifest, by = "virus_id_cut")

log_step(selected_orfs, "04_ORF selection")

# LOAD ORFS
load_orfs <- function(virus_id) {
  faa_path <- file.path(ORF_DIR, str_c(virus_id, ".faa"))
  fna_path <- file.path(ORF_DIR, str_c(virus_id, ".fna"))

  aa_seq <- if (file.exists(faa_path)) Biostrings::readAAStringSet(faa_path) else Biostrings::AAStringSet()
  nt_seq <- if (file.exists(fna_path)) Biostrings::readDNAStringSet(fna_path) else Biostrings::DNAStringSet()

  names(aa_seq) <- sub(" .*", "", names(aa_seq))
  names(nt_seq) <- sub(" .*", "", names(nt_seq))

  list(aa = aa_seq, nt = nt_seq)
}

ORF_CACHE <- list()

get_orfs <- function(virus_id) {
  if (is.null(ORF_CACHE[[virus_id]])) ORF_CACHE[[virus_id]] <<- load_orfs(virus_id)
  ORF_CACHE[[virus_id]]
}


# Remove terminal star
strip_terminal_star <- function(x) {
  x_chr <- as.character(x)
  x_chr2 <- sub("\\*$", "", x_chr)
  names(x_chr2) <- names(x)
  Biostrings::AAStringSet(x_chr2)
}


# Extract final marker sequebces and measure lengths
# The median threshold is applied to the sequences that will be exported.
aa_sequences <- vector("list", nrow(selected_orfs))
nt_sequences <- vector("list", nrow(selected_orfs))
selected_orfs$marker_length_AA <- NA_integer_
selected_orfs$marker_length_NT <- NA_integer_

for (i in seq_len(nrow(selected_orfs))) {
  r <- selected_orfs[i, ]
  orfs <- get_orfs(r$virus_id)
  aa_seq <- orfs$aa[r$orf_name]
  nt_seq <- orfs$nt[r$orf_name]

  # RdRp: export only the HMMER domain envelope, not the whole ORF.
  if (r$marker_group_id == "rdrp_orthornavirae") {
    aa_start <- max(1L, r$env_from)
    aa_end <- min(Biostrings::width(aa_seq), r$env_to)

    nt_start <- (aa_start - 1L) * 3L + 1L
    nt_end <- min(Biostrings::width(nt_seq), aa_end * 3L)

    aa_seq <- Biostrings::subseq(aa_seq, start = aa_start, end = aa_end)
    nt_seq <- Biostrings::subseq(nt_seq, start = nt_start, end = nt_end)
  }

  aa_seq <- strip_terminal_star(aa_seq)

  aa_sequences[[i]] <- aa_seq
  nt_sequences[[i]] <- nt_seq
  selected_orfs$marker_length_AA[i] <- Biostrings::width(aa_seq)
  selected_orfs$marker_length_NT[i] <- Biostrings::width(nt_seq)
}

# Filter markers shorter than 25% of their own median
selected_orfs <- selected_orfs %>%
  group_by(marker_group_id) %>%
  mutate(median_length_AA = median(marker_length_AA), 
         min_length_AA = 0.25 * median_length_AA) %>%
  ungroup()

keep <- selected_orfs$marker_length_AA >= selected_orfs$min_length_AA

selected_orfs <- selected_orfs %>%
  filter(keep) %>%
  select(-median_length_AA, -min_length_AA)

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


# Calculate logs after filtering
source("scripts/utils/hmm_search_logs.R")

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

  aa_m <- Biostrings::AAStringSet()
  nt_m <- Biostrings::DNAStringSet()

  for (j in seq_len(nrow(df_m))) {
    r <- df_m[j, ]
    hdr <- str_c(r$orf_name, " ", m, " ",
                 df_m[j, ] %>% unite("Taxonomy", Kingdom, Phylum, Class, Order, Family, Genus, Species, sep = ";") %>% pull(Taxonomy))
                 
    aa_seq <- aa_sequences[[r$sequence_index]]
    nt_seq <- nt_sequences[[r$sequence_index]]

    names(aa_seq) <- hdr
    names(nt_seq) <- hdr

    aa_m <- c(aa_m, aa_seq)
    nt_m <- c(nt_m, nt_seq)
  }

  Biostrings::writeXStringSet(aa_m, file.path(out_m, str_c(mg, "_protein.fasta")), format = "fasta", width = 60)
  Biostrings::writeXStringSet(nt_m, file.path(out_m, str_c(mg, "_nucleotide.fasta")), format = "fasta", width = 60)
}
