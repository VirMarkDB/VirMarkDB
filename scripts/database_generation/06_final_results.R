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


# INPUTS
manifest <- read_tsv(MANIFEST_TSV, show_col_types = FALSE)

marker_map <- read_tsv(MAP_TSV, show_col_types = FALSE) %>%
  mutate(marker_group_id = if_else(is.na(marker_group_id) | marker_group_id == "",
                                   str_c(group_id, "__", marker), marker_group_id)) %>%
  distinct(marker_group_id, group_id, marker)


# READ HMMER DOMTBLOUT TABLES
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


# RECOVER PRODIGAL DESCRIPTIONS FROM ORIGINAL ORF FASTA
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


# Keep best hit and credible additional copies
score_threshold_multicopy <- 0.8
length_threshold_multicopy <- 0.8

# First keep the best HMM/domain hit for each ORF.
# This is important for RdRp-scan because one ORF may match several HMM profiles.
# Then rank distinct ORFs within each virus x marker group.
tblout_ranked <- tblout %>%
  group_by(virus_id, marker_group_id, orf_name) %>%
  arrange(desc(domain_score), i_evalue, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  group_by(virus_id, marker_group_id) %>%
  arrange(desc(domain_score), i_evalue, .by_group = TRUE) %>%
  mutate(
    copy_rank = row_number(),
    best_score = first(domain_score),
    best_orf_length = first(target_length_aa),
    score_ratio = domain_score / best_score,
    length_ratio = target_length_aa / best_orf_length
  ) %>%
  rename(virus_id_tblout = virus_id) %>%
  ungroup()


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


# Calculate logs
source("scripts/utils/hmm_search_logs.R")


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


# EXPORT FASTA BY GROUP / MARKER
targets <- selected_orfs %>%
  distinct(group_id, marker, marker_group_id) %>%
  arrange(group_id, marker)

for (i in seq_len(nrow(targets))) {
  g  <- targets$group_id[i]
  m  <- targets$marker[i]
  mg <- targets$marker_group_id[i]

  out_m <- file.path(OUT_MARKERS, g, m)
  dir.create(out_m, recursive = TRUE, showWarnings = FALSE)

  df_m <- selected_orfs %>%
    filter(marker_group_id == mg) %>%
    arrange(virus_id, copy_rank, evalue, desc(score)) %>%
    mutate(description = str_remove(description, "^# ")) %>%
    separate(
      description,
      into = c("position_start", "position_end", "prodigal_strand", "prodigal_info"),
      sep = " # ",
      extra = "merge",
      fill = "right",
      convert = TRUE
    )

  df_m2 <- df_m %>%
    select(
      orf_name, virus_id, group_id, marker,
      evalue, score, target_length_aa,
      position_start, position_end,
      ali_from, ali_to, env_from, env_to,
      Species, Virus_names
    ) %>%
    rename(ORF_length_aa = target_length_aa)

  write_tsv(df_m2, file.path(out_m, str_c(mg, "_virus_orf_description.tsv")))

  aa_m <- Biostrings::AAStringSet()
  nt_m <- Biostrings::DNAStringSet()

  for (j in seq_len(nrow(df_m))) {
    r <- df_m[j, ]

    orfs <- get_orfs(r$virus_id)

    hdr <- str_c(
      r$orf_name, " ", m, " ",
      df_m[j, ] %>%
        unite("Taxonomy", Kingdom, Phylum, Class, Order, Family, Genus, Species, sep = ";") %>%
        pull(Taxonomy)
    )

    aa_seq <- orfs$aa[r$orf_name]
    nt_seq <- orfs$nt[r$orf_name]

    # For RdRp, extract only the region detected by RdRp-scan
    if (r$marker_group_id == "rdrp_orthornavirae") {
      aa_start <- max(1, r$env_from)
      aa_end <- min(Biostrings::width(aa_seq), r$env_to)

      nt_start <- max(1, ((aa_start - 1) * 3) + 1)
      nt_end <- min(Biostrings::width(nt_seq), aa_end * 3)

      aa_seq <- Biostrings::subseq(aa_seq, start = aa_start, end = aa_end)
      nt_seq <- Biostrings::subseq(nt_seq, start = nt_start, end = nt_end)
    }

    aa_m <- c(aa_m, aa_seq)
    aa_m <- strip_terminal_star(aa_m)
    names(aa_m)[length(aa_m)] <- hdr

    nt_m <- c(nt_m, nt_seq)
    names(nt_m)[length(nt_m)] <- hdr
  }

  Biostrings::writeXStringSet(
    aa_m,
    file.path(out_m, str_c(mg, "_protein.fasta")),
    format = "fasta",
    width = 60
  )

  Biostrings::writeXStringSet(
    nt_m,
    file.path(out_m, str_c(mg, "_nucleotide.fasta")),
    format = "fasta",
    width = 60
  )
}
