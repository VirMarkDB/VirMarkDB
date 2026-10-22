library(tidyverse)

# Config
ictv_xlsx <- "VMR_MSL41.v1.20260320.xlsx"
manual_gen_ncbi <- "input/manual_genomes_ncbi.tsv"
manual_gen_priv <- "input/manual_genomes_private.tsv"
map_file <- "input/mapfile.tsv"
out_tsv <- "output/config/manifest_genomes.tsv"
genomes_private_path <- "input/private_genomes"
out_genomes <- "output/genomes"

dir.create("output")
dir.create("output/config", recursive = TRUE)
dir.create(out_genomes, recursive = TRUE)

taxo_cols <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Read mapping
marker_map <- read_tsv(map_file, show_col_types = FALSE) %>%
  mutate(
    marker_group_id = if_else(
      is.na(marker_group_id) | marker_group_id == "",
      str_c(group_id, "__", marker),
      marker_group_id
    )
  ) %>%
  distinct()

wanted_taxa <- marker_map %>%
  distinct(taxo_rank, taxo_value)


# Handle ICTV standards
url <- str_c("https://ictv.global/sites/default/files/VMR/", ictv_xlsx)
dest <- tempfile(fileext = ".xlsx")
download.file(url, destfile = dest, mode = "wb")

ictv <- readxl::read_excel(dest, sheet = 2) %>%
  rename(
    virus_id = `Virus GENBANK accession`,
    Virus_names = `Virus name(s)`,
    Virus_names_abrv = `Virus name abbreviation(s)`,
    Host_source = `Host source`
  ) %>%
  select(
    virus_id, taxo_cols, Virus_names, Virus_names_abrv, Host_source, ICTV_ID
  ) %>%
  separate_rows(virus_id, sep = "\\s*;\\s*") %>%
  mutate(
    Source = "ICTV"
  ) %>%
  filter(!str_detect(virus_id, "\\(\\d+\\.\\d+\\)")) %>%
  mutate(
    virus_id = virus_id %>%
      str_replace("^partial\\s*:\\s*(.+)$", "\\1_partial") %>%
      str_remove("^.*:\\s*")
  )

# Load manual genomes on NCBI
manual_ncbi <- if (file.exists(manual_gen_ncbi)) {
  read_tsv(manual_gen_ncbi, show_col_types = FALSE)
  } else {
  tibble()
}

# Load manual private genomes
manual_priv <- if (file.exists(manual_gen_priv)) {
  read_tsv(manual_gen_priv, show_col_types = FALSE) 
} else {
  tibble()
}

# Merge
genomes_pre <- bind_rows(ictv, manual_ncbi, manual_priv) %>%
  distinct(virus_id, .keep_all = TRUE)

genomes <- tibble()

for(i in seq(nrow(wanted_taxa))){
  print(i)
  genomes <- rbind(
    genomes,
    genomes_pre %>%
      filter(.data[[wanted_taxa$taxo_rank[[i]]]] == wanted_taxa$taxo_value[[i]])
  )
}
genomes <- genomes %>%
  distinct()

# Check duplication
genomes %>%
  count(virus_id, name = "n") %>%
  filter(n > 1) %>%
  {
    if (nrow(.)) {
      write_tsv(., "output/config/dup_ids_check.tsv")
    }
  }

# Export manifest
write_tsv(genomes, out_tsv)

# Export cleaned map
write_tsv(marker_map, "output/config/marker_taxo_map.tsv")

# Genome -> marker_group_id
genome_marker_map <- tibble()

for (i in seq_len(nrow(marker_map))) {
  genome_marker_map <- bind_rows(
    genome_marker_map,
    genomes %>%
      filter(.data[[marker_map$taxo_rank[[i]]]] == marker_map$taxo_value[[i]]) %>%
      mutate(
        group_id = marker_map$group_id[[i]],
        marker = marker_map$marker[[i]],
        marker_group_id = marker_map$marker_group_id[[i]]
      ) %>%
      select(
        virus_id, taxo_cols, group_id, marker, marker_group_id
      )
  )
}

genome_marker_map <- genome_marker_map %>%
  distinct()

write_tsv(genome_marker_map, "output/config/genome_marker_map.tsv")



# Manage DNA files 

# Download genomes by accession
source("scripts/utils/down_ncbi.R") 

# Merge genomes sources
private_files <- list.files(genomes_private_path, full.names = TRUE)


if (length(private_files) > 0) {
  for (file in private_files) {
    file_name <- basename(file) 
    file.copy(from = file, to = file.path(out_genomes, file_name), overwrite = FALSE)
  }
}


# Manage references protein
source("scripts/utils/prepare_ref_protein.R")
