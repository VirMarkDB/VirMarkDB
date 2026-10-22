library(tidyverse)
library(rentrez)

guglielmini_dir <- "output/reference_sources/guglielmini"
guglielmini_files_path <- str_c(guglielmini_dir, "/Alignments/Sequences")
manual_ids_protein_path <- "input/manual_reference_protein"
merged_sequences_path <- "output/references_protein"
hmm_dir <- "output/hmm/hmms"

dir.create(guglielmini_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(merged_sequences_path, recursive = TRUE, showWarnings = FALSE)
dir.create(hmm_dir, recursive = TRUE, showWarnings = FALSE)


# Definition of marker groups
nucleocytoviricota_markers <- c("majorcapsid_nucleocytoviricota",
  "atpase_nucleocytoviricota",      
  "dnapol_nucleocytoviricota",
  "primase_nucleocytoviricota",
  "rnapol1_nucleocytoviricota",
  "rnapol2_nucleocytoviricota",
  "tf2s_nucleocytoviricota",
  "vltf3_nucleocytoviricota")

pantevenvirales_markers <- "g23_pantevenvirales"

## For Nucleocytoviricota and Polinton like viruses

# Download and unzip guglielmini references
zip_url  <- "https://zenodo.org/records/3368642/files/Additional%20data.zip?download=1"
zip_dest <- file.path(guglielmini_dir, "Additional_data.zip")

if (!file.exists(zip_dest)) {
  download.file(zip_url, destfile = zip_dest, mode = "wb")
  unzip(zip_dest, exdir = guglielmini_dir)
}


# Load guglielmini protein references
majorcapsid_nucleocytoviricota_guglielmini <- Biostrings::readAAStringSet(file.path(guglielmini_files_path, "capsid.prt"))
atpase_nucleocytoviricota_guglielmini       <- Biostrings::readAAStringSet(file.path(guglielmini_files_path, "pATPase.prt"))
dnapol_nucleocytoviricota_guglielmini       <- Biostrings::readAAStringSet(file.path(guglielmini_files_path, "dnapol.prt"))
primase_nucleocytoviricota_guglielmini      <- Biostrings::readAAStringSet(file.path(guglielmini_files_path, "primase.prt"))
rnapol1_nucleocytoviricota_guglielmini      <- Biostrings::readAAStringSet(file.path(guglielmini_files_path, "rnapol1.prt"))
rnapol2_nucleocytoviricota_guglielmini      <- Biostrings::readAAStringSet(file.path(guglielmini_files_path, "rnapol2.prt"))
tf2s_nucleocytoviricota_guglielmini         <- Biostrings::readAAStringSet(file.path(guglielmini_files_path, "tf2s.prt"))
vltf3_nucleocytoviricota_guglielmini        <- Biostrings::readAAStringSet(file.path(guglielmini_files_path, "vltf3.prt"))

atpase_polinto_guglielmini                  <- Biostrings::readAAStringSet(file.path(guglielmini_files_path, "pATPase-polinto.prt"))
majorcapsid_polinto_guglielmini             <- Biostrings::readAAStringSet(file.path(guglielmini_files_path, "capsid-polinto.prt"))



# For pantevenvirales

# Load G23 protein accessions from Sieradzki et al. 2019
g23_pantevenvirales_sieradzki_path <- "input/reference_data/sieradzki_et_al_2019_g23.tsv"
g23_pantevenvirales_sieradzki_outpath <- "output/reference_sources/sieradzki_g23.fasta"

g23_pantevenvirales_sieradzki_tsv <- read_tsv(
  g23_pantevenvirales_sieradzki_path,
  col_names = c("accession", "description"),
  show_col_types = FALSE
)

g23_pantevenvirales_sieradzki_txt <- entrez_fetch(
  db = "protein",
  id = g23_pantevenvirales_sieradzki_tsv$accession,
  rettype = "fasta",
  retmode = "text"
)

writeLines(
  g23_pantevenvirales_sieradzki_txt,
  g23_pantevenvirales_sieradzki_outpath
)

g23_pantevenvirales_sieradzki <- Biostrings::readAAStringSet(
  g23_pantevenvirales_sieradzki_outpath
)


# Load G23 proteins from Sullivan et al. 2010
g23_pantevenvirales_sullivan_path <- "input/reference_data/sullivan_et_al_2010_g23.fasta"

g23_pantevenvirales_sullivan <- Biostrings::readAAStringSet(
  g23_pantevenvirales_sullivan_path
)

g23_pantevenvirales_sullivan <- g23_pantevenvirales_sullivan[
  grep(
    "gp23 precursor of major head subunit",
    names(g23_pantevenvirales_sullivan),
    ignore.case = TRUE
  )
]

# Merge G23 sources
g23_pantevenvirales_sieradzki_sullivan <- c(g23_pantevenvirales_sieradzki,g23_pantevenvirales_sullivan)

# Load manual accessions and fetch NCBI protein sequences
manual_markers <- list.files(
  manual_ids_protein_path,
  pattern = "\\.tsv$",
  full.names = FALSE
) %>%
  tools::file_path_sans_ext()


for (marker in manual_markers) {
  print(marker)
  manual_file <- str_c(manual_ids_protein_path, "/", marker, ".tsv")

  if (file.info(manual_file)$size == 0) {
    assign(str_c(marker, "_manual"), tibble(
      accession = character(),
      description = character()
    ))
    assign(str_c(marker, "_seq"), Biostrings::AAStringSet())
    next
  }
  
  x <- manual_file |>
    read_tsv(col_names = c("accession", "description"),show_col_types = FALSE)

  assign(str_c(marker, "_manual"), x)
  
  if (nrow(x) == 0) {
    assign(str_c(marker, "_seq"), Biostrings::AAStringSet())
    next
  }
  
  fasta_txt <- entrez_fetch(
    db = "protein",
    id = x$accession,
    rettype = "fasta",
    retmode = "text"
  )
  
  tmp_fasta <- tempfile(fileext = ".fasta")
  writeLines(fasta_txt, tmp_fasta)
  
  seq_obj <- Biostrings::readAAStringSet(tmp_fasta)
  
  assign(str_c(marker, "_seq"), seq_obj)
}

# Merge Guglielmini / g23 (sieradzki and sullivan) + manual sequences, remove duplicates, write FASTA
all_markers <- c(nucleocytoviricota_markers,pantevenvirales_markers)

for (marker in all_markers) {
  
  if (marker %in% nucleocytoviricota_markers){
      reference_name <- str_c(marker, "_guglielmini")
  }
  if (marker %in% pantevenvirales_markers){
      reference_name <- str_c(marker, "_sieradzki_sullivan")
  }
  

  manual_name <- str_c(marker, "_seq")
  
  reference_seq <- if (exists(reference_name, inherits = FALSE)) get(reference_name) else Biostrings::AAStringSet()
  manual_seq <- if (exists(manual_name, inherits = FALSE)) get(manual_name) else Biostrings::AAStringSet()
  
  merged_seq <- c(reference_seq, manual_seq)
  merged_seq <- merged_seq[!duplicated(as.character(merged_seq))]
  
  assign(marker, merged_seq)
  
  Biostrings::writeXStringSet(
    merged_seq,
    filepath = file.path(merged_sequences_path, str_c(marker, ".fasta")),
    format = "fasta"
  )
  
  cat(marker, ":", length(merged_seq), "sequences written\n")
}



## For Orthornavirae, use the RdRp-scan HMM database

rdrpscan_dir <- file.path(hmm_dir, "rdrpscan")
dir.create(rdrpscan_dir, recursive = TRUE, showWarnings = FALSE)

rdrpscan_base_url <- str_c(
  "https://raw.githubusercontent.com/JustineCharon/",
  "RdRp-scan/06b514d79c77d34b240bfab74f04ed42d89808f8/",
  "Profile_db_and_alignments/"
)

rdrpscan_db_name <- "RdRp_HMM_profile_CLUSTALO.db"

for (ext in c("h3f", "h3i", "h3m", "h3p")) {

  file_name <- str_c(rdrpscan_db_name, ".", ext)
  dest_file <- file.path(rdrpscan_dir, file_name)

  if (!file.exists(dest_file)) {
    download.file(
      str_c(rdrpscan_base_url, file_name),
      destfile = dest_file,
      mode = "wb"
    )
  }
}

rdrpscan_db <- file.path(
  rdrpscan_dir,
  rdrpscan_db_name
)





