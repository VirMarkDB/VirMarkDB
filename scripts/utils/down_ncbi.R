library(tidyverse)
library(rentrez)

# Config
manifest <- "output/config/manifest_genomes.tsv"
outdir   <- "output/genomes"

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# NCBI credentials
creds <- read_lines("input/ncbi.creds")
options(entrez_email = creds[1])
options(entrez_key   = creds[2])

# Load manifest and filter out private genomes and alrady downloaded genomes
genomes <- read_tsv(manifest, show_col_types = FALSE) %>%
  filter(Source != "PRIVATE") %>%
  filter(!file.exists(file.path(outdir, str_c(virus_id, ".fna"))))


# Function to download genomes
# Log failed downloads
log_file <- str_c(outdir, "/ncbi_download_failures.tsv")
if (!file.exists(log_file)) {write_lines("virus_id\tncbi_id\terror", log_file)}

# Function to download genomes
download_genome <- function(virus_id) {
  outfile <- file.path(outdir, paste0(virus_id, ".fna"))
  ncbi_id <- gsub("_partial", "", virus_id)

  if (file.exists(outfile)) {
    return()
  }

  tryCatch({
    seq <- rentrez::entrez_fetch(
      db = "nuccore",
      id = ncbi_id,
      rettype = "fasta",
      retmode = "text"
    )

    write_lines(seq, outfile)
    message("Downloaded: ", virus_id)

  }, error = function(e) {
    msg <- str_replace_all(e$message, "[\r\n\t]+", " ")

    warning("Skipping ", virus_id, " / ", ncbi_id, ": ", msg)

    write_lines(
      str_c(virus_id, "\t", ncbi_id, "\t", msg),
      log_file,
      append = TRUE
    )

    return(NULL)
  })
}

# Apply funtion to download the genomes
genomes %>%
  select(virus_id) %>%
  pwalk(download_genome)
