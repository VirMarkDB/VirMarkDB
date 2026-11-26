library(readr)
library(rentrez)
library(Biostrings)

# Entrées fournies par Nextflow : correspondance des marqueurs et deux dossiers.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript collect_ref.R <marker_map.tsv> <manual_references> <reference_data>")
}
markers <- read_tsv(args[[1]], show_col_types = FALSE)
markers <- unique(markers[c("marker_group_id", "group_id", "marker", "manual_file")])
manual_dir <- args[[2]]
reference_data <- args[[3]]

dir.create("reference_sources", showWarnings = FALSE)
dir.create("references_protein", showWarnings = FALSE)

# Lit une liste d'accessions et récupère leurs protéines NCBI.
# Un fichier manuel absent ou vide signifie qu'il n'y a rien à ajouter.
fetch_proteins <- function(accession_file, fasta_file) {
  if (!file.exists(accession_file) || file.info(accession_file)$size == 0) {
    return(AAStringSet())
  }
  accessions <- read_tsv(accession_file, col_names = c("accession", "description"),
                         col_types = "cc")
  if (!nrow(accessions)) return(AAStringSet())
  # Trois essais avec une attente croissante, notamment si NCBI répond 429.
  for (attempt in 1:3) {
    fasta <- tryCatch(
      entrez_fetch(db = "protein", id = accessions$accession,
                   rettype = "fasta", retmode = "text"),
      error = function(e) e
    )
    Sys.sleep(0.4)
    if (!inherits(fasta, "error")) break
    if (attempt == 3L) stop(fasta)
    message("NCBI indisponible pour ", basename(accession_file), ", nouvel essai...")
    Sys.sleep(2^attempt)
  }
  writeLines(fasta, fasta_file)
  readAAStringSet(fasta_file)
}

# Les noms dans l'archive publiée diffèrent parfois des noms du mapfile.
guglielmini_files <- c(
  majorcapsid = "capsid.prt", atpase = "pATPase.prt", dnapol = "dnapol.prt",
  primase = "primase.prt", rnapol1 = "rnapol1.prt", rnapol2 = "rnapol2.prt",
  tf2s = "tf2s.prt", vltf3 = "vltf3.prt"
)
if (any(markers$group_id == "nucleocytoviricota")) {
  guglielmini_dir <- "reference_sources/guglielmini"
  dir.create(guglielmini_dir, recursive = TRUE, showWarnings = FALSE)
  archive <- file.path(guglielmini_dir, "Additional_data.zip")
  download.file("https://zenodo.org/records/3368642/files/Additional%20data.zip?download=1",
                destfile = archive, mode = "wb")
  unzip(archive, exdir = guglielmini_dir)
}

# Une ligne du mapfile = une référence, identifiée par marker_group_id.
for (i in seq_len(nrow(markers))) {
  ref_id <- markers$marker_group_id[[i]]
  group <- markers$group_id[[i]]
  marker <- markers$marker[[i]]
  sequences <- AAStringSet()

  # RdRp-scan fournit déjà des profils : pas de FASTA à construire ici.
  if (group == "orthornavirae" && marker == "rdrp") {
    dir.create("hmms/rdrpscan", recursive = TRUE, showWarnings = FALSE)
    base_url <- paste0(
      "https://raw.githubusercontent.com/JustineCharon/",
      "RdRp-scan/06b514d79c77d34b240bfab74f04ed42d89808f8/Profile_db_and_alignments/"
    )
    for (ext in c("h3f", "h3i", "h3m", "h3p")) {
      filename <- paste0("RdRp_HMM_profile_CLUSTALO.db.", ext)
      download.file(paste0(base_url, filename), file.path("hmms/rdrpscan", filename), mode = "wb")
    }
    next
  }

  # Source publiée pour les marqueurs Nucleocytoviricota.
  if (group == "nucleocytoviricota") {
    filename <- guglielmini_files[[marker]]
    sequences <- readAAStringSet(file.path(guglielmini_dir, "Alignments/Sequences", filename))
  }

  # Deux sources publiées pour G23 : Sieradzki et Sullivan.
  if (group == "pantevenvirales" && marker == "g23") {
    sequences <- fetch_proteins(
      file.path(reference_data, "sieradzki_et_al_2019_g23.tsv"),
      "reference_sources/sieradzki_g23.fasta"
    )
    sullivan <- readAAStringSet(file.path(reference_data, "sullivan_et_al_2010_g23.fasta"))
    sullivan <- sullivan[grepl("gp23 precursor of major head subunit", names(sullivan), ignore.case = TRUE)]
    sequences <- c(sequences, sullivan)
  }

  # Ajout des protéines manuelles, en utilisant le nom de fichier du mapfile.
  manual_file <- markers$manual_file[[i]]
  if (!is.na(manual_file) && nzchar(manual_file)) {
    manual <- fetch_proteins(file.path(manual_dir, manual_file),
                             file.path("reference_sources", paste0(ref_id, "_manual.fasta")))
    sequences <- c(sequences, manual)
  }

  # Une seule copie de chaque séquence, puis un FASTA portant l'identifiant de référence.
  sequences <- sequences[!duplicated(as.character(sequences))]
  if (!length(sequences)) stop("Aucune protéine de référence pour ", ref_id)
  writeXStringSet(sequences, file.path("references_protein", paste0(ref_id, ".fasta")))
  message(ref_id, " : ", length(sequences), " protéines")
}
