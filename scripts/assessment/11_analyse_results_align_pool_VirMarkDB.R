library(tidyverse)
library(patchwork)

ncbi_dir <- "output/assessment/pool_protein/pool_protein_filt"
hits_dir <- "output/assessment/alignment_VirMarkDB_pool/alignment"
out_dir  <- "output/assessment/alignment_VirMarkDB_pool/bench_results"
missing_prot_path <- "output/assessment/alignment_VirMarkDB_pool/missing_prot"
hits_prot_path <- "output/assessment/alignment_VirMarkDB_pool/hits_prot"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(missing_prot_path, recursive = TRUE, showWarnings = FALSE)
dir.create(hits_prot_path, recursive = TRUE, showWarnings = FALSE)

# read compo viro
meta <- read_tsv("output/VirMarkDB/virus_informations/virus_compo_taxo.tsv")


analyse_res <- function(marker,ncbi_dir,hits_dir,out_dir,missing_prot_path,hits_prot_path,meta){ 
  prot <- Biostrings::readAAStringSet(str_c(ncbi_dir,"/", marker,".fasta"))

  tab <- read.table(str_c(hits_dir,"/",marker,".m8"))
  colnames(tab) <- c("query","ref","ident","alnlen","mismatch","gapopen","qstart","qend","sstart","send","evalue","bits") 

  # Select best hit
  tab_best <- tab |>
    group_by(query)  |>
    slice_min(evalue, with_ties = FALSE) 

  # Write best match results
  write_tsv(tab_best, str_c(out_dir,"/",marker,"_best_hits.tsv"))

  # Let's check missing prot query 
  prot_id <- sub(" .*", "", names(prot))
  missing_prot_id <- setdiff(unique(prot_id), unique(tab$query))
  prot_missing <- prot[prot_id %in% missing_prot_id]
  

  # write in fasta missing prot 
  if (length(prot_missing) > 0) {
    Biostrings::writeXStringSet(prot_missing,str_c(missing_prot_path,"/",marker,".fasta"))
  }

  # Export hits prot 

  prot_hits <- prot[prot_id %in% tab_best$query]
  # write in fasta 
  if (length(prot_hits) > 0) {
    Biostrings::writeXStringSet(prot_hits,str_c(hits_prot_path,"/",marker,".fasta"))
  }

  # Visualise histogram of alignment
  # Create histograms for each alignment quality metric
  p1 <- ggplot(tab_best, aes(x = ident)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "white", alpha = 0.8) +
    labs(title = "Distribution of identity percentage (fident)",
        x = "Identity percentage",
        y = "Number of alignments") +
    theme_minimal()

  p2 <- ggplot(tab_best, aes(x = alnlen)) +
    geom_histogram(bins = 30, fill = "darkgreen", color = "white", alpha = 0.8) +
    labs(title = "Distribution of alignment lengths (alnlen)",
        x = "Alignment length (aa)",
        y = "Number of alignments") +
    theme_minimal()

  p3 <- ggplot(tab_best, aes(x = evalue)) +
    geom_histogram(bins = 50, fill = "orange", color = "white", alpha = 0.8) +
    scale_x_log10() +
    labs(title = "Distribution of E-values (evalue)",
        x = "E-value (log10 scale)",
        y = "Number of alignments") +
    theme_minimal()

  p4 <- ggplot(tab_best, aes(x = bits)) +
    geom_histogram(bins = 30, fill = "purple", color = "white", alpha = 0.8) +
    labs(title = "Distribution of bit scores (bits)",
        x = "Bit score",
        y = "Number of alignments") +
    theme_minimal()

  # Arrange and display plots
  # Combine plots
  combined_plot <- (p1 | p2) / (p3 | p4)
  
  # Save plot
  ggsave(
    filename = str_c(out_dir, "/", marker, ".svg"),
    plot = combined_plot,
    width = 12,
    height = 8,
    units = "in",
    dpi = 300
  )

  # Analyse part of the database witch is not covered by the pool
  marker_orf_col <- str_c(marker, "_orf_names")

  # Family of the database ref with have the marker 
  meta2 <- meta |>
    filter(!is.na(.data[[marker_orf_col]]))
  meta2$Family |> unique()


  # Analyse on many query has best match with wich family 
  meta_marker_long <- meta %>%
    select(virus_id, Family, all_of(marker_orf_col)) %>%
    filter(!is.na(.data[[marker_orf_col]])) %>%
    separate_rows(all_of(marker_orf_col), sep = ";") %>%
    rename(ref = all_of(marker_orf_col))

  family_count <- tab_best |>
    left_join(meta_marker_long, by = "ref") |>
    group_by(Family) |>
    summarise(
      nb_queries = n_distinct(query),
      .groups = "drop"
    ) |>
    arrange(desc(nb_queries))
  print(family_count)


}



# Launch
markers <- sub("\\.m8$", "", list.files(hits_dir, pattern = "\\.m8$", full.names = FALSE))

lapply(markers, function(marker) analyse_res(marker,ncbi_dir,hits_dir,out_dir,missing_prot_path,hits_prot_path,meta))
