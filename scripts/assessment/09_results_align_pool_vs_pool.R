library(tidyverse)

path_blast <- "output/assessment/pool_protein/auto_alignment/blast"
out <- "output/assessment/pool_protein/auto_alignment/blast/analysis"
pool_clean <- "output/assessment/pool_protein/pool_protein_filt"

if (!dir.exists(out)){dir.create(out)}
if (!dir.exists(pool_clean)){dir.create(pool_clean)}


analyse_res <- function(marker,path_blast,out) {
  # Load data
  aln <- read.delim(str_c(path_blast, "/", marker,".m8"), header = FALSE)
  colnames(aln) <- c("query", "target", "fident", "alnlen", "mismatch", "gapopen", "qstart", "qend", "tstart", "tend", "evalue", "bits")

  # Filter out auto_align
  aln2 <- aln |> filter(query != target)

  # Create histograms for each alignment quality metric
  p1 <- ggplot(aln2, aes(x = fident)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "white", alpha = 0.8) +
    labs(title = "Distribution of identity percentage (fident)",
        x = "Identity percentage",
        y = "Number of alignments") +
    theme_minimal()

  p2 <- ggplot(aln2, aes(x = alnlen)) +
    geom_histogram(bins = 30, fill = "darkgreen", color = "white", alpha = 0.8) +
    labs(title = "Distribution of alignment lengths (alnlen)",
        x = "Alignment length (aa)",
        y = "Number of alignments") +
    theme_minimal()

  p3 <- ggplot(aln2, aes(x = evalue)) +
    geom_histogram(bins = 50, fill = "orange", color = "white", alpha = 0.8) +
    scale_x_log10() +
    labs(title = "Distribution of E-values (evalue)",
        x = "E-value (log10 scale)",
        y = "Number of alignments") +
    theme_minimal()

  p4 <- ggplot(aln2, aes(x = bits)) +
    geom_histogram(bins = 30, fill = "purple", color = "white", alpha = 0.8) +
    labs(title = "Distribution of bit scores (bits)",
        x = "Bit score",
        y = "Number of alignments") +
    theme_minimal()

  # Arrange and display plots
  final_plot <- gridExtra::arrangeGrob(p1, p2, p3, p4, ncol = 2)

  # Save 
  ggsave(
    filename = str_c(out,"/" ,marker,".svg"),
    plot = final_plot,
    width = 12,
    height = 8,
    units = "in",
    dpi = 300)
}

# Filter out prot too far 
filter_prot <- function(marker, path_blast, pool_clean) {

  aln <- read.delim(str_c(path_blast, "/", marker,".m8"), header = FALSE)
  colnames(aln) <- c("query", "target", "fident", "alnlen", "mismatch", "gapopen", "qstart", "qend", "tstart", "tend", "evalue", "bits")

  # Filter out auto_align
  aln2 <- aln |> filter(query != target)

  # make rules
  fident_threshold <- 0.3
  alnlen_threshold <- 100
  evalue_threshold <- 1e-5

  #filter
  good_alignments <- aln2 |>
    filter(
      fident >= fident_threshold,
      alnlen >= alnlen_threshold,
      evalue <= evalue_threshold
    )

  # List prot with at least a good alignment
  queries_with_good_alignments <- unique(good_alignments$query)


  # Print logs
  cat(marker)
  cat("Total number of query protein :", length(unique(aln$query)), "\n")
  cat("Total number of with at least a good blast :", length(queries_with_good_alignments), "\n")
  cat("Total number of isolate prot :", length(unique(aln2$query))- length(queries_with_good_alignments), "\n")

  # filter out bad prot 
  pool <- "output/assessment/pool_protein/pool_protein_raw"

  prot <- Biostrings::readAAStringSet(str_c(pool,"/",marker,".fasta"))

  # Filter AA duplication
  # Initiate dataframe to store results
  prot_dedup <- prot[names(prot)[!duplicated(names(prot))]]
  
  # Filter by prot name
  prot_clean <- prot_dedup[sub(" .*", "", names(prot_dedup)) %in% queries_with_good_alignments]
  
  Biostrings::writeXStringSet(prot_clean, str_c(pool_clean,"/",marker,".fasta"))
}

# Launch
markers <- sub("\\.fasta$", "", list.files("output/assessment/pool_protein/pool_protein_raw", pattern = "\\.fasta$", full.names = FALSE))

lapply(markers, function(marker) analyse_res(marker, path_blast, out))
print("analyse results OK")

lapply(markers, function(marker) filter_prot(marker, path_blast, pool_clean))
print("filter protein OK")