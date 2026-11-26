library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
mode <- args[[1]]
marker_id <- args[[2]]
proteins <- Biostrings::readAAStringSet(args[[3]])
ids <- sub(" .*", "", names(proteins))
columns <- c("query", "ref", "ident", "alnlen", "mismatch", "gapopen",
             "qstart", "qend", "sstart", "send", "evalue", "bits")
# Empty MMseqs tables are valid when there are no matching proteins.
hits <- if (file.info(args[[4]])$size > 0) {
  read_tsv(args[[4]], col_names = columns, show_col_types = FALSE)
} else {
  tibble(query = character(), ref = character(), ident = double(), alnlen = integer(),
         mismatch = integer(), gapopen = integer(), qstart = integer(), qend = integer(),
         sstart = integer(), send = integer(), evalue = double(), bits = double())
}

if (mode == "filter") {
  # Identity is explicitly requested as a fraction from MMseqs, matching the 0.3 threshold.
  retained <- hits %>%
    filter(query != ref, ident >= 0.3, alnlen >= 100, evalue <= 1e-5) %>%
    pull(query) %>% unique()
  dir.create("filtered")
  selected <- proteins[!duplicated(ids) & ids %in% retained]
  Biostrings::writeXStringSet(selected, file.path("filtered", paste0(marker_id, ".fasta")))
  write_tsv(tibble(marker_group_id = marker_id, raw = length(proteins), retained = length(selected)),
            file.path("filtered", paste0(marker_id, "_summary.tsv")))
} else {
  dir.create("results")
  dir.create("missing")
  dir.create("hits")
  best <- hits %>% group_by(query) %>% slice_min(evalue, with_ties = FALSE) %>% ungroup()
  write_tsv(best, file.path("results", paste0(marker_id, "_best_hits.tsv")))
  Biostrings::writeXStringSet(proteins[!ids %in% best$query], file.path("missing", paste0(marker_id, ".fasta")))
  Biostrings::writeXStringSet(proteins[ids %in% best$query], file.path("hits", paste0(marker_id, ".fasta")))
  meta <- read_tsv(file.path(args[[5]], "virus_informations/virus_compo_taxo.tsv"), show_col_types = FALSE)
  orf_col <- paste0(marker_id, "_orf_names")
  families <- meta %>% select(Family, all_of(orf_col)) %>%
    separate_rows(all_of(orf_col), sep = ";") %>% rename(ref = all_of(orf_col))
  family_counts <- best %>% left_join(families, by = "ref") %>% count(Family, name = "nb_queries")
  write_tsv(family_counts, file.path("results", paste0(marker_id, "_family_counts.tsv")))
  write_tsv(tibble(marker_group_id = marker_id, queries = length(unique(ids)),
                   matched = nrow(best), missing = length(setdiff(unique(ids), best$query))),
            file.path("results", paste0(marker_id, "_summary.tsv")))
}

# The same quality panels are useful before filtering and after database matching.
plot_hits <- if (mode == "filter") filter(hits, query != ref) else best
if (nrow(plot_hits)) {
  metrics <- plot_hits %>% select(ident, alnlen, evalue, bits) %>%
    pivot_longer(everything(), names_to = "metric", values_to = "value")
  p <- ggplot(metrics, aes(value)) + geom_histogram(bins = 30) +
    facet_wrap(~metric, scales = "free") + theme_minimal()
  ggsave(file.path(if (mode == "filter") "filtered" else "results", paste0(marker_id, ".pdf")),
         p, width = 12, height = 8)
}
