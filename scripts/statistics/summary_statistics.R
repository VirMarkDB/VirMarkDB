library(tidyverse)

database <- commandArgs(trailingOnly = TRUE)[[1]]
dir.create("analysis/summary", recursive = TRUE)
dir.create("analysis/figures", recursive = TRUE)
manifest <- read_tsv("manifest_genomes.tsv", show_col_types = FALSE)
markers <- read_tsv("marker_taxo_map.tsv", show_col_types = FALSE) %>%
  distinct(marker_group_id, group_id, marker)
associations <- read_tsv("genome_marker_map.tsv", show_col_types = FALSE)
composition <- read_tsv(file.path(database, "virus_informations/virus_compo_taxo.tsv"), show_col_types = FALSE)
count_taxa <- function(x) n_distinct(x[!is.na(x) & x != ""])
write_summary <- function(data, name) write_tsv(data, file.path("analysis/summary", paste0(name, ".tsv")))
save_figure <- function(name, plot, width = 9, height = 6) {
  for (ext in c("pdf", "png")) {
    ggsave(file.path("analysis/figures", paste0(name, ".", ext)), plot,
           width = width, height = height, dpi = 300)
  }
}

# The genome-marker map defines the denominator, including genomes with no retained marker.
scope <- associations %>% left_join(markers, by = "marker_group_id") %>%
  distinct(virus_id, group_id) %>% left_join(manifest, by = "virus_id")
scope_summary <- scope %>% group_by(group_id) %>%
  summarise(genomes = n_distinct(virus_id), families = count_taxa(Family),
            genera = count_taxa(Genus), species = count_taxa(Species), .groups = "drop")
write_summary(scope_summary, "scope_summary")
write_summary(filter(scope, virus_id %in% composition$virus_id) %>%
  group_by(group_id) %>% summarise(genomes = n_distinct(virus_id), families = count_taxa(Family),
                                  .groups = "drop"), "database_summary")
write_summary(tibble(metric = c("Viral groups", "Genomes", "Families", "Genera", "Species", "Marker genes"),
  value = c(n_distinct(scope$group_id), n_distinct(scope$virus_id), count_taxa(scope$Family),
            count_taxa(scope$Genus), count_taxa(scope$Species), nrow(markers))), "release_summary")

summaries <- list()
lengths <- tibble(marker_group_id = character(), group_id = character(), protein_length = integer())
scores <- tibble(marker_group_id = character(), group_id = character(), orf_name = character(), score = double(), evalue = double())
coverage <- tibble(marker_group_id = character(), group_id = character(), Family = character(), genomes = integer(), detected = integer(), recovery = double())
diversity <- tibble(marker_group_id = character(), rank = character(), n_taxa = integer())
for (i in seq_len(nrow(markers))) {
  m <- markers[i, ]
  id <- m$marker_group_id
  orf_col <- paste0(id, "_orf_names")
  eligible <- associations %>% filter(marker_group_id == id) %>%
    select(virus_id) %>% left_join(manifest, by = "virus_id")
  # Missing columns and FASTAs mean zero retained hits, not a failed statistics task.
  counts <- if (orf_col %in% names(composition)) {
    composition %>% transmute(virus_id, copies = if_else(is.na(.data[[orf_col]]), 0L,
                                                      str_count(.data[[orf_col]], ";") + 1L))
  } else tibble(virus_id = character(), copies = integer())
  eligible <- eligible %>% left_join(counts, by = "virus_id") %>% mutate(copies = replace_na(copies, 0L))
  marker_dir <- file.path(database, "markers", m$group_id, m$marker)
  fasta <- file.path(marker_dir, paste0(id, "_protein.fasta"))
  description <- file.path(marker_dir, paste0(id, "_virus_orf_description.tsv"))
  aa <- if (file.exists(fasta)) Biostrings::readAAStringSet(fasta) else Biostrings::AAStringSet()
  desc <- if (file.exists(description)) read_tsv(description, show_col_types = FALSE) else
    tibble(orf_name = character(), score = double(), evalue = double())
  summaries[[i]] <- tibble(marker_group_id = id, marker = m$marker, group_id = m$group_id,
    genomes_in_scope = nrow(eligible), genomes_with_marker = sum(eligible$copies > 0),
    recovery = if (nrow(eligible)) 100 * mean(eligible$copies > 0) else 0,
    sequences = sum(eligible$copies), sequences_fasta = length(aa), multicopy_genomes = sum(eligible$copies > 1),
    median_length_aa = if (length(aa)) median(Biostrings::width(aa)) else NA_real_,
    median_hmm_score = if (nrow(desc)) median(desc$score) else NA_real_)
  lengths <- bind_rows(lengths, tibble(marker_group_id = rep(id, length(aa)),
    group_id = rep(m$group_id, length(aa)), protein_length = Biostrings::width(aa)))
  scores <- bind_rows(scores, desc %>% transmute(marker_group_id = id, group_id = m$group_id, orf_name, score, evalue))
  coverage <- bind_rows(coverage, eligible %>% filter(!is.na(Family)) %>% group_by(Family) %>%
    summarise(genomes = n(), detected = sum(copies > 0), recovery = 100 * detected / genomes, .groups = "drop") %>%
    mutate(marker_group_id = id, group_id = m$group_id))
  positives <- filter(eligible, copies > 0)
  diversity <- bind_rows(diversity, tibble(marker_group_id = id, rank = c("Family", "Genus", "Species"),
    n_taxa = c(count_taxa(positives$Family), count_taxa(positives$Genus), count_taxa(positives$Species))))
  if (sum(eligible$copies) != length(aa) || length(aa) != nrow(desc)) stop("Inconsistent sequence counts for ", id)
}
summary <- bind_rows(summaries)
write_summary(summary, "marker_summary")
write_summary(scores, "marker_hmm_scores")
write_summary(diversity, "marker_taxonomic_diversity")
write_summary(coverage, "family_marker_coverage")

if (nrow(summary)) {
  p <- ggplot(summary, aes(marker_group_id, recovery, fill = group_id)) + geom_col() +
    coord_flip() + labs(x = NULL, y = "Genomes with marker (%)") + theme_classic()
  save_figure("Figure2A_marker_recovery", p)
  p <- summary %>% select(marker_group_id, genomes_with_marker, sequences) %>%
    pivot_longer(-marker_group_id) %>% ggplot(aes(marker_group_id, value, fill = name)) +
    geom_col(position = "dodge") + coord_flip() + theme_classic()
  save_figure("Figure2B_marker_positive_vs_sequences", p)
  p <- ggplot(summary, aes(marker_group_id, multicopy_genomes, fill = group_id)) +
    geom_col() + coord_flip() + theme_classic()
  save_figure("Figure2D_multicopy_genomes", p)
  p <- ggplot(diversity, aes(marker_group_id, n_taxa, fill = rank)) +
    geom_col(position = "dodge") + coord_flip() + theme_classic()
  save_figure("Figure2E_marker_taxonomic_diversity", p)
}
if (nrow(lengths)) {
  p <- ggplot(lengths, aes(marker_group_id, protein_length, fill = group_id)) +
    geom_boxplot() + coord_flip() + theme_classic()
  save_figure("Figure2C_protein_length_distribution", p)
}
if (nrow(scores)) {
  p <- ggplot(scores, aes(score, marker_group_id, fill = group_id)) + geom_boxplot() + theme_classic()
  save_figure("FigureS1_HMM_scores", p)
}
if (nrow(coverage)) {
  # One generic heatmap supports every group present in the mapfile.
  p <- ggplot(coverage, aes(marker_group_id, Family, fill = recovery)) + geom_tile() +
    facet_wrap(~group_id, scales = "free") + scale_fill_viridis_c(limits = c(0, 100)) +
    theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_figure("family_marker_coverage", p, height = max(6, n_distinct(coverage$Family) * 0.25))
}
