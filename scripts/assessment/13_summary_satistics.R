#!/usr/bin/env Rscript

library(tidyverse)

# Paths
vmdb_path <- "output/VirMarkDB"
manifest_file <- "output/config/manifest_genomes.tsv"
out_dir <- "website/data/analysis"
fig_dir <- file.path(out_dir, "figures")
summary_dir <- file.path(out_dir, "summary")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

# Colors
colors <- c(Nucleocytoviricota = "#4778A8", Pantevenvirales = "#459475", Orthornavirae = "#C47C55")
theme_article <- theme_classic(base_size = 11) +
  theme(axis.text = element_text(color = "#333333"), legend.position = "top",
        legend.title = element_blank(), plot.title = element_text(face = "bold"))

save_figure <- function(name, plot, width = 8, height = 5){
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), plot, width = width, height = height)
  ggsave(file.path(fig_dir, paste0(name, ".png")), plot, width = width, height = height, dpi = 400)
}

count_marker <- function(x) ifelse(is.na(x) | trimws(x) == "", 0L, str_count(x, ";") + 1L)
count_taxa <- function(x) n_distinct(x[!is.na(x) & x != ""])

# Data
manifest <- read_tsv(manifest_file, show_col_types = FALSE) %>%
  distinct(virus_id, .keep_all = TRUE) %>%
  mutate(scope = case_when(Phylum == "Nucleocytoviricota" ~ "Nucleocytoviricota",
                           Order == "Pantevenvirales" ~ "Pantevenvirales",
                           Kingdom == "Orthornavirae" ~ "Orthornavirae")) %>%
  filter(!is.na(scope))

virus <- read_tsv(file.path(vmdb_path, "virus_informations/virus_compo_taxo.tsv"), show_col_types = FALSE)

markers <- tibble(
  marker_id = c("atpase_nucleocytoviricota", "dnapol_nucleocytoviricota", "majorcapsid_nucleocytoviricota",
                "primase_nucleocytoviricota", "rnapol1_nucleocytoviricota", "rnapol2_nucleocytoviricota",
                "tf2s_nucleocytoviricota", "vltf3_nucleocytoviricota", "g23_pantevenvirales", "rdrp_orthornavirae"),
  marker = c("Packaging ATPase", "DNA polymerase", "Major capsid protein", "Primase",
             "RNA polymerase subunit 1", "RNA polymerase subunit 2", "TFIIS", "VLTF3", "gp23", "RdRp"),
  scope = c(rep("Nucleocytoviricota", 8), "Pantevenvirales", "Orthornavirae")
) %>% mutate(orf_col = paste0(marker_id, "_orf_names"))

missing_cols <- setdiff(markers$orf_col, names(virus))
if(length(missing_cols) > 0) stop("Missing marker columns: ", paste(missing_cols, collapse = ", "))
if(anyDuplicated(virus$virus_id)) stop("Duplicated virus_id in virus_compo_taxo.tsv")
if(any(!virus$virus_id %in% manifest$virus_id)) stop("Virus IDs not found in manifest")

data <- manifest %>% left_join(virus %>% select(virus_id, all_of(markers$orf_col)), by = "virus_id")

database_summary <- data %>%
  filter(if_any(ends_with("_orf_names"), ~ !is.na(.x) & .x != "")) %>%
  group_by(scope) %>%
  summarise(
    genomes = n_distinct(virus_id),
    families = n_distinct(Family[!is.na(Family) & Family != ""]),
    .groups = "drop"
  )

write_tsv(database_summary, file.path(summary_dir, "database_summary.tsv"))

# Global summary
scope_summary <- manifest %>% group_by(scope) %>% summarise(
  genomes = n_distinct(virus_id), families = count_taxa(Family), genera = count_taxa(Genus),
  species = count_taxa(Species), .groups = "drop")
write_tsv(scope_summary, file.path(summary_dir, "scope_summary.tsv"))

release_summary <- tibble(metric = c("Viral groups", "Genomes", "Families", "Genera", "Species", "Marker genes"),
                          value = c(nrow(scope_summary), n_distinct(manifest$virus_id), count_taxa(manifest$Family),
                                    count_taxa(manifest$Genus), count_taxa(manifest$Species), nrow(markers)))
write_tsv(release_summary, file.path(summary_dir, "release_summary.tsv"))

# Marker statistics
marker_summary <- tibble()
length_data <- tibble()
score_data <- tibble()
taxo_data <- tibble()

for(i in seq_len(nrow(markers))){
  id <- markers$marker_id[i]
  marker <- markers$marker[i]
  group <- markers$scope[i]
  col <- markers$orf_col[i]
  sub <- data %>% filter(scope == group)
  copies <- count_marker(sub[[col]])

  fasta <- list.files(file.path(vmdb_path, "markers"), pattern = paste0("^", id, "_protein\\.fasta$"),
                      recursive = TRUE, full.names = TRUE)
  desc_file <- list.files(file.path(vmdb_path, "markers"), pattern = paste0("^", id, "_virus_orf_description\\.tsv$"),
                          recursive = TRUE, full.names = TRUE)
  if(length(fasta) != 1 || length(desc_file) != 1) stop("Missing or duplicated marker files: ", id)

  aa <- Biostrings::readAAStringSet(fasta)
  desc <- read_tsv(desc_file, show_col_types = FALSE)
  if(!all(c("orf_name", "score", "evalue") %in% names(desc))) stop("Missing columns in: ", desc_file)

  marker_summary <- bind_rows(marker_summary, tibble(
    marker_id = id, marker = marker, scope = group,
    genomes_in_scope = nrow(sub), genomes_with_marker = sum(copies > 0),
    recovery = 100 * mean(copies > 0), sequences = sum(copies),
    sequences_fasta = length(aa), multicopy_genomes = sum(copies > 1),
    median_length_aa = median(BiocGenerics::width(aa)), median_hmm_score = median(desc$score, na.rm = TRUE)))

  length_data <- bind_rows(length_data, tibble(marker = marker, scope = group, protein_length = BiocGenerics::width(aa)))
  score_data <- bind_rows(score_data, desc %>% transmute(marker = .env$marker, scope = .env$group, orf_name,
                                                        score = as.numeric(score), evalue = as.numeric(evalue)))

  positives <- sub[copies > 0, ]
  taxo_data <- bind_rows(taxo_data, tibble(marker = marker, scope = group,
    rank = c("Family", "Genus", "Species"),
    n_taxa = c(count_taxa(positives$Family), count_taxa(positives$Genus), count_taxa(positives$Species))))

  if(sum(copies) != length(aa) || length(aa) != nrow(desc))
    warning("Inconsistent sequence counts for ", id, ": table=", sum(copies),
            ", fasta=", length(aa), ", description=", nrow(desc))
}

marker_summary$marker <- factor(marker_summary$marker, levels = markers$marker)
length_data$marker <- factor(length_data$marker, levels = markers$marker)
score_data$marker <- factor(score_data$marker, levels = markers$marker)
taxo_data$marker <- factor(taxo_data$marker, levels = markers$marker)

write_tsv(marker_summary, file.path(summary_dir, "marker_summary.tsv"))
write_tsv(taxo_data, file.path(summary_dir, "marker_taxonomic_diversity.tsv"))
if (anyNA(score_data$marker)) stop("Missing marker names in HMM score data")
write_tsv(score_data, file.path(summary_dir, "marker_hmm_scores.tsv"))

# Figure 2A - Marker recovery
p <- ggplot(marker_summary, aes(reorder(marker, recovery), recovery, fill = scope)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", recovery)), hjust = -0.15, size = 3) +
  coord_flip() + scale_fill_manual(values = colors) +
  scale_y_continuous(limits = c(0, 112), breaks = seq(0, 100, 20)) +
  labs(x = NULL, y = "Genomes with marker (%)") + theme_article
save_figure("Figure2A_marker_recovery", p)

# Figure 2B - Genomes and retained sequences
plot_data <- marker_summary %>% select(marker, genomes_with_marker, sequences) %>%
  pivot_longer(-marker, names_to = "type", values_to = "n") %>%
  mutate(type = recode(type, genomes_with_marker = "Marker-positive genomes", sequences = "Retained sequences"))
p <- ggplot(plot_data, aes(marker, n, fill = type)) +
  geom_col(position = "dodge", width = 0.7) + coord_flip() +
  scale_fill_manual(values = c("Marker-positive genomes" = "#5C85AA", "Retained sequences" = "#CF9168")) +
  labs(x = NULL, y = "Number") + theme_article
save_figure("Figure2B_marker_positive_vs_sequences", p)

# Figure 2C - Protein lengths
p <- ggplot(length_data, aes(marker, protein_length, fill = scope)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.7, color = NA) +
  geom_boxplot(width = 0.13, fill = "white", linewidth = 0.3, outlier.size = 0.3) +
  coord_flip() + scale_fill_manual(values = colors) +
  labs(x = NULL, y = "Protein length (aa)") + theme_article
save_figure("Figure2C_protein_length_distribution", p, height = 5.5)

# Figure 2D - Multicopy genomes
p <- ggplot(marker_summary, aes(reorder(marker, multicopy_genomes), multicopy_genomes, fill = scope)) +
  geom_col(width = 0.7) + coord_flip() + scale_fill_manual(values = colors) +
  labs(x = NULL, y = "Genomes with more than one retained copy") + theme_article
save_figure("Figure2D_multicopy_genomes", p)

# Figure 2E - Taxonomic representation
p <- ggplot(taxo_data, aes(marker, n_taxa, fill = rank)) +
  geom_col(position = "dodge", width = 0.7) + coord_flip() +
  scale_fill_manual(values = c(Family = "#587EA9", Genus = "#55A087", Species = "#CE8A62")) +
  labs(x = NULL, y = "Number of represented taxa") + theme_article
save_figure("Figure2E_marker_taxonomic_diversity", p)

# Figure S1 - HMM bit scores
p <- ggplot(filter(score_data, is.finite(score)), aes(score, marker, fill = scope)) +
  geom_boxplot(width = 0.65, outlier.size = 0.35, linewidth = 0.35) +
  facet_grid(scope ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = colors) +
  labs(x = "Full-sequence HMM bit score", y = NULL) +
  theme_article + theme(legend.position = "none", strip.text.y = element_text(angle = 0))
save_figure("FigureS1_HMM_scores", p, width = 9, height = 6)

# Figure S2 - Nucleocytoviricota marker recovery by family
ncldv_cols <- markers %>% filter(scope == "Nucleocytoviricota") %>% pull(orf_col)
ncldv_heat <- data %>% filter(scope == "Nucleocytoviricota", !is.na(Family), Family != "") %>%
  select(Family, all_of(ncldv_cols)) %>%
  pivot_longer(-Family, names_to = "orf_col", values_to = "orf") %>%
  mutate(marker = markers$marker[match(orf_col, markers$orf_col)], detected = count_marker(orf) > 0) %>%
  group_by(Family, marker) %>%
  summarise(genomes = n(), detected = sum(detected), recovery = 100 * detected / genomes, .groups = "drop")
stopifnot(all(ncldv_heat$recovery >= 0 & ncldv_heat$recovery <= 100))
write_tsv(ncldv_heat, file.path(summary_dir, "ncldv_family_heatmap.tsv"))

family_stats <- ncldv_heat %>% group_by(Family) %>%
  summarise(mean_recovery = mean(recovery), genomes = genomes[1], .groups = "drop") %>%
  arrange(mean_recovery)
family_order <- family_stats$Family
family_labels <- paste0(family_stats$Family, " (n=", family_stats$genomes, ")")
names(family_labels) <- family_stats$Family
ncldv_heat <- ncldv_heat %>%
  mutate(Family = factor(Family, levels = family_order),
         label = if_else(recovery == 0, "0", sprintf("%.0f", recovery)))

p <- ggplot(ncldv_heat, aes(marker, Family, fill = recovery)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = label, color = recovery >= 65), size = 2.7) +
  scale_color_manual(values = c("FALSE" = "#293849", "TRUE" = "white"), guide = "none") +
  scale_fill_gradientn(colors = c("#ECEFF1", "#BFD9E9", "#75A9CC", "#285989"),
                       values = scales::rescale(c(0, 1, 35, 100)), limits = c(0, 100),
                       breaks = c(0, 25, 50, 75, 100)) +
  scale_y_discrete(labels = family_labels) +
  labs(x = NULL, y = NULL, fill = "Recovery (%)") + theme_article +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right",
        legend.title = element_text())
save_figure("FigureS2_ncldv_family_heatmap", p, width = 10, height = max(5, length(family_order) * 0.42))

# Figure S3 - gp23 recovery by family
panteven <- data %>% filter(scope == "Pantevenvirales", !is.na(Family), Family != "") %>%
  group_by(Family) %>% summarise(genomes = n(), detected = sum(count_marker(g23_pantevenvirales_orf_names) > 0),
                                recovery = 100 * detected / genomes, .groups = "drop")
write_tsv(panteven, file.path(summary_dir, "pantevenvirales_family_coverage.tsv"))
p <- ggplot(panteven, aes(reorder(Family, recovery), recovery)) +
  geom_col(fill = colors["Pantevenvirales"], width = 0.7) + coord_flip() +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = NULL, y = "Genomes with gp23 (%)") + theme_article
save_figure("FigureS3_pantevenvirales_gp23", p, width = 7, height = max(4, nrow(panteven) * 0.35))

# Figure S4 - RdRp recovery by phylum
rna <- data %>% filter(scope == "Orthornavirae", !is.na(Phylum), Phylum != "") %>%
  group_by(Phylum) %>% summarise(genomes = n(), detected = sum(count_marker(rdrp_orthornavirae_orf_names) > 0),
                                recovery = 100 * detected / genomes, .groups = "drop")
write_tsv(rna, file.path(summary_dir, "orthornavirae_phylum_coverage.tsv"))
p <- ggplot(rna, aes(reorder(Phylum, recovery), recovery)) +
  geom_col(fill = colors["Orthornavirae"], width = 0.7) + coord_flip() +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = NULL, y = "Genomes with RdRp (%)") + theme_article
save_figure("FigureS4_orthornavirae_rdrp", p, width = 8, height = max(4, nrow(rna) * 0.4))

cat("VirMarkDB statistics and figures generated.\n")