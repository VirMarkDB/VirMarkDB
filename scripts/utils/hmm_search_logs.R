library(tidyverse)

dir.create(OUT_LOGS, recursive = TRUE, showWarnings = FALSE)

selected_orfs_log <- selected_orfs 

# ORFs shared across marker groups

orf_shared_log <- selected_orfs_log %>%
  group_by(orf_name) %>%
  summarise(
    n_marker_groups = n_distinct(marker_group_id),
    marker_group_ids = str_c(sort(unique(marker_group_id)), collapse = ";"),
    .groups = "drop"
  ) %>%
  filter(n_marker_groups > 1)

write_tsv(
  orf_shared_log,
  file.path(OUT_LOGS, "log_orf_shared_across_markers.tsv")
)


# ORF length plot

orf_length_raw <- selected_orfs_log %>%
  select(
    virus_id,
    marker_group_id,
    group_id,
    marker,
    orf_name,
    target_length_aa
  )

orf_length_summary <- orf_length_raw %>%
  group_by(group_id, marker) %>%
  summarise(
    n_orfs = n(),
    min_length = min(target_length_aa, na.rm = TRUE),
    median_length = median(target_length_aa, na.rm = TRUE),
    mean_length = mean(target_length_aa, na.rm = TRUE),
    max_length = max(target_length_aa, na.rm = TRUE),
    .groups = "drop"
  )

p_orf_length <- orf_length_raw %>%
  ggplot(aes(x = marker_group_id, y = target_length_aa)) +
  geom_boxplot(outlier.size = 0.5) +
  coord_flip() +
  facet_wrap(~marker, scales = "free_y") +
  labs(
    x = "Marker group",
    y = "ORF length",
    title = "ORF length distributions"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 8),
    legend.position = "none"
  )

ggsave(
  file.path(OUT_LOGS, "plot_orf_length_boxplot.png"),
  p_orf_length,
  width = 12,
  height = 8
)


# Marker summary for README

marker_summary <- selected_orfs_log %>%
  distinct(Family, marker, virus_id) %>%
  count(Family, marker, name = "n_virus_genomes") %>%
  pivot_wider(
    names_from = marker,
    values_from = n_virus_genomes,
    values_fill = 0
  )%>%
  arrange(Family) %>%
  bind_rows(
    summarise(., across(where(is.numeric), sum), Family = "Total genomes")
  )

marker_summary_md <- knitr::kable(
  marker_summary,
  format = "pipe"
)

write_lines(
  marker_summary_md,
  "output/hmm/logs/family_marker_summary_wide.md"
)

# Family x marker barplot
marker_order <- c(  "atpase",  "dnapol",  "majorcapsid",  "primase",  "rnapol1",  "rnapol2",  "tf2s",  "vltf3","rdrp", "g23")
marker_colors <- c(atpase= "#4E79A7",dnapol= "#F28E2B",majorcapsid= "#59A14F",primase= "#E15759",rnapol1= "#B07AA1",rnapol2= "#9C755F",tf2s= "#76B7B2",vltf3= "#EDC948")

p_family_marker_summary  <- selected_orfs_log %>%
  distinct(Family, marker, virus_id) %>%
  count(Family, marker, name = "n_virus_genomes") %>%
  mutate(
    Family = fct_reorder(Family, n_virus_genomes, .fun = sum),
    marker = factor(marker, levels = marker_order)
  ) %>%
  ggplot(aes(x = Family, y = n_virus_genomes, fill = marker)) +
  geom_col(position = "dodge") +
  facet_wrap(~marker, scales = "free_y") +
  scale_fill_manual(values = marker_colors, drop = FALSE) +
  coord_flip() +
  labs(
    x = "Family",
    y = "Virus genomes",
    fill = "Marker",
    title = "Number of virus genomes by family and marker"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.y = element_text(size = 8),
    legend.position = "right"
  )

ggsave(
  file.path(OUT_LOGS, "plot_family_marker_summary_barplot.png"),
  p_family_marker_summary,
  width = 12,
  height = 8
)