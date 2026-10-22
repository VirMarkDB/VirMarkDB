#!/usr/bin/env Rscript

library(tidyverse)

# ============================================================
# Paths
# ============================================================

vmdb_path <- "output/VirMarkDB"
out_dir <- "docs/data/analysis"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dirs <- c(
  "summary",
  "nucleocytoviricota",
  "pantevenvirales",
  "orthornavirae",
  "figures"
)

for(d in dirs){
  dir.create(
    file.path(out_dir,d),
    recursive=TRUE,
    showWarnings=FALSE
  )
}


# ============================================================
# Load metadata
# ============================================================

virus <- read_tsv(
  file.path(
    vmdb_path,
    "virus_informations/virus_compo_taxo.tsv"
  ),
  show_col_types = FALSE
)


# ============================================================
# Marker definition
# ============================================================

ncldv_markers <- c(
  "atpase_nucleocytoviricota_orf_names",
  "dnapol_nucleocytoviricota_orf_names",
  "majorcapsid_nucleocytoviricota_orf_names",
  "primase_nucleocytoviricota_orf_names",
  "rnapol1_nucleocytoviricota_orf_names",
  "rnapol2_nucleocytoviricota_orf_names",
  "tf2s_nucleocytoviricota_orf_names",
  "vltf3_nucleocytoviricota_orf_names"
)

rdrp_marker <- "rdrp_orthornavirae_orf_names"

g23_marker <- "g23_pantevenvirales_orf_names"


marker_labels <- c(
  atpase_nucleocytoviricota_orf_names="Packaging ATPase",
  dnapol_nucleocytoviricota_orf_names="DNA polymerase",
  majorcapsid_nucleocytoviricota_orf_names="Major capsid protein",
  primase_nucleocytoviricota_orf_names="Primase",
  rnapol1_nucleocytoviricota_orf_names="RNA polymerase subunit 1",
  rnapol2_nucleocytoviricota_orf_names="RNA polymerase subunit 2",
  tf2s_nucleocytoviricota_orf_names="TFIIS",
  vltf3_nucleocytoviricota_orf_names="VLTF3",
  rdrp_orthornavirae_orf_names="RdRp",
  g23_pantevenvirales_orf_names="g23"
)


count_marker <- function(x){

  ifelse(
    is.na(x) | x=="",
    0,
    str_count(x,";")+1
  )

}


has_marker <- function(x){

  !is.na(x) & x!=""

}


# ============================================================
# VirMarkDB scope
# ============================================================

virus <- virus %>%
  mutate(
    virmark_scope = case_when(
      Phylum=="Nucleocytoviricota" ~ "Nucleocytoviricota",
      Kingdom=="Orthornavirae" ~ "Orthornavirae",
      Order=="Pantevenvirales" ~ "Pantevenvirales",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(virmark_scope))


# ============================================================
# Global summary
# ============================================================

release_summary <- tibble(
  metric=c(
    "Viral groups",
    "Genomes",
    "Families",
    "Genera",
    "Species",
    "Marker genes"
  ),
  value=c(
    n_distinct(virus$virmark_scope),
    n_distinct(virus$virus_id),
    n_distinct(virus$Family),
    n_distinct(virus$Genus),
    n_distinct(virus$Species),
    10
  )
)


write_tsv(
 release_summary,
 file.path(out_dir,"summary/release_summary.tsv")
)


# ============================================================
# Genome distribution among groups
# ============================================================

scope_summary <- virus %>%
 group_by(virmark_scope)%>%
 summarise(
   genomes=n_distinct(virus_id),
   families=n_distinct(Family),
   species=n_distinct(Species),
   .groups="drop"
 )

write_tsv(
 scope_summary,
 file.path(out_dir,"summary/scope_summary.tsv")
)


p <- ggplot(
 scope_summary,
 aes(
  reorder(virmark_scope,genomes),
  genomes
 )
)+
geom_col()+
coord_flip()+
labs(
 x=NULL,
 y="Number of genomes"
 )

ggsave(
 file.path(out_dir,"figures/release_groups.png"),
 p,
 width=6,
 height=4
)


# ============================================================
# NUCLEOCYTOVIRICOTA
# ============================================================

ncldv <- virus %>%
 filter(virmark_scope=="Nucleocytoviricota")


# Marker summary

ncldv_marker_summary <- map_dfr(
 ncldv_markers,
 function(m){

 tibble(
 marker=marker_labels[m],
 genomes=n_distinct(
  ncldv$virus_id[
   has_marker(ncldv[[m]])
  ]
 ),
 sequences=sum(count_marker(ncldv[[m]]))
 )

})


write_tsv(
 ncldv_marker_summary,
 file.path(
 out_dir,
 "nucleocytoviricota/marker_summary.tsv"
 )
)


p <- ggplot(
 ncldv_marker_summary,
 aes(
  reorder(marker,genomes),
  genomes
 )
)+
geom_col()+
coord_flip()+
labs(
 x=NULL,
 y="Genomes containing marker"
 )

ggsave(
 file.path(out_dir,"figures/ncldv_marker_coverage.png"),
 p,
 width=7,
 height=5
)



# Family marker coverage

family_marker <- ncldv %>%
 group_by(Family)%>%
 summarise(
 across(
  all_of(ncldv_markers),
  ~sum(has_marker(.x))
 ),
 genomes=n_distinct(virus_id),
 .groups="drop"
 )


write_tsv(
 family_marker,
 file.path(
 out_dir,
 "nucleocytoviricota/family_marker_coverage.tsv"
 )
)


heat <- family_marker %>%
 pivot_longer(
  cols = all_of(ncldv_markers),
  names_to = "marker",
  values_to = "marker_count"
 ) %>%
 mutate(
  marker = marker_labels[marker]
 )


p <- ggplot(heat, aes( marker, Family, fill=marker_count))+
geom_tile()+
theme( axis.text.x=element_text(  angle=45,  hjust=1 ))+
labs( x=NULL, y=NULL )

ggsave(
 file.path(out_dir,"figures/ncldv_family_heatmap.png"),
 p,
 width=9,
 height=7
)



# Marker completeness

ncldv_complete <- ncldv %>%
 mutate(
 marker_number=rowSums(
  across(
   all_of(ncldv_markers),
   has_marker
  )
 )
)


write_tsv(
 ncldv_complete %>% count(marker_number),
 file.path(
 out_dir,
 "nucleocytoviricota/genome_marker_completeness.tsv"
 )
)


p <- ggplot(
 ncldv_complete,
 aes(marker_number)
)+
geom_histogram(
 bins=10
)+
labs(
 x="Number of detected markers",
 y="Number of genomes"
 )

ggsave(
 file.path(out_dir,"figures/ncldv_marker_completeness.png"),
 p,
 width=6,
 height=4
)



# ============================================================
# PANTEVENVIRALES
# ============================================================

panteven <- virus %>%
 filter(virmark_scope=="Pantevenvirales")


g23_family <- panteven %>%
 group_by(Family)%>%
 summarise(
 genomes=n_distinct(virus_id),
 g23_detected=sum(
  has_marker(.data[[g23_marker]])
 ),
 .groups="drop"
 )


write_tsv(
 g23_family,
 file.path(
 out_dir,
 "pantevenvirales/family_g23_coverage.tsv"
 )
)


p <- ggplot(
 g23_family,
 aes(
 reorder(Family,g23_detected),
 g23_detected
 )
)+
geom_col()+
coord_flip()+
labs(
 x=NULL,
 y="Genomes with g23"
 )

ggsave(
 file.path(
 out_dir,
 "figures/pantevenvirales_family_g23.png"
 ),
 p,
 width=7,
 height=5
)



# ============================================================
# ORTHORNAVIRAE
# ============================================================

rna <- virus %>%
 filter(virmark_scope=="Orthornavirae")


phylum_rdrp <- rna %>%
 group_by(Phylum)%>%
 summarise(
 genomes=n_distinct(virus_id),
 families=n_distinct(Family),
 rdrp=sum(
  has_marker(.data[[rdrp_marker]])
 ),
 .groups="drop"
 )


write_tsv(
 phylum_rdrp,
 file.path(
 out_dir,
 "orthornavirae/phylum_rdrp_coverage.tsv"
 )
)


p <- ggplot(
 phylum_rdrp,
 aes(
 reorder(Phylum,genomes),
 genomes
 )
)+
geom_col()+
coord_flip()+
labs(
 x=NULL,
 y="Number of genomes"
 )

ggsave(
 file.path(
 out_dir,
 "figures/orthornavirae_phylum.png"
 ),
 p,
 width=7,
 height=5
)



family_rdrp <- rna %>%
 group_by(Phylum,Family)%>%
 summarise(
 genomes=n_distinct(virus_id),
 rdrp=sum(
  has_marker(.data[[rdrp_marker]])
 ),
 .groups="drop"
 )


write_tsv(
 family_rdrp,
 file.path(
 out_dir,
 "orthornavirae/family_rdrp_coverage.tsv"
 )
)


writeLines(
"VirMarkDB analysis statistics generated successfully."
)

