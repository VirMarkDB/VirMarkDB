summary_dir <- "data/analysis/summary"

release <- read.delim(file.path(summary_dir, "release_summary.tsv"))
scope <- read.delim(file.path(summary_dir, "scope_summary.tsv"))
database <- read.delim(file.path(summary_dir, "database_summary.tsv"))

stat <- function(name) {
  x <- switch(name,
    "Viral groups" = nrow(database),
    "Marker genes" = as.numeric(release$value[release$metric == "Marker genes"]),
    "Genomes" = sum(database$genomes)
  )
  format(x, big.mark = ",")
}

scope_stat <- function(name, column) {
  format(scope[scope$scope == name, column], big.mark = ",")
}

database_stat <- function(name, column) {
  format(database[database$scope == name, column], big.mark = ",")
}