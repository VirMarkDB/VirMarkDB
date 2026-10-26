summary_dir <- "data/analysis/summary"

scope <- read.delim(file.path(summary_dir, "scope_summary.tsv"))

scope_stat <- function(name, column) {
  format(scope[scope$scope == name, column], big.mark = ",")
}

database <- read.delim(file.path(summary_dir, "database_summary.tsv"))

database_stat <- function(name, column) {
  format(database[database$scope == name, column], big.mark = ",")
}