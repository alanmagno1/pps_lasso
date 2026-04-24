#!/usr/bin/env Rscript

cran_packages <- c(
  "caret",
  "data.table",
  "doParallel",
  "fst",
  "getopt",
  "glmnet",
  "readr",
  "stringr"
)
bioc_packages <- c("impute")

install_missing <- function(packages, installer) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    installer(missing)
  }
}

install_missing(
  cran_packages,
  function(pkgs) install.packages(pkgs, repos = "https://cloud.r-project.org")
)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
install_missing(bioc_packages, BiocManager::install)

message("All required packages are installed.")
