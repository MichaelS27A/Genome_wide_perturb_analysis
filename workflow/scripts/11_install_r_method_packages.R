#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option("--marker", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$marker)) stop("--marker is required")

dir.create(dirname(opt$marker), recursive = TRUE, showWarnings = FALSE)

cran_repo <- "https://cloud.r-project.org"

ensure_cran <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = cran_repo)
  }
}

ensure_bioc <- function(pkg) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = cran_repo)
  }
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

ensure_github <- function(pkg, repo) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools", repos = cran_repo)
  }
  if (!requireNamespace(pkg, quietly = TRUE)) {
    devtools::install_github(repo, upgrade = "never")
  }
}

# Base dependencies used by Mixscale/PS workflows.
for (p in c("Seurat", "PMA", "protoclust", "jsonlite", "Matrix", "reticulate", "devtools")) {
  ensure_cran(p)
}

ensure_bioc("glmGamPoi")

# Method packages from official repositories.
ensure_github("Mixscale", "longmanz/Mixscale")
ensure_github("scMAGeCK", "weililab/scMAGeCK")

writeLines("ok", con = opt$marker)
cat("[ok] installed/verified Mixscale + scMAGeCK and dependencies\n")
