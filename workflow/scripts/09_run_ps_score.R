#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
  library(Matrix)
  library(Seurat)
  library(reticulate)
  library(scMAGeCK)
})

option_list <- list(
  make_option("--h5ad", type = "character"),
  make_option("--outdir", type = "character"),
  make_option("--chunk-cells", type = "character", default = NULL, dest = "chunk_cells"),
  make_option("--chunk-id", type = "character", default = NULL, dest = "chunk_id"),
  make_option("--pert-col", type = "character", default = "gene_target", dest = "pert_col"),
  make_option("--control-label", type = "character", default = "Non-Targeting", dest = "control_label"),
  make_option("--min-cells-per-perturbation", type = "integer", default = 30, dest = "min_cells_per_perturbation"),
  make_option("--max-perturbations", type = "integer", default = 0, dest = "max_perturbations"),
  make_option("--max-cells", type = "integer", default = 0, dest = "max_cells"),
  make_option("--random-seed", type = "integer", default = 0, dest = "random_seed")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$h5ad) || is.null(opt$outdir)) stop("--h5ad and --outdir are required")
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
set.seed(opt$random_seed)

anndata <- import("anndata")
sp <- import("scipy.sparse")
np <- import("numpy")
# Keep AnnData on disk and only materialize per-perturbation subsets.
ad <- anndata$read_h5ad(opt$h5ad, backed = "r")

if (!is.null(opt$chunk_cells) && nzchar(opt$chunk_cells)) {
  chunk_df <- read.table(opt$chunk_cells, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  if (!("cell_barcode" %in% colnames(chunk_df))) {
    stop("chunk-cells file must include 'cell_barcode' column")
  }
  barcodes <- unique(as.character(chunk_df$cell_barcode))
  if (length(barcodes) == 0) stop("No barcodes found in chunk-cells file")
  keep <- ad$obs_names$astype("str")$isin(r_to_py(barcodes))
  idx <- np$where(keep)[[1]]
  if (py_to_r(idx$size) == 0) stop("No overlap between chunk barcodes and adata.obs_names")
  ad <- ad[idx, ]
}

obs <- py_to_r(ad$obs)
if (!(opt$pert_col %in% colnames(obs))) stop(sprintf("perturbation column '%s' not found", opt$pert_col))

if (opt$max_cells > 0 && nrow(obs) > opt$max_cells) {
  keep <- sample(seq_len(nrow(obs)), opt$max_cells)
  ad <- ad[r_to_py(as.integer(keep - 1L)), ]
  obs <- py_to_r(ad$obs)
}

cell_ids <- py_to_r(ad$obs_names$to_list())
gene_ids <- py_to_r(ad$var_names$to_list())
rownames(obs) <- cell_ids

if (!(opt$control_label %in% unique(as.character(obs[[opt$pert_col]])))) {
  stop(sprintf("control label '%s' not present in '%s'", opt$control_label, opt$pert_col))
}

counts_per <- sort(table(as.character(obs[[opt$pert_col]])), decreasing = TRUE)
perts <- names(counts_per[counts_per >= opt$min_cells_per_perturbation])
perts <- perts[perts != opt$control_label]
if (opt$max_perturbations > 0 && length(perts) > opt$max_perturbations) {
  perts <- perts[seq_len(opt$max_perturbations)]
}

score_parts <- list()
summary_rows <- list()

for (pg in perts) {
  keep_pg <- as.character(obs[[opt$pert_col]]) %in% c(pg, opt$control_label)
  if (!any(keep_pg)) next

  idx_pg <- as.integer(which(keep_pg) - 1L)
  ad_pg <- ad[r_to_py(idx_pg), ]
  x_pg <- ad_pg$X
  shape_pg <- as.integer(py_to_r(x_pg$shape))

  if (isTRUE(sp$issparse(x_pg))) {
    # Convert via CSC to avoid the large extra row/col COO vectors.
    csc <- x_pg$tocsc()
    mat_pg <- sparseMatrix(
      i = as.integer(py_to_r(csc$indices)) + 1L,
      p = as.integer(py_to_r(csc$indptr)),
      x = as.numeric(py_to_r(csc$data)),
      dims = shape_pg,
      index1 = TRUE
    )
  } else {
    mat_pg <- as.matrix(py_to_r(x_pg))
  }

  obs_pg <- obs[keep_pg, , drop = FALSE]
  cell_ids_pg <- cell_ids[keep_pg]
  rownames(obs_pg) <- cell_ids_pg
  rownames(mat_pg) <- cell_ids_pg
  colnames(mat_pg) <- gene_ids

  so_pg <- CreateSeuratObject(counts = t(mat_pg), meta.data = obs_pg)
  bc_pg <- data.frame(
    cell = cell_ids_pg,
    barcode = as.character(obs_pg[[opt$pert_col]]),
    sgrna = as.character(obs_pg[[opt$pert_col]]),
    gene = as.character(obs_pg[[opt$pert_col]]),
    read_count = 1L,
    umi_count = 1L,
    stringsAsFactors = FALSE
  )

  eff_obj <- scmageck_eff_estimate(
    rds_object = so_pg,
    bc_frame = bc_pg,
    perturb_gene = pg,
    non_target_ctrl = opt$control_label
  )
  rds_sub <- eff_obj$rds
  col_name <- paste0(pg, "_eff")
  if (!(col_name %in% colnames(rds_sub@meta.data))) {
    rm(list = intersect(c("ad_pg", "x_pg", "csc", "mat_pg", "so_pg", "bc_pg", "eff_obj", "rds_sub"), ls()))
    gc()
    next
  }
  v <- rds_sub@meta.data[[col_name]]
  score_parts[[pg]] <- data.frame(
    cell_barcode = rownames(rds_sub@meta.data),
    perturbation = pg,
    ps_score = as.numeric(v),
    stringsAsFactors = FALSE
  )
  summary_rows[[pg]] <- data.frame(
    perturbation = pg,
    n_cells = as.integer(counts_per[[pg]]),
    mean_ps_score = mean(v, na.rm = TRUE),
    sd_ps_score = stats::sd(v, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  rm(list = intersect(c("ad_pg", "x_pg", "csc", "mat_pg", "so_pg", "bc_pg", "eff_obj", "rds_sub"), ls()))
  gc()
}

scores <- if (length(score_parts) > 0) do.call(rbind, score_parts) else data.frame()
summ <- if (length(summary_rows) > 0) do.call(rbind, summary_rows) else data.frame()

write.table(scores, gzfile(file.path(opt$outdir, "cell_scores.tsv.gz")), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summ, gzfile(file.path(opt$outdir, "perturbation_summary.tsv.gz")), sep = "\t", quote = FALSE, row.names = FALSE)

meta_out <- list(
  method = "PS_scMAGeCK",
  h5ad = opt$h5ad,
  chunk_cells = opt$chunk_cells,
  chunk_id = opt$chunk_id,
  pert_col = opt$pert_col,
  control_label = opt$control_label,
  n_cells = nrow(obs),
  n_perturbations_tested = length(perts),
  max_cells = opt$max_cells,
  package_versions = list(
    R = as.character(getRversion()),
    Seurat = as.character(utils::packageVersion("Seurat")),
    scMAGeCK = as.character(utils::packageVersion("scMAGeCK"))
  )
)
writeLines(toJSON(meta_out, pretty = TRUE, auto_unbox = TRUE), con = file.path(opt$outdir, "method_meta.json"))
writeLines("ok", con = file.path(opt$outdir, "done.txt"))
