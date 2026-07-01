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
ad <- anndata$read_h5ad(opt$h5ad)

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
  ad <- ad[keep, ]
  obs <- py_to_r(ad$obs)
}

x <- ad$X
shape <- py_to_r(x$shape)
if (isTRUE(sp$issparse(x))) {
  coo <- x$tocoo()
  mat <- sparseMatrix(
    i = as.integer(py_to_r(coo$row)) + 1L,
    j = as.integer(py_to_r(coo$col)) + 1L,
    x = as.numeric(py_to_r(coo$data)),
    dims = as.integer(shape)
  )
} else {
  mat <- as.matrix(py_to_r(x))
}

cell_ids <- py_to_r(ad$obs_names$to_list())
gene_ids <- py_to_r(ad$var_names$to_list())
rownames(obs) <- cell_ids
rownames(mat) <- cell_ids
colnames(mat) <- gene_ids

if (!(opt$control_label %in% unique(as.character(obs[[opt$pert_col]])))) {
  stop(sprintf("control label '%s' not present in '%s'", opt$control_label, opt$pert_col))
}

so <- CreateSeuratObject(counts = t(mat), meta.data = obs)

bc_frame <- data.frame(
  cell = rownames(obs),
  barcode = as.character(obs[[opt$pert_col]]),
  sgrna = as.character(obs[[opt$pert_col]]),
  gene = as.character(obs[[opt$pert_col]]),
  read_count = 1L,
  umi_count = 1L,
  stringsAsFactors = FALSE
)

counts_per <- sort(table(bc_frame$gene), decreasing = TRUE)
perts <- names(counts_per[counts_per >= opt$min_cells_per_perturbation])
perts <- perts[perts != opt$control_label]
if (opt$max_perturbations > 0 && length(perts) > opt$max_perturbations) {
  perts <- perts[seq_len(opt$max_perturbations)]
}

score_parts <- list()
summary_rows <- list()

for (pg in perts) {
  eff_obj <- scmageck_eff_estimate(
    rds_object = so,
    bc_frame = bc_frame,
    perturb_gene = pg,
    non_target_ctrl = opt$control_label
  )
  rds_sub <- eff_obj$rds
  col_name <- paste0(pg, "_eff")
  if (!(col_name %in% colnames(rds_sub@meta.data))) next
  v <- rds_sub@meta.data[[col_name]]
  score_parts[[pg]] <- data.frame(
    cell_barcode = rownames(rds_sub@meta.data),
    perturbation = pg,
    ps_score = as.numeric(v),
    stringsAsFactors = FALSE
  )
  summary_rows[[pg]] <- data.frame(
    perturbation = pg,
    n_cells = sum(as.character(obs[[opt$pert_col]]) == pg),
    mean_ps_score = mean(v, na.rm = TRUE),
    sd_ps_score = stats::sd(v, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
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
