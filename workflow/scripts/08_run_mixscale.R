#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
  library(Matrix)
  library(Seurat)
  library(Mixscale)
  library(reticulate)
})

option_list <- list(
  make_option("--h5ad", type = "character"),
  make_option("--outdir", type = "character"),
  make_option("--pert-col", type = "character", default = "gene_target", dest = "pert_col"),
  make_option("--control-label", type = "character", default = "Non-Targeting", dest = "control_label"),
  make_option("--min-cells-per-perturbation", type = "integer", default = 30, dest = "min_cells_per_perturbation"),
  make_option("--max-perturbations", type = "integer", default = 100, dest = "max_perturbations"),
  make_option("--max-cells", type = "integer", default = 0, dest = "max_cells"),
  make_option("--random-seed", type = "integer", default = 0, dest = "random_seed")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$h5ad) || is.null(opt$outdir)) stop("--h5ad and --outdir are required")
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
set.seed(opt$random_seed)

anndata <- import("anndata")
sp <- import("scipy.sparse")
ad <- anndata$read_h5ad(opt$h5ad)

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
so <- NormalizeData(so)
so <- FindVariableFeatures(so)
so <- ScaleData(so, verbose = FALSE)
so <- RunPCA(so, npcs = 40, verbose = FALSE)

so <- CalcPerturbSig(
  object = so,
  assay = "RNA",
  slot = "data",
  gd.class = opt$pert_col,
  nt.cell.class = opt$control_label,
  reduction = "pca",
  ndims = 40,
  num.neighbors = 20,
  new.assay.name = "PRTB",
  split.by = NULL
)

so <- RunMixscale(
  object = so,
  assay = "PRTB",
  slot = "scale.data",
  labels = opt$pert_col,
  nt.class.name = opt$control_label,
  min.de.genes = 5,
  logfc.threshold = 0.2,
  de.assay = "RNA",
  max.de.genes = 100,
  new.class.name = "mixscale_score",
  fine.mode = FALSE,
  verbose = FALSE,
  split.by = NULL
)

meta <- so@meta.data
cell_scores <- data.frame(
  cell_barcode = rownames(meta),
  perturbation = as.character(meta[[opt$pert_col]]),
  mixscale_score = as.numeric(meta$mixscale_score),
  stringsAsFactors = FALSE
)

counts_per <- sort(table(cell_scores$perturbation), decreasing = TRUE)
perts <- names(counts_per[counts_per >= opt$min_cells_per_perturbation])
perts <- perts[perts != opt$control_label]
if (length(perts) > opt$max_perturbations) perts <- perts[seq_len(opt$max_perturbations)]

de_tbl <- data.frame()
if (length(perts) > 0) {
  de_res <- Run_wmvRegDE(
    object = so,
    assay = "RNA",
    slot = "counts",
    labels = opt$pert_col,
    nt.class.name = opt$control_label,
    PRTB_list = perts,
    logfc.threshold = 0.2,
    split.by = NULL
  )
  parts <- lapply(names(de_res), function(nm) {
    df <- de_res[[nm]]
    if (!is.data.frame(df)) return(NULL)
    df$perturbation <- nm
    df
  })
  parts <- Filter(Negate(is.null), parts)
  if (length(parts) > 0) {
    de_tbl <- do.call(rbind, parts)
  }
}

write.table(cell_scores, gzfile(file.path(opt$outdir, "cell_scores.tsv.gz")), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(de_tbl, gzfile(file.path(opt$outdir, "perturbation_de.tsv.gz")), sep = "\t", quote = FALSE, row.names = FALSE)

meta_out <- list(
  method = "Mixscale",
  h5ad = opt$h5ad,
  pert_col = opt$pert_col,
  control_label = opt$control_label,
  n_cells = nrow(cell_scores),
  n_perturbations_tested = length(perts),
  max_cells = opt$max_cells,
  package_versions = list(
    R = as.character(getRversion()),
    Seurat = as.character(utils::packageVersion("Seurat")),
    Mixscale = as.character(utils::packageVersion("Mixscale"))
  )
)
writeLines(toJSON(meta_out, pretty = TRUE, auto_unbox = TRUE), con = file.path(opt$outdir, "method_meta.json"))
writeLines("ok", con = file.path(opt$outdir, "done.txt"))
