#!/usr/bin/env Rscript
# PCA module (Seurat-backed) for omnibenchmark.
#
# Solvers:
#   approximate — Seurat::RunPCA(approx = TRUE)  [IRLBA]
#   exact       — Seurat::RunPCA(approx = FALSE) [full SVD]
#
# Phases (obkit-events.jsonl):
#   load      h5 -> dgCMatrix -> Seurat object (everything that gets the input
#             into the framework's expected in-memory form).
#   densify   ScaleData(do.scale = FALSE) — allocates a dense
#             (n_genes × n_cells) double matrix in the `scale.data` layer
#             because that's what RunPCA reads from. This is the memory
#             elephant for Seurat: ~n_genes * n_cells * 8 bytes.
#   pca       RunPCA (IRLBA or full SVD).
#   write     TSV serialization.
#
# Mean-center only (do.scale = FALSE) to match scanpy/scrapper/rapids,
# which don't divide by per-gene std. Seurat's default ScaleData does both.

suppressPackageStartupMessages({
  library(Matrix)
  library(HDF5Array)
  library(Seurat)
  library(data.table)
  library(jsonlite)  # logger emits via toJSON; satisfied transitively by Seurat
})

script_dir <- (function() {
  cargs <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", cargs)
  if (length(m) > 0) dirname(sub("^--file=", "", cargs[[m]])) else getwd()
})()
# Vendored shared-CLI engine, namespaced in `cli` (see src/common/cli.R). We own
# the parser; the helpers add the benchmark's base + five-pca I/O args.
cli <- new.env()
source(file.path(script_dir, "src", "common", "cli.R"), local = cli)
source(file.path(script_dir, "src", "obkit_logger.R"))
source(file.path(script_dir, "src", "phases.R"))


parse_pca_args <- function() {
  p <- arg_parser("PCA module (Seurat)")
  p <- cli$add_base_args(p)               # --output_dir, --name
  p <- cli$add_stage_args(p, "five-pca")  # --normalized_selected.h5
  p <- add_argument(p, "--solver",       type = "character", help = "PCA solver (exact, approximate)")
  p <- add_argument(p, "--n_components", type = "integer",   help = "Number of principal components")
  p <- add_argument(p, "--random_seed",  type = "integer",   help = "Seed for reproducibility")
  raw <- parse_args(p)

  args <- list(
    output_dir = raw$output_dir, name = raw$name,
    input_h5 = raw[["normalized_selected.h5"]],
    solver = raw$solver, n_components = raw$n_components, random_seed = raw$random_seed
  )
  # The R helpers don't enforce required/choices — check by hand.
  missing <- names(args)[vapply(args, function(v) is.null(v) || is.na(v), logical(1))]
  if (length(missing) > 0) stop("Missing required argument(s): ", paste(missing, collapse = ", "))
  if (!(args$solver %in% c("exact", "approximate")))
    stop("Invalid --solver: ", args$solver, " (valid: exact, approximate)")
  args
}


run_pca <- function(so, X, args, phase) {
  # so: Seurat object with `counts` + `data` layers populated from X.
  # X:  the same gene-by-cell sparse dgCMatrix (kept for total_var).
  set.seed(args$random_seed)

  approx <- switch(args$solver,
    approximate = TRUE,
    exact       = FALSE,
    stop("unknown solver: ", args$solver)
  )

  so <- phase("densify", function(attrs) {
    attrs$do_scale        <- FALSE
    attrs$dense_bytes_est <- as.numeric(nrow(X)) * as.numeric(ncol(X)) * 8
    ScaleData(so, features = rownames(X), do.scale = FALSE, verbose = FALSE)
  })

  so <- phase("pca", function(attrs) {
    attrs$solver       <- args$solver
    attrs$approx       <- approx
    attrs$n_components <- args$n_components
    RunPCA(so,
           features = rownames(X),
           npcs = args$n_components,
           approx = approx,
           seed.use = args$random_seed,
           verbose = FALSE)
  })

  red       <- so[["pca"]]
  embedding <- Embeddings(red)              # (n_cells, n_components)
  loadings  <- Loadings(red)                # (n_genes, n_components)
  stdev     <- Stdev(red)
  variance  <- as.numeric(stdev^2)

  # total variance: sum of per-gene variances on the centered matrix
  # (matches Seurat's PCA input). Computed sparse-safe on X (genes x cells).
  n   <- ncol(X)
  rs2 <- Matrix::rowSums(X^2)
  rs1 <- Matrix::rowSums(X)
  total_var <- sum((rs2 - rs1^2 / n) / (n - 1))

  variance_ratio <- variance / total_var
  rownames(embedding) <- colnames(X)
  colnames(embedding) <- paste0("PC", seq_len(ncol(embedding)))

  list(
    embedding      = embedding,
    loadings       = matrix(as.double(loadings), nrow = nrow(loadings)),
    variance       = as.double(variance),
    variance_ratio = as.double(variance_ratio)
  )
}

main <- function() {
  args <- parse_pca_args()
  cat(sprintf("Full command: %s\n", paste(commandArgs(trailingOnly = FALSE), collapse = " ")))
  for (k in c("output_dir", "name", "input_h5",
              "solver", "n_components", "random_seed")) {
    cat(sprintf("  %s: %s\n", k, args[[k]]))
  }

  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)
  logger_init(args$output_dir)

  loaded <- phase("load", function(attrs) {
    # TENxMatrix is lazy/disk-backed; as(., dgCMatrix) materializes it
    # in memory. CreateSeuratObject + SetAssayData wrap the sparse matrix
    # in Seurat's container — sparse-pointer-cheap but allocates metadata
    # (cell IDs as factor levels, etc.). All three steps are "ingest into
    # framework-native form" so they belong in `load`.
    m_lazy <- TENxMatrix(args$input_h5, group = "matrix")
    m_mem  <- as(m_lazy, "dgCMatrix")
    so     <- CreateSeuratObject(counts = m_mem)
    so     <- SetAssayData(so, layer = "data", new.data = m_mem)
    attrs$n_genes <- nrow(m_mem)
    attrs$n_cells <- ncol(m_mem)
    attrs$nnz     <- as.integer(length(m_mem@x))
    list(X = m_mem, so = so)
  })
  cat(sprintf("  matrix (genes x cells): %d x %d\n",
              nrow(loaded$X), ncol(loaded$X)))

  res <- run_pca(loaded$so, loaded$X, args, phase)
  cat(sprintf("  embedding: %d x %d, loadings: %d x %d\n",
              nrow(res$embedding), ncol(res$embedding),
              nrow(res$loadings),  ncol(res$loadings)))

  phase("write", function(attrs) {
    out <- file.path(args$output_dir, sprintf("%s_pcas.tsv", args$name))
    cat("output_file:", out, "\n")
    fwrite(data.frame(cell_id = rownames(res$embedding), res$embedding), out,
           sep = "\t", quote = FALSE, row.names = FALSE)
    attrs$path <- out
    cat(sprintf("  wrote: %s\n", out))
  })
}

if (sys.nframe() == 0L) {
  main()
}
