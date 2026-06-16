#!/usr/bin/env Rscript
# PCA module (Seurat-backed) for omnibenchmark.
#
# Solvers:
#   approximate          Seurat::RunPCA(approx = TRUE)  [IRLBA]
#   exact                Seurat::RunPCA(approx = FALSE) [full SVD]
#   bpcells-approximate  same as approximate, but input is a tarball'd
#                        BPCells dir; the Assay layers stay BPCells-backed
#                        so ScaleData skips the dense materialization.
#   bpcells-exact        same as exact, but BPCells-backed input.
#
# Phases (obkit-events.jsonl):
#   load      ingest input into framework form.
#             - in-memory path: h5 -> dgCMatrix -> Seurat object
#             - bpcells path:   untar -> open_matrix_dir() -> Seurat v5 Assay
#   densify   ScaleData(do.scale = FALSE).
#             - in-memory path: allocates a dense (n_genes × n_cells) double
#               matrix in scale.data — the Seurat memory elephant
#             - bpcells path:   keeps scale.data as a BPCells-wrapped transform
#               (no dense materialization). Phase still recorded; expected
#               to be much shorter / lower-RSS.
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

# BPCells is optional — only required for bpcells-* solver tokens. Load lazily
# so the in-memory solvers don't pay the import cost (and so the module still
# runs in environments without BPCells installed).
.maybe_load_bpcells <- function() {
  if (!requireNamespace("BPCells", quietly = TRUE)) {
    stop("BPCells is not installed but a bpcells-* solver was requested. ",
         "Add r-bpcells to the seurat conda env.")
  }
}

script_dir <- (function() {
  cargs <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", cargs)
  if (length(m) > 0) dirname(sub("^--file=", "", cargs[[m]])) else getwd()
})()
source(file.path(script_dir, "src", "cli.R"))
source(file.path(script_dir, "src", "obkit_logger.R"))
source(file.path(script_dir, "src", "phases.R"))


run_pca <- function(so, X, args, phase) {
  # so: Seurat object with `counts` + `data` layers populated from X.
  # X:  the gene-by-cell input (dgCMatrix for in-memory solvers, IterableMatrix
  #     for bpcells-* solvers). Used only for total_var; rowSums/rowSums(.^2)
  #     work on both representations.
  set.seed(args$random_seed)

  approx <- switch(args$solver,
    approximate         = TRUE,
    exact               = FALSE,
    `bpcells-approximate` = TRUE,
    `bpcells-exact`       = FALSE,
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

  uses_bpcells <- startsWith(args$solver, "bpcells-")
  if (uses_bpcells) .maybe_load_bpcells()

  loaded <- phase("load", function(attrs) {
    attrs$backend <- if (uses_bpcells) "bpcells" else "in-memory"
    if (uses_bpcells) {
      # Untar the BPCells dir under a job-scoped tempdir, open_matrix_dir()
      # gives a streaming IterableMatrix. Wrap it in a v5 Assay; Seurat keeps
      # the counts/data layers BPCells-backed, so ScaleData later wraps a
      # transform rather than materializing a dense matrix.
      bp_root <- file.path(tempdir(), sprintf("%s_bpcells_in", args$name))
      if (dir.exists(bp_root)) unlink(bp_root, recursive = TRUE)
      dir.create(bp_root, recursive = TRUE)
      untar(args$bpcells_tar, exdir = bp_root)
      sub <- list.dirs(bp_root, recursive = FALSE)
      if (length(sub) != 1L) {
        stop("expected exactly one BPCells dir in tarball, found: ", length(sub))
      }
      m <- BPCells::open_matrix_dir(sub)
      assay <- CreateAssay5Object(counts = m)
      so    <- CreateSeuratObject(counts = assay)
      so    <- SetAssayData(so, layer = "data", new.data = m)
      attrs$n_genes <- nrow(m)
      attrs$n_cells <- ncol(m)
      list(X = m, so = so)
    } else {
      m_lazy <- TENxMatrix(args$input_h5, group = "matrix")
      m_mem  <- as(m_lazy, "dgCMatrix")
      so     <- CreateSeuratObject(counts = m_mem)
      so     <- SetAssayData(so, layer = "data", new.data = m_mem)
      attrs$n_genes <- nrow(m_mem)
      attrs$n_cells <- ncol(m_mem)
      attrs$nnz     <- as.integer(length(m_mem@x))
      list(X = m_mem, so = so)
    }
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
