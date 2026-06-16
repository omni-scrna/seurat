#!/usr/bin/env Rscript
# BPCells format conversion stage.
#
# Reads a normalized_selected.h5 (TENx-format sparse) and writes it as a
# tar.gz'd BPCells on-disk directory format. The tarball is consumed by
# the bpcells-* solver tokens in pca.R, which untar to a temp dir and
# open_matrix_dir() into a streaming on-disk Assay — sidesteps the dense
# scale.data materialization that pins ~1 GB of RAM in the standard
# Seurat path on this dataset.
#
# Phases (obkit-events.jsonl):
#   load   h5 -> dgCMatrix
#   pack   BPCells::write_matrix_dir — packed binary format
#   tar    archive the directory into the single stage output file

suppressPackageStartupMessages({
  library(Matrix)
  library(HDF5Array)
  library(BPCells)
  library(jsonlite)
})

script_dir <- (function() {
  cargs <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", cargs)
  if (length(m) > 0) dirname(sub("^--file=", "", cargs[[m]])) else getwd()
})()
source(file.path(script_dir, "src", "cli.R"))
source(file.path(script_dir, "src", "obkit_logger.R"))
source(file.path(script_dir, "src", "phases.R"))


main <- function() {
  args <- parse_bpcells_convert_args()
  cat(sprintf("Full command: %s\n", paste(commandArgs(trailingOnly = FALSE), collapse = " ")))
  for (k in c("output_dir", "name", "input_h5")) {
    cat(sprintf("  %s: %s\n", k, args[[k]]))
  }
  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)
  logger_init(args$output_dir)

  m <- phase("load", function(attrs) {
    m_lazy <- TENxMatrix(args$input_h5, group = "matrix")
    m_mem  <- as(m_lazy, "dgCMatrix")
    attrs$n_genes <- nrow(m_mem)
    attrs$n_cells <- ncol(m_mem)
    attrs$nnz     <- as.integer(length(m_mem@x))
    m_mem
  })
  cat(sprintf("  matrix (genes x cells): %d x %d\n", nrow(m), ncol(m)))

  bp_dir <- file.path(args$output_dir, sprintf("%s_bpcells", args$name))
  phase("pack", function(attrs) {
    if (dir.exists(bp_dir)) unlink(bp_dir, recursive = TRUE)
    # write_matrix_dir() chooses sensible defaults (LZ4 compression).
    write_matrix_dir(m, bp_dir)
    attrs$dir <- bp_dir
    files <- list.files(bp_dir, full.names = TRUE)
    attrs$packed_bytes <- as.integer(sum(file.info(files)$size))
  })

  out_tar <- file.path(args$output_dir, sprintf("%s_bpcells.tar.gz", args$name))
  phase("tar", function(attrs) {
    # tar from the parent dir so the archive uses a stable relative path
    cwd <- getwd()
    on.exit(setwd(cwd), add = TRUE)
    setwd(dirname(bp_dir))
    tar(out_tar, files = basename(bp_dir), compression = "gzip",
        tar = "internal")
    unlink(bp_dir, recursive = TRUE)
    attrs$path        <- out_tar
    attrs$tar_bytes   <- as.integer(file.info(out_tar)$size)
  })
  cat(sprintf("  wrote: %s\n", out_tar))
}

if (sys.nframe() == 0L) main()
