#!/usr/bin/env Rscript
# Seurat-based normalization module for omnibenchmark.
#
# Supported flavor values:
#   sctransformv2      

suppressPackageStartupMessages({
    library(Seurat)
    library(anndataR)
    library(HDF5Array)
    library(future)
})

# arg parsing
source("src/common/cli.R")
p <- arg_parser("NORM module")
p <- add_base_args(p)                      # --output_dir, --name
p <- add_stage_args(p, "NORM")  # the stage I/O contract
# your own method params — argparser directly (its add_argument requires `help`):
p <- add_argument(p, "--flavor", type = "character", help = "Normalization type")
p <- add_argument(p, "--random_seed", type = "integer", help = "random seed")
args <- parse_args(p)                      # argparser's own parser

# tracking stuff
source("src/obkit_logger.R")
source("src/phases.R")

if (!requireNamespace("glmGamPoi", quietly = TRUE)) {
  stop("glmGamPoi is required for SCTransform(method = 'glmGamPoi')")
}

run_normalize <- function(args) {
  options(future.globals.maxSize = 4 * 1024^3)
  plan(sequential)

  so <- read_h5ad(args$rawdata_h5ad, as = "Seurat")
  # potential: Feature names were replaced with dashes ('-')

  cellids <- readLines(gzfile(args$filtered_cellids))
  so <- subset(so, cells = cellids)
  cat(sprintf("  dim(so) after filtering: %d x %d\n", nrow(so), ncol(so)))

  if (args$flavor == "sctransformv2") {
    set.seed(args$random_seed)
    so <- SCTransform(so, vst.flavor = "v2", method = "glmGamPoi", 
                      verbose = FALSE, return.only.var.genes = FALSE)
    m <- GetAssayData(so, assay = "SCT", layer = "data")
    # layer = "data" for log1p(corrected UMI)
    # layer = "scale.data" for Pearson residuals
  } else {
    stop("Unsupported flavor: ", args$flavor)
  }

  return(m)
}

main <- function() {
  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

  m <- run_normalize(args)
  out <- file.path(args$output_dir, paste0(args$name, "_normalized.h5"))
  cat("output_file:", out, "\n")
  writeTENxMatrix(m, out, group = "matrix")
  cat(sprintf("  wrote: %s\n", out))
  print(file.info(out)[, c("size", "ctime")])
}

if (sys.nframe() == 0L) {
  main()
}
