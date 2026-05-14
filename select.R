#!/usr/bin/env Rscript
# Seurat-based HVG selection module for omnibenchmark.
#
# Supported selection_type values:
#   seurat_vst       — per-dataset VST (FindVariableFeatures)
#   seurat_vst_batch — per-batch VST aggregated with SelectIntegrationFeatures

suppressPackageStartupMessages({
  library(Seurat)
  library(anndata)
  library(HDF5Array)
  library(data.table)
})

script_dir <- (function() {
  cargs <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", cargs)
  if (length(m) > 0) dirname(sub("^--file=", "", cargs[[m]])) else getwd()
})()
source(file.path(script_dir, "src", "cli.R"))

run_select <- function(args) {
  so <- read_h5ad(args$rawdata_h5ad, as = "Seurat")
  cellids <- readLines(gzfile(args$filtered_cellids))
  so <- subset(so, cells = cellids)
  cat(sprintf("  dim(so) after filtering: %d x %d\n", nrow(so), ncol(so)))

  if (args$selection_type == "seurat_vst") {
    so <- FindVariableFeatures(so, selection.method = "vst",
      nfeatures = args$number_selected)
    sel_feats <- VariableFeatures(so)

  } else if (args$selection_type == "seurat_vst_batch") {
    batch_col <- args$batch_variable
    batches <- unique(so[[batch_col, drop = TRUE]])
    cat("  batches:", paste(batches, collapse = ", "), "\n")
    seurat_list <- lapply(batches, function(b) {
      cells_b <- colnames(so)[so[[batch_col, drop = TRUE]] == b]
      sub_so <- subset(so, cells = cells_b)
      FindVariableFeatures(sub_so, selection.method = "vst",
                           nfeatures = args$number_selected)
    })
    sel_feats <- SelectIntegrationFeatures(
      seurat_list,
      nfeatures = args$number_selected
    )
  } else {
    stop("Unsupported selection_type: ", args$selection_type)
  }
}

main <- function() {
  args <- parse_select_args()
  cat(sprintf("Full command: %s\n", paste(commandArgs(trailingOnly = FALSE), collapse = " ")))
  for (k in c("output_dir", "name", "input_h5", "rawdata_h5ad",
              "filtered_cellids", "selection_type", "number_selected", "batch_variable")) {
    cat(sprintf("  %s: %s\n", k, args[[k]]))
  }

  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

  sel_feats <- run_select(args)

  m <- TENxMatrix(args$input_h5, group = "matrix")
  m <- as(m, "dgCMatrix")

  out <- file.path(args$output_dir, paste0(args$name, "_normalized_selected.h5"))
  cat("output_file:", out, "\n")
  writeTENxMatrix(m[sel_feats, ], out, group = "matrix")
  cat(sprintf("  wrote: %s\n", out))
  print(file.info(out)[, c("size", "ctime")])
}

if (sys.nframe() == 0L) {
  main()
}
