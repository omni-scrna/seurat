#!/usr/bin/env Rscript
# Seurat-based HVG selection module for omnibenchmark.
#
# Supported selection_type values:
#   seurat_vst       — per-dataset VST (FindVariableFeatures)
#   seurat_vst_batch — per-batch VST aggregated with SelectIntegrationFeatures

suppressPackageStartupMessages({
  library(Seurat)
  library(anndataR)
  library(data.table)
  library(HDF5Array)
  library(yaml)
})

# arg parsing
source("src/common/cli.R")
p <- arg_parser("FEAT module")
p <- add_base_args(p)                    # --output_dir, --name
p <- add_stage_args(p, "FEAT")     # the stage I/O contract
# your own method params — argparser directly (its add_argument requires `help`):
p <- add_argument(p, "--selection_type", type = "character", help = "type of feature selection")
p <- add_argument(p, "--number_selected", type = "integer", help = "number of PCs")
args <- parse_args(p)                    # argparser's own parser



# from properties input, get batch variable
props <- yaml::read_yaml(args$properties.info)
if (is.null(props$batch_var) || props$batch_var == "") {
  stop("batch_var is required in properties.info for selection_type 'seurat_vst_batch'")
}
args$batch_variable <- props$batch_var

# logging
cat(sprintf("Full command: %s\n", paste(commandArgs(trailingOnly = FALSE), collapse = " ")))
cat(sprintf("LOG: command line args\n----------------------------------\n"))
for (i in 1:length(args)) {
  cat(sprintf("  %s: %s\n", names(args)[i], args[[i]]))
}
cat(sprintf("----------------------------------\n"))


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
