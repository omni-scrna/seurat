#!/usr/bin/env Rscript
# Argument parser for omnibenchmark seurat modules.

suppressPackageStartupMessages(library(optparse))

build_select_parser <- function() {
  option_list <- list(
    make_option("--output_dir",        type = "character",
                help = "Output directory for results"),
    make_option("--name",              type = "character",
                help = "Module name/identifier"),
    make_option("--normalized.h5",     type = "character",
                help = "TENx-format HDF5 of normalized expression (genes x cells)"),
    make_option("--rawdata.h5ad",      type = "character",
                help = "AnnData h5ad with raw counts (used by seurat_vst methods)"),
    make_option("--filtered.cellids",  type = "character",
                help = "Gzipped text file of filtered cell IDs"),
    make_option("--selection_type",    type = "character",
                help = "Gene selection method (seurat_vst, seurat_vst_batch)"),
    make_option("--number_selected",   type = "integer",
                help = "Number of highly variable genes to select"),
    make_option("--batch_variable",    type = "character", default = NULL,
                help = "obs column name used as batch variable (seurat_vst_batch only)")
  )

  OptionParser(
    option_list = option_list,
    description = "OmniBenchmark gene selection module (Seurat)"
  )
}

parse_select_args <- function() {
  parser <- build_select_parser()
  raw <- parse_args(parser)

  args <- list(
    output_dir      = raw$output_dir,
    name            = raw$name,
    input_h5        = raw[["normalized.h5"]],
    rawdata_h5ad    = raw[["rawdata.h5ad"]],
    filtered_cellids = raw[["filtered.cellids"]],
    selection_type  = raw$selection_type,
    number_selected = raw$number_selected,
    batch_variable  = raw$batch_variable
  )

  required <- c("output_dir", "name", "input_h5", "rawdata_h5ad",
                "filtered_cellids", "selection_type", "number_selected")
  missing <- required[vapply(args[required], function(v) is.null(v) || is.na(v),
                             logical(1))]
  if (length(missing) > 0) {
    stop("Missing required argument(s): ", paste(missing, collapse = ", "))
  }

  valid_types <- c("seurat_vst", "seurat_vst_batch")
  if (!(args$selection_type %in% valid_types)) {
    stop("Invalid --selection_type: ", args$selection_type,
         " (valid: ", paste(valid_types, collapse = ", "), ")")
  }

  if (args$selection_type == "seurat_vst_batch" && is.null(args$batch_variable)) {
    stop("--batch_variable is required for selection_type 'seurat_vst_batch'")
  }

  args
}
