#!/usr/bin/env Rscript
# Argument parser for omnibenchmark seurat modules.

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
})

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
    make_option("--properties.info",   type = "character", default = NULL,
                help = "YAML file with batch_var, sample_var, labels_var fields")
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
    output_dir       = raw$output_dir,
    name             = raw$name,
    input_h5         = raw[["normalized.h5"]],
    rawdata_h5ad     = raw[["rawdata.h5ad"]],
    filtered_cellids = raw[["filtered.cellids"]],
    selection_type   = raw$selection_type,
    number_selected  = raw$number_selected,
    properties_info  = raw[["properties.info"]],
    batch_variable   = NULL
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

  if (args$selection_type == "seurat_vst_batch") {
    if (is.null(args$properties_info)) {
      stop("--properties.info is required for selection_type 'seurat_vst_batch'")
    }
    props <- yaml::read_yaml(args$properties_info)
    if (is.null(props$batch_var) || props$batch_var == "") {
      stop("batch_var is required in properties.info for selection_type 'seurat_vst_batch'")
    }
    args$batch_variable <- props$batch_var
  }

  args
}

build_integrate_parser <- function() {
  option_list <- list(
    make_option("--output_dir",             type = "character",
                help = "Output directory for results"),
    make_option("--name",                   type = "character",
                help = "Module name/identifier"),
    make_option("--pcas.tsv",               type = "character",
                help = "Global PCA embeddings (cell_id, PC1..PCn)"),
    make_option("--loadings.tsv",           type = "character",
                help = "Global PCA loadings (gene, PC1..PCn)"),
    make_option("--normalized_selected.h5", type = "character",
                help = "TENx-format HDF5 of normalized, selected expression"),
    make_option("--rawdata.h5ad",           type = "character",
                help = "AnnData h5ad (obs read for batch labels)"),
    make_option("--properties.info",        type = "character",
                help = "YAML file with batch_var field"),
    make_option("--method",                 type = "character",
                help = "Integration method (rpca, fastmnn)"),
    make_option("--k_anchor",               type = "integer", default = 5L,
                help = "Number of anchors per batch pair (RPCA only)")
  )
  OptionParser(
    option_list = option_list,
    description = "OmniBenchmark integration module (Seurat IntegrateLayers)"
  )
}

parse_integrate_args <- function() {
  parser <- build_integrate_parser()
  raw <- parse_args(parser)

  args <- list(
    output_dir      = raw$output_dir,
    name            = raw$name,
    pcas_tsv        = raw[["pcas.tsv"]],
    loadings_tsv    = raw[["loadings.tsv"]],
    input_h5        = raw[["normalized_selected.h5"]],
    rawdata_h5ad    = raw[["rawdata.h5ad"]],
    properties_info = raw[["properties.info"]],
    method          = raw$method,
    k_anchor        = raw$k_anchor
  )

  required <- c("output_dir", "name", "pcas_tsv", "loadings_tsv",
                "input_h5", "rawdata_h5ad", "properties_info", "method")
  missing <- required[vapply(args[required], function(v) is.null(v) || is.na(v),
                             logical(1))]
  if (length(missing) > 0) {
    stop("Missing required argument(s): ", paste(missing, collapse = ", "))
  }

  valid_methods <- c("rpca", "fastmnn")
  if (!(args$method %in% valid_methods)) {
    stop("Invalid --method: ", args$method,
         " (valid: ", paste(valid_methods, collapse = ", "), ")")
  }

  props <- yaml::read_yaml(args$properties_info)
  if (is.null(props$batch_var) || props$batch_var == "") {
    stop("batch_var is required in properties.info for integration")
  }
  args$batch_variable <- props$batch_var

  args
}

build_pca_parser <- function() {
  option_list <- list(
    make_option("--output_dir", type = "character",
                help = "Output directory for results"),
    make_option("--name", type = "character",
                help = "Module name/identifier"),
    make_option("--normalized_selected.h5", type = "character",
                help = "TENx-format HDF5 of normalized expression (genes x cells)"),
    make_option("--solver", type = "character",
                help = "PCA solver (seurat: exact, approximate)"),
    make_option("--n_components", type = "integer",
                help = "Number of principal components to compute"),
    make_option("--random_seed", type = "integer",
                help = "Seed for randomized solvers (and for reproducibility)")
  )

  OptionParser(
    option_list = option_list,
    description = "OmniBenchmark PCA module (Seurat)"
  )
}

parse_pca_args <- function() {
  parser <- build_pca_parser()
  raw <- parse_args(parser)

  args <- list(
    output_dir      = raw$output_dir,
    name            = raw$name,
    input_h5        = raw[["normalized_selected.h5"]],
    solver          = raw$solver,
    n_components    = raw$n_components,
    random_seed     = raw$random_seed
  )

  required <- names(args)
  missing <- required[vapply(args[required], function(v) is.null(v) || is.na(v),
                             logical(1))]
  if (length(missing) > 0) {
    stop("Missing required argument(s): ", paste(missing, collapse = ", "))
  }

  valid_solvers <- c("exact", "approximate")
  if (!(args$solver %in% valid_solvers)) {
    stop("Invalid --solver: ", args$solver,
         " (valid: ", paste(valid_solvers, collapse = ", "), ")")
  }

  args
}
