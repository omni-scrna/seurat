#!/usr/bin/env Rscript
# Seurat v5 batch integration module for omnibenchmark.
#
# Supports RPCA and FastMNN via IntegrateLayers.
# Per-batch PCA is computed internally by Seurat after layer splitting.

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  #library(SeuratWrappers) # skipping because object 'FastMNNIntegration' not found
  library(rhdf5)
  library(HDF5Array)
  library(data.table)
  library(yaml)
})


# arg parsing
source("src/common/cli.R")
p <- arg_parser("INTG8 module")
p <- add_base_args(p)                      # --output_dir, --name
p <- add_stage_args(p, "INTG8")  # the stage I/O contract
# your own method params — argparser directly (its add_argument requires `help`):
p <- add_argument(p, "--method", type = "integer", help = "number of PCs")
p <- add_argument(p, "--k_anchor", type = "integer", help = "number of PCs")
args <- parse_args(p)                      # argparser's own parser

# logging
cat(sprintf("Full command: %s\n", paste(commandArgs(trailingOnly = FALSE), collapse = " ")))
cat(sprintf("LOG: command line args\n----------------------------------\n"))
for (i in 1:length(args))
  cat(sprintf("  %s: %s\n", names(args)[i], args[[i]]))
cat(sprintf("----------------------------------\n"))


main <- function() {
  args <- parse_integrate_args()
  cat(sprintf("Full command: %s\n", paste(commandArgs(trailingOnly = FALSE), collapse = " ")))
  for (k in names(args)) cat(sprintf("  %s: %s\n", k, args[[k]]))

  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

  batch_var <- args$batch_variable
  cat(sprintf("  batch_variable: %s\n", batch_var))

  # read normalized selected expression matrix
  message("  Reading expression matrix")
  m <- TENxMatrix(args$input_h5, group = "matrix")
  m <- as(m, "dgCMatrix")
 
  # read batch labels from h5ad obs
  cell_ids_h5ad <- as.character(h5read(args$rawdata_h5ad, "obs/_index"))
  batch_raw <- h5read(args$rawdata_h5ad, paste0("obs/", batch_var)) |>
    as.character()
  
  # align batch labels to cells in the expression matrix
  cell_ids_mat <- colnames(m)
  idx <- match(cell_ids_mat, cell_ids_h5ad)
  batch_aligned <- batch_raw[idx]

  # create Seurat object with normalized data
  message("  Creating Seurat object")
  rna_assay <- CreateAssayObject(data = m)
  obj <- CreateSeuratObject(counts = rna_assay)
  obj[[batch_var]] <- batch_aligned

  # split into layers by batch
  message("  Splitting Seurat object into layers by batch")
  obj[["RNA"]] <- split(obj[["RNA"]], f = obj[[batch_var, drop = TRUE]])

  # read global PCA embeddings
  pca_df <- fread(args$pcas_tsv)
  pc_cols <- colnames(pca_df)[grep("^PC", colnames(pca_df))]
  embedding <- as.matrix(pca_df[, ..pc_cols])
  rownames(embedding) <- pca_df$cell_id

  # read global PCA loadings
  loadings_df <- fread(args$loadings_tsv)
  loadings_mat <- as.matrix(loadings_df[, ..pc_cols])
  rownames(loadings_mat) <- loadings_df$gene

  # rename PC columns to match Seurat convention
  pc_cols_seurat <- gsub("^PC", "PC_", pc_cols)
  colnames(embedding) <- pc_cols_seurat
  colnames(loadings_mat) <- pc_cols_seurat

  # add precomputed PCA to the object
  message("  Adding precomputed PCs to Seurat object")
  obj[["pca"]] <- CreateDimReducObject(
    embeddings = embedding,
    loadings   = loadings_mat,
    key        = "PC_",
    assay      = "RNA"
  )

  # integrate
    message("  Integrating layers using RPCA")
    batch_sizes <- table(batch_aligned)
    # set inital k_weight to prevent Error: k.weight (100) is set larger than the number of cells in the smallest object
    k_weight <- min(100, (min(batch_sizes) - 1L))

    run_integrate <- function(obj, k_weight, k_anchor) {
      IntegrateLayers(
        object         = obj,
        method         = RPCAIntegration,
        orig.reduction = "pca",
        new.reduction  = "integrated.rpca",
        k.anchor       = k_anchor,
        verbose        = TRUE,
        features       = rownames(obj),
        k.weight       = k_weight
      )
    }

    message(sprintf("  Trying k.weight = %d", k_weight))
    obj <- tryCatch(
      run_integrate(obj, k_weight, args$k_anchor),
      error = function(e) {
        msg <- conditionMessage(e)
        m <- regmatches(msg, regexpr("less than \\d+", msg))
        if (length(m) > 0) {
          new_k.weight <- as.integer(sub("less than ", "", m)) - 1L
          message(sprintf("  FindWeights failed; retrying with k.weight = %d", new_k.weight))
          run_integrate(obj, new_k.weight, args$k_anchor)
        } else {
          stop(e)
        }
      }
    )
    reduction_name <- "integrated.rpca"
  #  skipping because object 'FastMNNIntegration' not found
  #} else if (args$method == "fastmnn") { 
  #  message("  Integrating layers using FastMNN")
  #  obj <- IntegrateLayers(
  #    object        = obj,
  #    method        = FastMNNIntegration,
  #    new.reduction = "integrated.mnn",
  #    verbose       = TRUE, 
  #    features      = rownames(obj)
  #  )
  #  reduction_name <- "integrated.mnn"
  #}


  # extract corrected embeddings
  message("  Extracting corrected embeddings")
  corrected <- Embeddings(obj, reduction_name)
  colnames(corrected) <- paste0("corrected_dim", seq_len(ncol(corrected)))

  # write output
  out <- file.path(args$output_dir, paste0(args$name, "_corrected.tsv"))
  fwrite(data.table(cell_id = rownames(corrected), corrected), out,
         sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  wrote: %s\n", out))
}

if (sys.nframe() == 0L) {
  main()
}
