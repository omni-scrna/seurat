#!/usr/bin/env Rscript
# Seurat RPCA batch correction module for omnibenchmark.
#
# Reads per-batch PCA embeddings + loadings, reconstructs per-batch Seurat objects,
# finds integration anchors vith RPCA, and writes corrected embeddings.

suppressPackageStartupMessages({
  library(Seurat)
  library(rhdf5)
  library(HDF5Array)
  library(data.table)
  library(yaml)
  library(assertthat)
})

# files for testing
args <- list()
args$loadings_per_batch_tsv <- "/projects/site/pred/neurogenomics/users/kodermam/omni-scrna/test_files/dataset_name-sc-mix_loadings_per_batch.tsv"
args$pcas_per_batch_tsv <- "/projects/site/pred/neurogenomics/users/kodermam/omni-scrna/test_files/dataset_name-sc-mix_pcas_per_batch.tsv"
args$output_dir <- "/projects/site/pred/neurogenomics/users/kodermam/omni-scrna/test_files"
args$input_h5 <- "/projects/site/pred/neurogenomics/users/kodermam/omni-scrna/split-stages-plan/out/one-data/datasets/dataset_name-sc-mix/two-filter/filtering-r/filter_type-scrapper-auto/three-normalize/normalization-r/normalization_type-seurat_log1pCP10k/datasets_normalized.h5"

script_dir <- (function() {
  cargs <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", cargs)
  if (length(m) > 0) dirname(sub("^--file=", "", cargs[[m]])) else getwd()
})()
source(file.path(script_dir, "src", "cli.R"))


main <- function() {
  args <- parse_rpca_args()
  cat(sprintf("Full command: %s\n", paste(commandArgs(trailingOnly = FALSE), collapse = " ")))
  for (k in names(args)) cat(sprintf("  %s: %s\n", k, args[[k]]))

  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # read matrix with normalized data
  m <- TENxMatrix(args$input_h5, group = "matrix")
  m <- as(m, "dgCMatrix")

  # read per-batch PCA embeddings and loadings
  pca_df  <- fread(args$pcas_per_batch_tsv)
  loadings_df <- fread(args$loadings_per_batch_tsv)

  # rename PC columns to match Seurat's convention  
  pc_cols <- colnames(pca_df)[grep("^PC", colnames(pca_df))]
  new_pc_cols <- gsub("^PC", "PC_", pc_cols)
  setnames(pca_df, pc_cols, new_pc_cols)
  setnames(loadings_df, pc_cols, new_pc_cols)

  batches <- unique(pca_df$batch_id)
  cat(sprintf("  batches: %s\n", paste(batches, collapse = ", ")))

  #  build per-batch Seurat objects
  object_list <- lapply(batches, function(b) {
    batch_cells <- pca_df[batch_id == b, cell_id]
    batch_embedding <- as.matrix(pca_df[batch_id == b, ..new_pc_cols])
    rownames(batch_embedding) <- batch_cells
    
    batch_genes <- loadings_df[batch_id == b, gene]
    batch_loadings <- as.matrix(loadings_df[batch_id == b, ..new_pc_cols])
    rownames(batch_loadings) <- batch_genes

    batch_norm_mat <- m[, batch_cells]
    assert_that(setequal(batch_genes, rownames(batch_norm_mat)))

    # initialize Seurat object
    rna_assay <- CreateAssayObject(data = batch_norm_mat)
    so <- CreateSeuratObject(counts = rna_assay)

    # insert custom PCA embeddings to the dummny Seurat object
    so[["pca"]] <- CreateDimReducObject(
      embeddings = batch_embedding,
      loadings   = batch_loadings,
      key        = "PC_",
      assay      = "RNA"
    )
    so
  })
  
  names(object_list) <- batches
  cat(sprintf("  built %d Seurat objects\n", length(object_list)))

  # find integration anchors with rpca
  anchors <- FindIntegrationAnchors(
    object.list = object_list,
    reduction   = "rpca",
    dims        = seq_len(args$dims),
    k.anchor    = args$k_anchor,
    anchor.features = rownames(m)
  )
  cat(sprintf("  found %d anchors\n", nrow(anchors@anchors)))
  
  # global embeddings and global loadings are wrong
  global_embeddings <- do.call(rbind, lapply(object_list, function(so) {
    Embeddings(so, "pca")
  }))
  global_loadings_df <- loadings_df[, lapply(.SD, mean), by = gene, .SDcols = new_pc_cols]
  global_loadings <- as.matrix(global_loadings_df[match(rownames(so), gene), ..new_pc_cols])
  rownames(global_loadings) <- rownames(so)
  global_loadings[is.na(global_loadings)] <- 0

  # 2. wrap the global matrix into a single unified DimReduc object
  global_pca <- CreateDimReducObject(
    embeddings = global_embeddings,
    loadings   = global_loadings,
    key        = "PC_",
    assay      = "RNA"
  )

  # 3. pass this global reduction explicitly to the function
  integrated_obj <- IntegrateEmbeddings(
    anchorset          = anchors,
    reductions         = global_pca,        # <-- Crucial fix
    new.reduction.name = "integrated.rpca",
    dims.to.integrate  = seq_len(args$dims)
  )

  corrected <- Embeddings(integrated, "rpca")
  colnames(corrected) <- pc_cols

  # write output
  out <- file.path(args$output_dir, paste0(args$name, "_corrected.tsv"))
  fwrite(data.table(cell_id = rownames(corrected), corrected), out,
         sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  wrote: %s\n", out))
}

if (sys.nframe() == 0L) {
  main()
}
