#!/usr/bin/env Rscript
# kNN module (Seurat-backed) for omnibenchmark: Seurat::FindNeighbors on an embedding.
#
# Input:  --embedding_tsv  (cell_id, PC1..PCn), as the PCA modules write it.
# Output: {name}_neighbors.h5 (src/write_neighbors.R):
#   distances       the k-nearest-neighbour graph (each cell's k - 1 neighbours, itself
#                   excluded), weighted by Euclidean distance in the embedding. Seurat
#                   returns the neighbour sets, not their distances, alongside the SNN
#                   graph; the distances are recomputed for those pairs in `write`, so
#                   the search runs once.
#   connectivities  --connectivities native or snn (default native, the same graph):
#                   Seurat's shared-nearest-neighbour graph (Jaccard of the neighbour
#                   sets, each cell itself included, pruned below --prune_snn,
#                   diagonal 1), what FindClusters clusters on in a default Seurat
#                   workflow, and the snn definition every kNN module writes.
#                   --connectivities umap: UMAP's fuzzy simplicial set on the same
#                   neighbours (src/connectivities.R, uwot), the weights scanpy and
#                   rapids-singlecell write.
#
# Defaults are Seurat's (FindNeighbors): k = 20, annoy with 50 trees, Euclidean,
# prune.SNN = 1/15. --n_dims picks the first PCs (default: all the embedding has;
# Seurat's own default is 1:10).
#
# Phases (obkit-events.jsonl): load (TSV -> matrix), knn (FindNeighbors), write.

suppressPackageStartupMessages({
  library(Matrix)
  library(rhdf5)
  library(Seurat)
  library(data.table)
  library(jsonlite)
})

source("src/common/cli.R")
source("src/obkit_logger.R")
source("src/phases.R")
source("src/write_neighbors.R")
source("src/connectivities.R")

parse_knn_args <- function(argv = commandArgs(trailingOnly = TRUE)) {
  p <- arg_parser("kNN module")
  p <- add_base_args(p)  # --output_dir, --name
  p <- add_argument(p, "--embedding_tsv", type = "character", help = "Embedding TSV (cell_id, PC1..PCn)")
  p <- add_argument(p, "--k", type = "integer", default = 20L, help = "neighbours per cell, itself included (Seurat's k.param)")
  p <- add_argument(p, "--n_dims", type = "integer", default = 0L, help = "first PCs to use (0: all)")
  p <- add_argument(p, "--nn_method", type = "character", default = "annoy", help = "annoy or rann")
  p <- add_argument(p, "--n_trees", type = "integer", default = 50L, help = "annoy trees")
  p <- add_argument(p, "--prune_snn", type = "numeric", default = 1 / 15, help = "SNN edges below this are removed")
  p <- add_argument(p, "--random_seed", type = "integer", default = 42L, help = "random seed")
  p <- add_argument(p, "--connectivities", type = "character", default = "native",
                    help = "edge weights: native or snn (Seurat's SNN), or umap (UMAP's fuzzy weights)")
  parse_args_outputs(p, argv)  # argparser, minus --output id=PATH (api 0.8)
}

main <- function() {
  args <- parse_knn_args()
  if (!args$connectivities %in% c("native", "snn", "umap"))
    stop("--connectivities: native, snn or umap, not ", args$connectivities)
  cat(sprintf("Full command: %s\n", paste(commandArgs(trailingOnly = FALSE), collapse = " ")))
  for (k in c("output_dir", "name", "embedding_tsv", "k", "n_dims", "nn_method", "n_trees", "prune_snn", "random_seed", "connectivities"))
    cat(sprintf("  %s: %s\n", k, format(args[[k]])))

  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)
  logger_init(args$output_dir)
  set.seed(args$random_seed)

  emb <- phase("load", function(attrs) {
    t <- fread(args$embedding_tsv, sep = "\t")
    m <- as.matrix(t[, -1, with = FALSE])
    storage.mode(m) <- "double"
    rownames(m) <- as.character(t[[1]])
    if (args$n_dims > 0) m <- m[, seq_len(min(args$n_dims, ncol(m))), drop = FALSE]
    attrs$n_cells <- nrow(m)
    attrs$n_dims  <- ncol(m)
    m
  })

  graphs <- phase("knn", function(attrs) {
    attrs$k <- args$k
    attrs$nn_method <- args$nn_method
    FindNeighbors(emb, k.param = args$k, nn.method = args$nn_method, n.trees = args$n_trees,
                  annoy.metric = "euclidean", prune.SNN = args$prune_snn, compute.SNN = TRUE,
                  verbose = FALSE)
  })

  conn <- if (args$connectivities == "umap") {
    phase("umap_weights", function(attrs) umap_connectivities(graphs$nn, emb))
  } else graphs$snn

  phase("write", function(attrs) {
    oid <- if (!is.null(args$outputs[["neighbors_corrected_h5"]])) "neighbors_corrected_h5" else "neighbors_h5"
    out <- output_path(args, oid, "_neighbors.h5")  # --output id=PATH (api 0.8)
    dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
    dist <- knn_distances(graphs$nn, emb)
    write_neighbors(out, rownames(emb), dist, conn)
    attrs$path <- out
    attrs$connectivities <- args$connectivities
    attrs$distances_nnz <- length(dist@x)
    attrs$connectivities_nnz <- length(conn@x)
    cat(sprintf("  wrote: %s\n", out))
  })
}

if (sys.nframe() == 0L) {
  main()
}
