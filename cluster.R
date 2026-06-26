#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(HDF5Array)
  library(rhdf5)
  library(Seurat)
  library(Matrix)
})

# Read a CSR-stored graph from a *_neighbors.h5 file into a (symmetric) sparse
# matrix. The graph lives under `group` as scipy CSR triples (data/indices/indptr,
# 0-based) with cell barcodes at /cell_ids. Matrix::sparseMatrix reads the
# row-compressed layout directly via p + j with index1 = FALSE and returns a
# column-compressed dgCMatrix, which is what Seurat::as.Graph expects.
# (HDF5Array::H5SparseMatrix won't work here: the group carries no shape
# dataset/attr, so it can't infer the matrix dimensions.)
read_neighbors <- function(path, group = "connectivities") {
  ids <- as.character(h5read(path, "cell_ids"))
  n <- length(ids)

  sparseMatrix(
    p    = h5read(path, paste0(group, "/indptr")),
    j    = h5read(path, paste0(group, "/indices")),  # 0-based column indices
    x    = as.numeric(h5read(path, paste0(group, "/data"))),
    dims = c(n, n),
    index1   = FALSE,
    dimnames = list(ids, ids)
  )
}

# arg parsing
source("src/common/cli.R")
p <- arg_parser("CLUST module")
p <- add_base_args(p)                      # --output_dir, --name
p <- add_stage_args(p, "CLUST")  # the stage I/O contract
# your own method params — argparser directly (its add_argument requires `help`):
p <- add_argument(p, "--modularity_algorithm", type = "character", help = "Clustering algorithm")
p <- add_argument(p, "--resolution", type = "numeric", help = "Clustering resolution")
p <- add_argument(p, "--random_seed", type = "integer", help = "Random seed")
args <- parse_args(p)                      # argparser's own parser

# logging
cat(sprintf("Full command: %s\n", paste(commandArgs(trailingOnly = FALSE), collapse = " ")))
cat(sprintf("LOG: command line args\n----------------------------------\n"))
for (i in 1:length(args)) {
  cat(sprintf("  %s: %s\n", names(args)[i], args[[i]]))
}
cat(sprintf("----------------------------------\n"))


# Reproducibility
set.seed(args$random_seed)

# Load neighbors graph into Seurat Object. Cluster on the connectivities graph
# (UMAP-style affinities); the distances graph is the flat root layout.
neighbors_mat <- read_neighbors(args$neighbors_h5, group = "connectivities")

neighbors_graph <- as.Graph(neighbors_mat)

cat("Neighbors graph dimensions:\n")
print(dim(neighbors_graph))
cat("\n")

so <- CreateSeuratObject(
  counts = neighbors_graph,
  assay = "RNA"
)

so@graphs$neighbors <- neighbors_graph


# Algorithm name mapping to Seurat ID (1,2,3, or 4)
alg_map <- list(
  louvain_original = 1,
  louvain_multilevel_refinement = 2,
  slm = 3,
  leiden = 4
)

if (!(args$modularity_algorithm %in% names(alg_map))) {
  stop("Invalid modularity_algorithm specified")
}

algorithm_seurat_id <- alg_map[[args$modularity_algorithm]]


# Run clustering
so <- FindClusters(
  so,
  algorithm = algorithm_seurat_id,
  resolution = args$resolution,
  graph.name = "neighbors",
  random.seed = args$random_seed,
  verbose = TRUE
)

cat("Running clustering...\n")
cat("Selected algorithm:", args$modularity_algorithm, "\n")
cat("Algorithm ID:", algorithm_seurat_id, "\n\n")


# Extract clusters matrix
m_clusters <- matrix(
  as.character(so$seurat_clusters),
  nrow = 1,
  dimnames = list(
    "clusters",
    colnames(so)
  )
)

cat("Cluster matrix dimensions:\n")
print(dim(m_clusters))


# Save cluster matrix as .tsv
output_file <- file.path(
  args$output_dir,
  paste0(args$name, "_clusters.tsv")
)

cat("Writing output to:\n")
cat(output_file, "\n\n")

write.table(
  m_clusters,
  file = output_file,
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)

print(file.info(output_file)[, c("size", "ctime")])
