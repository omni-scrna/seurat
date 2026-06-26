#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(HDF5Array)
  library(Seurat)
  library(Matrix)
})

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
set.seed(args$seed)

# Load neighbors graph into Seurat Object
neighbors_mat <- TENxMatrix(
  args$neighbors_h5,
  group = "matrix"
)

neighbors_mat <- as(neighbors_mat, "dgCMatrix")

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
  random.seed = args$seed,
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
