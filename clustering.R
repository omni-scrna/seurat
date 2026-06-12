#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(HDF5Array)
  library(Seurat)
  library(Matrix)
})

# Parse command line arguments
parser <- ArgumentParser(description = "Seurat Clustering Module")

# Required by OmniBenchmark
parser$add_argument(
  "--output_dir",
  dest = "output_dir",
  type = "character",
  required = TRUE,
  help = "Output directory for results"
)

parser$add_argument(
  "--name",
  dest = "name",
  type = "character",
  required = TRUE,
  help = "Module name/identifier"
)

# Stage-specific inputs
parser$add_argument(
  "--neighbors.h5",
  dest = "neighbors_h5",
  type = "character",
  nargs = "+",
  required = TRUE,
  help = "Input neighbors graph"
)

parser$add_argument(
  "--modularity_algorithm",
  dest = "modularity_algorithm",
  type = "character",
  required = TRUE,
  help = "Clustering algorithm"
)

parser$add_argument(
  "--resolution",
  dest = "resolution",
  type = "double",
  help = "Clustering resolution"
)

parser$add_argument(
  "--seed",
  dest = "seed",
  type = "integer",
  default = 42,
  help = "Random seed"
)

args <- parser$parse_args()


# Reproducibility
set.seed(args$seed)

# Parameters
cat("Full command:\n")
cat(paste0(commandArgs(), collapse = " "), "\n\n")

cat("neighbors_h5:", args$neighbors_h5, "\n")
cat("output_dir:", args$output_dir, "\n")
cat("name:", args$name, "\n")
cat("modularity_algorithm:", args$modularity_algorithm, "\n")
cat("resolution:", args$resolution, "\n")
cat("seed:", args$seed, "\n\n")


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
