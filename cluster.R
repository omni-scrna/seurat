#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(HDF5Array)
  library(rhdf5)
  library(Seurat)
  library(Matrix)
})

# read_neighbors(): load a CSR neighbors graph from a *_neighbors.h5 file.
source("src/read_neighbors.R")
# flavor_algorithm_id(): map a --flavor name to a Seurat FindClusters algorithm id.
source("src/flavor.R")

# arg parsing
source("src/common/cli.R")
p <- arg_parser("CLUST module")
p <- add_base_args(p)                      # --output_dir, --name
p <- add_stage_args(p, "CLUST")  # the stage I/O contract
# your own method params — argparser directly (its add_argument requires `help`):
p <- add_argument(p, "--flavor", type = "character", help = "Clustering algorithm (louvain_original|louvain_multilevel_refinement|slm|leiden)")
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


# Map the requested --flavor to a Seurat FindClusters algorithm id (errors if
# --flavor is missing or unknown).
algorithm_seurat_id <- flavor_algorithm_id(args$flavor)


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
cat("Selected algorithm:", args$flavor, "\n")
cat("Algorithm ID:", algorithm_seurat_id, "\n\n")


# Extract clusters matrix
m_clusters <- data.frame(cell_id = colnames(so),
                         cluster = as.character(so$seurat_clusters))

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
  row.names = FALSE
)

print(file.info(output_file)[, c("size", "ctime")])
