# Read a CSR-stored graph from a *_neighbors.h5 file into a (symmetric) sparse
# matrix. The graph lives under `group` as scipy CSR triples (data/indices/indptr,
# 0-based) with cell barcodes at /cell_ids. Matrix::sparseMatrix reads the
# row-compressed layout directly via p + j with index1 = FALSE and returns a
# column-compressed dgCMatrix, which is what Seurat::as.Graph expects.
# (HDF5Array::H5SparseMatrix won't work here: the group carries no shape
# dataset/attr, so it can't infer the matrix dimensions.)
#
# Depends on rhdf5 (h5read) and Matrix (sparseMatrix); callers attach both.
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
