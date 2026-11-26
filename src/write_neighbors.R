# Write a *_neighbors.h5 file in the layout the other omni-scrna modules write and
# read_neighbors() reads: /cell_ids (fixed-length strings), the distances CSR flat at
# the file root (data, indices, indptr; 0-based, scipy's layout) and the
# connectivities CSR under /connectivities. No `shape`: both are square over cell_ids.
#
# Depends on rhdf5 and Matrix; callers attach both.

# A sparse matrix's CSR triples, rows in order and columns sorted within each row.
.csr <- function(m) {
  t_m <- as(t(m), "CsparseMatrix")  # CSC of t(m) is CSR of m
  t_m <- as(t_m, "dgCMatrix")
  list(data = as.numeric(t_m@x), indices = as.integer(t_m@i), indptr = as.integer(t_m@p))
}

.write_csr <- function(path, group, m) {
  csr <- .csr(m)
  prefix <- if (nzchar(group)) paste0(group, "/") else ""
  if (nzchar(group)) h5createGroup(path, group)
  for (k in names(csr)) h5write(csr[[k]], path, paste0(prefix, k))
}

write_neighbors <- function(path, cell_ids, distances, connectivities) {
  if (file.exists(path)) file.remove(path)
  h5createFile(path)
  # fixed-length strings: what h5py writes (dtype "S") and what Python readers expect
  ids <- as.character(cell_ids)
  h5createDataset(path, "cell_ids", dims = length(ids), storage.mode = "character",
                  size = max(nchar(ids, type = "bytes")) + 1L, chunk = NULL, level = 0)
  h5write(ids, path, "cell_ids")
  .write_csr(path, "", distances)
  .write_csr(path, "connectivities", connectivities)
  h5closeAll()
  invisible(path)
}

# Euclidean distances for the nonzero pairs of a kNN indicator graph, from the embedding
# rows the search ran on; the diagonal (each cell itself) is dropped.
knn_distances <- function(nn, emb) {
  s <- summary(as(nn, "dgCMatrix"))  # i, j, x (1-based)
  s <- s[s$i != s$j, , drop = FALSE]
  d <- sqrt(rowSums((emb[s$i, , drop = FALSE] - emb[s$j, , drop = FALSE])^2))
  d[d == 0] <- .Machine$double.xmin  # a duplicate cell stays an edge, not a structural zero
  sparseMatrix(i = s$i, j = s$j, x = d, dims = dim(nn), dimnames = dimnames(nn))
}
