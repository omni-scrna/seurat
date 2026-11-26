# Edge weights for a kNN graph, on its neighbours (knn.R --connectivities).
#
# umap_connectivities(nn, emb): UMAP's fuzzy simplicial set (uwot::similarity_graph,
# the algorithm umap-learn's fuzzy_simplicial_set implements and scanpy and
# rapids-singlecell write as their connectivities), built from the neighbours in `nn`
# (a kNN indicator graph, each cell itself included, as FindNeighbors returns it) and
# their Euclidean distances in `emb`. The search is not rerun.
#
# Depends on Matrix and uwot (>= 0.1.15, which has similarity_graph; Seurat imports uwot).

umap_connectivities <- function(nn, emb) {
  if (utils::packageVersion("uwot") < "0.1.15")
    stop("--connectivities umap needs uwot >= 0.1.15 (similarity_graph), found ", utils::packageVersion("uwot"))
  s <- summary(as(nn, "dgCMatrix"))  # i, j (1-based)
  d <- sqrt(rowSums((emb[s$i, , drop = FALSE] - emb[s$j, , drop = FALSE])^2))
  k <- tabulate(s$i, nbins = nrow(nn))
  if (any(k != k[1])) stop("umap_connectivities: every cell needs the same number of neighbours")
  ord <- order(s$i, d)  # each cell's row, nearest first (itself, at 0)
  idx <- matrix(s$j[ord], ncol = k[1], byrow = TRUE)
  dst <- matrix(d[ord], ncol = k[1], byrow = TRUE)
  g <- uwot::similarity_graph(nn_method = list(idx = idx, dist = dst))
  g <- as(g, "dgCMatrix")
  dimnames(g) <- dimnames(nn)
  g
}
