suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
})

source("../../src/connectivities.R")

.toy <- function(n = 60, d = 5) {
  set.seed(1)
  m <- rbind(matrix(rnorm(n / 3 * d, 0, 0.3), ncol = d), matrix(rnorm(n / 3 * d, 4, 0.3), ncol = d),
             matrix(rnorm(n / 3 * d, 8, 0.3), ncol = d))
  rownames(m) <- paste0("c", seq_len(nrow(m)))
  m
}

test_that("Seurat's SNN is the Jaccard index of neighbour sets, itself included, pruned at 1/15", {
  # the snn definition the scanpy, rapids-singlecell and brisc kNN modules implement
  emb <- .toy()
  g <- FindNeighbors(emb, k.param = 10, nn.method = "rann", prune.SNN = 1 / 15, verbose = FALSE)
  nn <- as(g$nn, "dgCMatrix")
  sets <- lapply(seq_len(nrow(nn)), function(i) which(nn[i, ] != 0))
  want <- outer(seq_along(sets), seq_along(sets), Vectorize(function(i, j)
    length(intersect(sets[[i]], sets[[j]])) / length(union(sets[[i]], sets[[j]]))))
  want[want < 1 / 15] <- 0
  expect_equal(unname(as.matrix(g$snn)), want, tolerance = 1e-6)
  expect_true(all(diag(as.matrix(g$snn)) == 1))
})

test_that("umap_connectivities is uwot's fuzzy simplicial set on the same neighbours", {
  emb <- .toy()
  g <- FindNeighbors(emb, k.param = 10, nn.method = "rann", verbose = FALSE)
  got <- umap_connectivities(g$nn, emb)
  # the same neighbours, ranked by brute-force distance, straight into uwot
  dm <- as.matrix(dist(emb))
  idx <- t(apply(dm, 1, function(r) order(r)[1:10]))
  dst <- t(sapply(seq_len(nrow(dm)), function(i) dm[i, idx[i, ]]))
  want <- uwot::similarity_graph(nn_method = list(idx = idx, dist = dst))
  expect_equal(unname(as.matrix(got)), unname(as.matrix(want)), tolerance = 1e-6)
  expect_equal(as.matrix(got), t(as.matrix(got)), tolerance = 1e-6, ignore_attr = TRUE)
  expect_true(all(got@x > 0 & got@x <= 1))
  expect_identical(rownames(got), rownames(emb))
})
