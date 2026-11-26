suppressPackageStartupMessages({
  library(rhdf5)
  library(Matrix)
})

source("../../src/read_neighbors.R")
source("../../src/write_neighbors.R")

test_that("write_neighbors round-trips through read_neighbors, rows as CSR", {
  ids <- c("AAAC-1", "AAAG-1", "AACT-1", "AAGT-1")
  # an asymmetric distances graph: row i holds cell i's neighbours
  d <- sparseMatrix(i = c(1, 1, 2, 3, 4), j = c(2, 3, 1, 4, 1), x = c(0.5, 1.5, 0.5, 2, 3),
                    dims = c(4, 4), dimnames = list(ids, ids))
  snn <- sparseMatrix(i = c(1, 2, 3, 4, 1, 2), j = c(1, 2, 3, 4, 2, 1), x = c(1, 1, 1, 1, 0.3, 0.3),
                      dims = c(4, 4), dimnames = list(ids, ids))
  path <- tempfile(fileext = ".h5")
  write_neighbors(path, ids, d, snn)

  expect_identical(as.character(h5read(path, "cell_ids")), ids)
  # read_neighbors reads CSR rows: reading back gives the same matrix, not its transpose
  back_d <- read_neighbors(path, group = "")
  back_c <- read_neighbors(path, group = "connectivities")
  expect_equal(as.matrix(back_d), as.matrix(d), ignore_attr = TRUE)
  expect_equal(as.matrix(back_c), as.matrix(snn), ignore_attr = TRUE)
  # 0-based indptr, first entry 0, last nnz
  ip <- h5read(path, "indptr")
  expect_equal(ip[1], 0)
  expect_equal(ip[length(ip)], length(d@x))
})

test_that("knn_distances: Euclidean distance per kNN pair, the diagonal dropped", {
  emb <- matrix(c(0, 3, 0, 0, 4, 1), nrow = 3, dimnames = list(c("a", "b", "c"), NULL))
  nn <- sparseMatrix(i = c(1, 1, 2, 2, 3, 3), j = c(1, 2, 2, 1, 3, 1), x = 1, dims = c(3, 3),
                     dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
  d <- knn_distances(nn, emb)
  expect_equal(unname(diag(as.matrix(d))), c(0, 0, 0))
  expect_equal(d[1, 2], 5)
  expect_equal(d[3, 1], 1)
  expect_equal(length(d@x), 3)
})
