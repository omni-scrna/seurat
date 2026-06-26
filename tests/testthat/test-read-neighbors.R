suppressPackageStartupMessages({
  library(rhdf5)
  library(Matrix)
})

# read_neighbors() is defined in src/read_neighbors.R. test_dir() runs with the
# module root as the working directory.
source("../../src/read_neighbors.R")

fixture <- "fixtures/data_neighbors.h5"

test_that("fixture is present", {
  expect_true(file.exists(fixture))
})

test_that("read_neighbors loads the connectivities graph as a dgCMatrix", {
  m <- read_neighbors(fixture, group = "connectivities")

  # column-compressed sparse, which is what Seurat::as.Graph expects
  expect_s4_class(m, "dgCMatrix")

  # square, one row/col per cell barcode in /cell_ids
  ids <- as.character(h5read(fixture, "cell_ids"))
  expect_identical(dim(m), c(length(ids), length(ids)))
  expect_identical(rownames(m), ids)
  expect_identical(colnames(m), ids)

  # nnz matches the stored CSR data vector exactly (no values dropped/added)
  expect_identical(length(m@x), length(h5read(fixture, "connectivities/data")))

  # connectivities is an undirected affinity graph -> symmetric, no NA/negatives
  expect_true(isSymmetric(m))
  expect_false(anyNA(m@x))
  expect_true(all(m@x > 0))
})

test_that("read_neighbors reconstructs the CSR triples correctly", {
  m <- read_neighbors(fixture, group = "connectivities")

  # rebuild independently from the raw CSR triples (row-compressed, 0-based) and
  # compare against the reader's output value-for-value.
  data    <- as.numeric(h5read(fixture, "connectivities/data"))
  indices <- as.integer(h5read(fixture, "connectivities/indices"))
  indptr  <- as.integer(h5read(fixture, "connectivities/indptr"))
  n       <- nrow(m)

  ref <- new("dgRMatrix", x = data, j = indices, p = indptr, Dim = c(n, n))
  ref <- as(ref, "CsparseMatrix")

  # identical sparsity pattern and values (drop dimnames for the comparison)
  expect_equal(unname(m), ref)
})

test_that("the default group is connectivities", {
  expect_identical(
    read_neighbors(fixture),
    read_neighbors(fixture, group = "connectivities")
  )
})

test_that("read_neighbors can also read the flat root distances graph", {
  # the file also stores a flat distances graph at the root (group = "")
  d <- read_neighbors(fixture, group = "")
  expect_s4_class(d, "dgCMatrix")
  expect_identical(dim(d), dim(read_neighbors(fixture)))
  # distances graph is sparser than the connectivities graph in this fixture
  expect_lt(length(d@x), length(read_neighbors(fixture)@x))
})
