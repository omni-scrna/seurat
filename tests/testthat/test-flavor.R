# flavor_algorithm_id() is defined in src/flavor.R. test_dir() runs with the
# module root as the working directory.
source("../../src/flavor.R")

test_that("each flavor maps to its Seurat algorithm id", {
  expect_identical(flavor_algorithm_id("louvain_original"), 1L)
  expect_identical(flavor_algorithm_id("louvain_multilevel_refinement"), 2L)
  expect_identical(flavor_algorithm_id("slm"), 3L)
  expect_identical(flavor_algorithm_id("leiden"), 4L)
})

test_that("the returned id is an unnamed integer scalar", {
  id <- flavor_algorithm_id("slm")
  expect_type(id, "integer")
  expect_length(id, 1L)
  expect_null(names(id))
})

test_that("every mapped id is one Seurat::FindClusters accepts (1-4)", {
  expect_setequal(unname(FLAVOR_ALGORITHMS), 1:4)
})

test_that("--flavor is required: missing/NA value errors", {
  expect_error(flavor_algorithm_id(NA_character_), "--flavor must be one of")
  expect_error(flavor_algorithm_id(character(0)), "--flavor must be one of")
})

test_that("unknown flavor errors and lists the valid choices", {
  expect_error(flavor_algorithm_id("louvain"), "louvain_original")   # native alias is not a flavor
  expect_error(flavor_algorithm_id("bogus"), "--flavor must be one of")
})
