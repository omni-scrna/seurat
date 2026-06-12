#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(anndataR)
  library(optparse)
  library(data.table)
  library(yaml)
  library(rhdf5)
  library(HDF5Array)
})
cat("OK\n")
