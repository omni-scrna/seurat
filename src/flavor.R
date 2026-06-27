# Map a clustering --flavor name to the integer algorithm id that
# Seurat::FindClusters() accepts natively (1 = original Louvain, 2 = Louvain with
# multilevel refinement, 3 = SLM, 4 = Leiden). Seurat only takes these integers
# (plus "louvain"/"leiden" aliases), so the readable names are ours; the ids are
# Seurat's fixed public enum. --flavor is required: a missing/unknown value errors.
FLAVOR_ALGORITHMS <- c(
  louvain_original              = 1L,
  louvain_multilevel_refinement = 2L,
  slm                           = 3L,
  leiden                        = 4L
)

flavor_algorithm_id <- function(flavor) {
  id <- unname(FLAVOR_ALGORITHMS[flavor])
  if (length(id) != 1L || is.na(id)) {
    stop("--flavor must be one of: ", paste(names(FLAVOR_ALGORITHMS), collapse = ", "))
  }
  id
}
