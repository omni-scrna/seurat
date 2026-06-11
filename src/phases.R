# Phase boundary helper around obkit::logger_emit.
#
# Mirrors the Python rapids-singlecell module's src/phases.py: a small wrapper
# that emits start/end events and lets the callee mutate an attrs env that
# gets flushed on the end-event. The deferred-mutation pattern matters because
# end-event attrs (shapes, nnz, output paths, ...) are only known *after* the
# work runs.
#
# Usage:
#
#   obkit::logger_init(output_dir)   # once, before any phase()
#
#   X <- phase("load", function(attrs) {
#     m <- TENxMatrix(h5_path)
#     attrs$n_genes <- nrow(m)
#     attrs$n_cells <- ncol(m)
#     m
#   })
#
# The end-event is emitted even if `fn` raises, so partial timings still land
# in obkit-events.jsonl for downstream alignment with denet profiling traces.

phase <- function(name, fn) {
  attrs <- new.env(parent = emptyenv())
  logger_emit(name, "start")
  on.exit(
    logger_emit(name, "end", attrs = as.list(attrs)),
    add = TRUE
  )
  fn(attrs)
}
