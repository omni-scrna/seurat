# Vendored from omnibenchmark/obkit @ r/obkit/R/logger.R (MIT).
# Reason: r-obkit has no CRAN/conda-forge/bioconda release — the only
# published build lives on prefix.dev/edge. Vendoring this one small
# logger avoids pulling that extra channel into the solve. Drop the
# vendored copy and switch back to library(obkit) once a conda-forge or
# bioconda build exists. Same wire format — denet alignment is unaffected.

.obkit_state <- new.env(parent = emptyenv())
.obkit_state$log_file <- NULL
.obkit_state$warned_uninit <- FALSE

.iso_ts <- function(t = Sys.time()) {
  ms <- sprintf("%03d", as.integer((as.numeric(t) %% 1) * 1000))
  paste0(format(t, "%Y-%m-%dT%H:%M:%S", tz = "UTC"), ".", ms, "Z")
}

#' Initialize the event logger.
#'
#' Creates `path` if it does not exist and configures subsequent
#' `logger_emit()` calls to append to `<path>/obkit-events.jsonl`.
#'
#' @param path Directory that will hold the event log.
#' @return The resolved log file path, invisibly.
#' @export
logger_init <- function(path) {
  if (!is.character(path) || length(path) != 1L || !nzchar(path)) {
    stop("logger_init(): 'path' must be a non-empty string")
  }
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  log_file <- file.path(normalizePath(path, mustWork = TRUE),
                        "obkit-events.jsonl")
  .obkit_state$log_file <- log_file
  .obkit_state$warned_uninit <- FALSE
  invisible(log_file)
}

#' Emit a lifecycle event.
#'
#' Appends a single JSON-lines record to the configured log file. See
#' `SPEC.md` (obkit-events 0.1) for the wire format.
#'
#' @param event Caller-supplied event name.
#' @param phase Either `"start"` or `"end"`.
#' @param attrs Optional named list of scalar attributes.
#' @return `TRUE` if written, `FALSE` if the logger was not initialized
#'   or the write failed. Returned invisibly.
#' @export
logger_emit <- function(event, phase, attrs = NULL) {
  if (!is.character(event) || length(event) != 1L || !nzchar(event)) {
    stop("logger_emit(): 'event' must be a non-empty string")
  }
  if (!identical(phase, "start") && !identical(phase, "end")) {
    stop("logger_emit(): 'phase' must be \"start\" or \"end\"")
  }
  log_file <- .obkit_state$log_file
  if (is.null(log_file)) {
    if (!isTRUE(.obkit_state$warned_uninit)) {
      warning("obkit: logger_emit() called before logger_init(); events discarded.")
      .obkit_state$warned_uninit <- TRUE
    }
    return(invisible(FALSE))
  }

  rec <- list(
    ts    = .iso_ts(),
    event = event,
    phase = phase,
    pid   = Sys.getpid(),
    host  = as.character(Sys.info()[["nodename"]])
  )
  if (!is.null(attrs)) {
    if (!is.list(attrs) || (length(attrs) > 0 && is.null(names(attrs)))) {
      stop("logger_emit(): 'attrs' must be a named list")
    }
    rec$attrs <- attrs
  }

  line <- jsonlite::toJSON(rec, auto_unbox = TRUE, null = "null")
  ok <- tryCatch({
    con <- file(log_file, open = "a")
    on.exit(close(con), add = TRUE)
    writeLines(line, con, sep = "\n", useBytes = TRUE)
    TRUE
  }, error = function(e) {
    message("obkit: logger_emit() write failed: ", conditionMessage(e))
    FALSE
  })
  invisible(ok)
}
