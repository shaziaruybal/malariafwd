## load.R -- the `vivax_run` object + tidy loaders for a run's outputs.

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a[1])) b else a

# Accept either a `vivax_run` object or a bare directory path.
.run_dir <- function(x) if (inherits(x, "vivax_run")) x$outdir else as.character(x)

#' @keywords internal
read_params_tsv <- function(path) {
  if (!file.exists(path)) return(list())
  d <- data.table::fread(path, sep = "\t", header = TRUE, showProgress = FALSE)
  stats::setNames(as.list(d$value), d$param)
}

#' Load a completed run as a `vivax_run` object
#'
#' Wraps a run's output directory; use [load_events()], [load_state()], [trees_file()]
#' and [run_params()] to access its contents.
#' @param outdir a run output directory (or an existing `vivax_run`).
#' @return a `vivax_run` object.
#' @export
load_run <- function(outdir) {
  outdir <- .run_dir(outdir)
  if (!dir.exists(outdir)) stop("no such run directory: ", outdir, call. = FALSE)
  p <- read_params_tsv(file.path(outdir, "params.tsv"))
  seed <- suppressWarnings(as.integer(p[["ACTUAL_SEED"]] %||% p[["SEED"]] %||% NA))
  structure(list(outdir = outdir, params = p, seed = seed), class = "vivax_run")
}

#' Load the event log of a run
#'
#' Reads `events.tsv`. A schema shim standardises the older `brood_id` column to
#' `lineage_id` and back-fills `oocyst_id`/`detected` so downstream code is
#' generation-agnostic.
#' @param run a `vivax_run` or a run directory path.
#' @return a `data.table` of events.
#' @export
load_events <- function(run) {
  ev <- data.table::fread(file.path(.run_dir(run), "events.tsv"), sep = "\t", showProgress = FALSE)
  if (!"lineage_id" %in% names(ev) && "brood_id" %in% names(ev))
    data.table::setnames(ev, "brood_id", "lineage_id")
  if (!"oocyst_id" %in% names(ev)) ev[, oocyst_id := NA_integer_]
  if (!"detected" %in% names(ev)) ev[, detected := 1L]
  ev[]
}

#' Load the per-tick state table of a run
#' @param run a `vivax_run` or a run directory path.
#' @return a `data.table` of state records.
#' @export
load_state <- function(run) {
  data.table::fread(file.path(.run_dir(run), "state.tsv"), sep = "\t", showProgress = FALSE)[]
}

#' Path to a run's tree-sequence file
#' @param run a `vivax_run` or a run directory path.
#' @return the path to `plasmo.trees`.
#' @export
trees_file <- function(run) file.path(.run_dir(run), "plasmo.trees")

#' Resolved parameters of a run (from `params.tsv`)
#' @param run a `vivax_run` or a run directory path.
#' @return a named list of resolved parameter values.
#' @export
run_params <- function(run) {
  if (inherits(run, "vivax_run")) run$params
  else read_params_tsv(file.path(.run_dir(run), "params.tsv"))
}

#' Detect the model generation that produced a run
#' @param run a `vivax_run` or a run directory path.
#' @return `"v2 (meiosis)"` or `"v1 (pre-meiosis)"`.
#' @export
model_generation <- function(run) {
  h <- names(data.table::fread(file.path(.run_dir(run), "events.tsv"), sep = "\t", nrows = 0))
  if ("oocyst_id" %in% h) "v2 (meiosis)" else "v1 (pre-meiosis)"
}
