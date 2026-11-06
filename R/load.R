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
#' Wraps a run's output directory; use [load_events()], [load_state()], [trees_path()]
#' and [load_params()] to access its contents.
#'
#' A run goes extinct on the first tick on which no parasite is left in any host or mosquito,
#' hypnozoites included, and no importation is still to come. It ends on that tick: its
#' `events.tsv` and `state.tsv` stop there, and its `params.tsv` ends with the tick as an
#' `EXTINCT_TICK` row, also when the tick is the run's last day.
#'
#' A founder-free run whose imports all fail to establish a parasite never starts, whether by chance
#' or because its settings let no import take (`P_PRIMARY = 0` with a hypnozoite batch mean of 0,
#' or an import window closed throughout, as by `INTERVENTION_EFFICACY = 1`). It ends on the first
#' tick after its last import day, but not before tick 10, or on its last day if that comes first,
#' and its `params.tsv` ends with a `NEVER_STARTED` row. It holds no parasite: its tree sequence
#' has no nodes and its events are at most imports that did not take, and several readers, such as
#' [arg_summary()] and [moi_vs_relapse()], stop on it.
#' @param outdir a run output directory (or an existing `vivax_run`).
#' @return a `vivax_run` object: `outdir`, `params` (from `params.tsv`, as [load_params()] gives
#'   them), `seed`, `extinct_tick`, the tick on which the run went extinct and ended, or `NA`
#'   for a run that did not go extinct, and `never_started`, `TRUE` for a run that never started.
#' @export
load_run <- function(outdir) {
  outdir <- .run_dir(outdir)
  if (!dir.exists(outdir)) stop("no such run directory: ", outdir, call. = FALSE)
  p <- read_params_tsv(file.path(outdir, "params.tsv"))
  seed <- suppressWarnings(as.integer(p[["ACTUAL_SEED"]] %||% p[["SEED"]] %||% NA))
  extinct <- suppressWarnings(as.integer(p[["EXTINCT_TICK"]] %||% NA))
  never <- !is.null(p[["NEVER_STARTED"]])
  p[["EXTINCT_TICK"]] <- NULL   # outcomes, not parameters: see load_params()
  p[["NEVER_STARTED"]] <- NULL
  structure(list(outdir = outdir, params = p, seed = seed, extinct_tick = extinct,
                 never_started = never), class = "vivax_run")
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
trees_path <- function(run) file.path(.run_dir(run), "plasmo.trees")

#' Resolved parameters of a run (from `params.tsv`)
#'
#' An extinct run's `params.tsv` ends with an `EXTINCT_TICK` row, and one that never started with a
#' `NEVER_STARTED` row (see [load_run()]). They are outcomes of the run rather than parameters, so
#' they are left out here, and every run's list has the same names; [load_run()] reads them into
#' the run's `extinct_tick` and `never_started`.
#' @param run a `vivax_run` or a run directory path.
#' @return a named list of the resolved parameter values, with `ACTUAL_SEED`, `SLIM_VERSION`,
#'   `MALARIAFWD_VERSION` and the derived `GENOME_LENGTH`, `N_CHROMOSOMES` and `RECOMB_WITHIN`.
#' @export
load_params <- function(run) {
  p <- if (inherits(run, "vivax_run")) run$params
       else read_params_tsv(file.path(.run_dir(run), "params.tsv"))
  p[["EXTINCT_TICK"]] <- NULL
  p[["NEVER_STARTED"]] <- NULL
  p
}

#' Detect the model generation that produced a run
#' @param run a `vivax_run` or a run directory path.
#' @return `"v2 (meiosis)"` or `"v1 (pre-meiosis)"`.
#' @export
model_version <- function(run) {
  h <- names(data.table::fread(file.path(.run_dir(run), "events.tsv"), sep = "\t", nrows = 0))
  if ("oocyst_id" %in% h) "v2 (meiosis)" else "v1 (pre-meiosis)"
}
