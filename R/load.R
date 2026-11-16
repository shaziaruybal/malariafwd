## load.R -- the `vivax_run` object + tidy loaders for a run's outputs.

# Accept either a `vivax_run` object or a bare directory path; a path to a directory holding a run
# that is not whole warns, once a session (D6).
.run_dir <- function(x) {
  if (inherits(x, "vivax_run")) return(x$outdir)
  if (is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)) {
    .warn_if_incomplete(x)
    return(x)
  }
  stop("`run` must be a vivax_run (from run_sim() or load_run()) or the path of a run directory, got ",
       .what_table(x), ".", call. = FALSE)
}

## What a table `x` is, for the messages that refuse it.
.what_table <- function(x) {
  if (is.data.frame(x) && "event_type" %in% names(x)) return("an events table (from load_events())")
  if (is.data.frame(x) && "entity_type" %in% names(x)) return("a state table (from load_state())")
  paste("an object of class", class(x)[1])
}

## Warnings given once per directory per session, keyed by class and path.
.warned <- new.env(parent = emptyenv())

.warn_once <- function(dir, class, msg) {
  key <- paste(class, normalizePath(dir, winslash = "/", mustWork = FALSE))
  if (!is.null(.warned[[key]])) return(invisible(FALSE))
  .warned[[key]] <- TRUE
  warning(structure(class = c(class, "warning", "condition"), list(message = msg, call = NULL)))
  invisible(TRUE)
}

## A directory with a run's outputs but no whole run (.run_complete()).
.warn_if_incomplete <- function(dir) {
  if (!any(file.exists(file.path(dir, c("params.tsv", "events.tsv"))))) return(invisible(FALSE))
  rc <- .run_complete(dir)
  if (rc$complete) return(invisible(FALSE))
  .warn_once(dir, "malariafwd_incomplete_run",
             paste(dir, rc$reason, "This is said once per directory per session."))
}

## A params.tsv that does not hold this version's parameters (r from .read_params_tsv()).
.warn_params_mismatch <- function(r, dir) {
  if (!length(r$missing) && !length(r$dropped)) return(invisible(FALSE))
  one <- function(x, a, b) if (length(x) == 1L) a else b
  why <- c(if (length(r$missing))
             paste0("it lacks ", paste(r$missing, collapse = ", "), " (",
                    one(r$missing, "its default is", "their defaults are"), " used)"),
           if (length(r$dropped))
             paste0(paste(r$dropped, collapse = ", "),
                    one(r$dropped, " is not a parameter", " are not parameters"), " (left out)"))
  .warn_once(dir, "malariafwd_params_mismatch",
             paste0("params.tsv in ", dir, " does not match this version's parameters: ",
                    paste(why, collapse = ", and "), ". This is said once per directory per session."))
}

## A run object whose parameter set has this version's names is returned as it is; any other is
## read again from its directory.
.has_current_set <- function(run)
  inherits(run, "vivax_run") && inherits(run$params, "vivax_params") &&
    identical(names(run$params), names(.malariafwd_defaults))

## An events table, from a table or from a run (a vivax_run or its directory).
.as_events <- function(x, arg = "events") {
  if (inherits(x, "vivax_run") || (is.character(x) && length(x) == 1L)) return(load_events(x))
  if (is.data.frame(x) && "event_type" %in% names(x)) return(data.table::as.data.table(x))
  stop("`", arg, "` must be an events table from load_events(), a vivax_run or a run directory",
       if (is.data.frame(x) && "entity_type" %in% names(x)) "; this is a state table (from load_state())", ".",
       call. = FALSE)
}

## A state table, from a run or its directory, or as given.
.as_state <- function(x, arg = "state") {
  if (inherits(x, "vivax_run") || (is.character(x) && length(x) == 1L)) return(load_state(x))
  if (is.data.frame(x) && "entity_type" %in% names(x)) return(data.table::as.data.table(x))
  stop("`", arg, "` must be a state table from load_state(), a vivax_run or a run directory",
       if (is.data.frame(x) && "event_type" %in% names(x)) "; this is an events table (from load_events())", ".",
       call. = FALSE)
}

## The table `fun()` returns: `x` itself when it is such a table, with the columns `cols` a plot
## reads, or fun() of `x`, a run or its directory.
.as_result <- function(x, fun, cols, arg = "x") {
  if (!is.data.frame(x)) return(fun(x))
  miss <- setdiff(cols, names(x))
  if (length(miss))
    stop("`", arg, "` must be a run, or the table ", deparse(substitute(fun)), "() returns; it has no ",
         paste0("`", miss, "`", collapse = " or "), ".", call. = FALSE)
  data.table::as.data.table(x)
}

## A day argument `name`: a single whole number of days, 0 or more.
.check_day <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) || x < 0 || x != round(x))
    stop("`", name, "` must be a single whole number of days, 0 or more.", call. = FALSE)
  invisible(as.integer(x))
}

## A function's flag argument `name`: a single TRUE or FALSE.
.check_flag <- function(x, name) {
  if (!isTRUE(x) && !isFALSE(x)) stop("`", name, "` must be TRUE or FALSE.", call. = FALSE)
  invisible(x)
}

## A run's length, SIM_DAYS, from which every function that summarises a window of the run takes its
## end, so that the days after an extinction count as days with no parasite (D49).
.sim_days <- function(run, fun = "this function") {
  if (is.data.frame(run))
    stop(fun, "() needs the run for its length, SIM_DAYS: pass the vivax_run or its directory, got ",
         .what_table(run), ".", call. = FALSE)
  n <- .params_of(run)$SIM_DAYS
  if (length(n) != 1L || is.na(n) || n < 1 || n != round(n))
    stop("the run's params.tsv gives no SIM_DAYS: ", .run_dir(run), call. = FALSE)
  as.integer(n)
}

## A run's parameter set, for package code that needs a few of its values: what load_params()
## gives, without the warning.
.params_of <- function(run) {
  if (.has_current_set(run)) return(run$params)
  .read_params_tsv(file.path(.run_dir(run), "params.tsv"))$params
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
#' A founder-free run whose imports all fail to establish a parasite never starts, whether by
#' chance or because its settings let no import take (`P_PRIMARY = 0` with a hypnozoite batch mean
#' of 0, or an `IMPORT_RATE_PER_TICK` window whose rate is 0 throughout, as under
#' `IMPORT_RATE_BOOST = 0` from a trigger before the window). It ends on the first tick after its
#' last import day, but not before tick 10, or on its last day if that comes first, and its
#' `params.tsv` ends with a `NEVER_STARTED` row. It holds no parasite: its tree sequence has no
#' nodes and its events are at most imports that did not take. Every function that reads a run
#' takes such a run: its
#' tables have their usual columns and no rows, its counts are 0, and a quantity that needs a
#' parasite is `NA` ([arg_summary()] reports no samples).
#'
#' The functions that summarise a window of a run take its length from `SIM_DAYS`, so the days
#' after an extinction count as days with no parasite.
#'
#' A directory without `params.tsv` or `events.tsv` is not a run and is refused. One whose `DONE`
#' file is missing, or whose files have changed since the run finished (see [run_sim()]), loads
#' with `complete = FALSE` and a warning of class `malariafwd_incomplete_run`, given once per
#' directory per session by every function that reads the directory. A run started with `slim`
#' directly
#' has no `DONE`.
#' `suppressWarnings(..., classes = "malariafwd_incomplete_run")` silences it.
#' @param outdir a run output directory, or a `vivax_run`, whose directory is read again.
#' @return a `vivax_run` object, a list of:
#'   * `outdir`, the run directory (character);
#'   * `params`, the parameter set, as [load_params()] gives it;
#'   * `seed`, the seed SLiM ran with (integer; a double beyond the 32-bit range, which only a run
#'     started with `slim` directly, with `-d SEED`, can have);
#'   * `extinct_tick`, the tick on which the run went extinct and ended, or `NA` (integer);
#'   * `never_started`, `TRUE` for a run that never started (logical);
#'   * `slim_version` and `malariafwd_version` (character; `malariafwd_version` is `"unknown"`
#'     for a run started with `slim` directly);
#'   * `genome_length`, the simulated genome in bp, and `n_chromosomes` (integers), and
#'     `recomb_within`, the crossover rate within a chromosome per bp per meiosis (double);
#'   * `complete`, `TRUE` when the directory holds the whole run its `DONE` file describes
#'     (logical).
#'
#'   A field the run's `params.tsv` does not record is `NA`.
#' @examples
#' run <- load_run(example_run()$outdir)   # a run directory, read again
#' run$seed
#' run$params$N_HOSTS
#' @export
load_run <- function(outdir) {
  dir <- .run_dir(outdir)
  if (!dir.exists(dir)) stop("no such run directory: ", dir, call. = FALSE)
  miss <- c("params.tsv", "events.tsv")[!file.exists(file.path(dir, c("params.tsv", "events.tsv")))]
  if (length(miss))
    stop(dir, " is not a malariafwd run: it has no ", paste(miss, collapse = " or "), ".", call. = FALSE)
  r <- .read_params_tsv(file.path(dir, "params.tsv"))
  .warn_params_mismatch(r, dir)
  rec <- r$record
  or <- function(x, y) if (is.null(x)) y else x
  structure(list(outdir = dir, params = r$params,
                 seed = or(rec$ACTUAL_SEED, or(rec$SEED, NA_integer_)),
                 extinct_tick = or(rec$EXTINCT_TICK, NA_integer_),
                 never_started = !is.null(rec$NEVER_STARTED),
                 slim_version = or(rec$SLIM_VERSION, NA_character_),
                 malariafwd_version = or(rec$MALARIAFWD_VERSION, NA_character_),
                 genome_length = or(rec$GENOME_LENGTH, NA_integer_),
                 n_chromosomes = or(rec$N_CHROMOSOMES, NA_integer_),
                 recomb_within = or(rec$RECOMB_WITHIN, NA_real_),
                 complete = .run_complete(dir)$complete),
            class = "vivax_run")
}

#' Load the event log of a run
#'
#' Reads `events.tsv`. Every run gives the same 21 columns, in the order and classes the model
#' writes them ([events.tsv][malariafwd-outputs] describes them), whatever the run holds: a run with no event gives
#' a table with no rows. A column the file lacks is `NA`, a column it adds follows the 21, and a
#' value that does not fit its column stops with the column named.
#' @param run a run: a `vivax_run` (from [run_sim()] or [load_run()]) or the path of a run directory.
#' @return a `data.table` with one row per event, in the order they were logged.
#' @examples
#' ev <- load_events(example_run())
#' ev[, .N, by = event_type]
#' @export
load_events <- function(run) .read_output(run, "events.tsv", .EVENTS_SCHEMA)

#' Load the per-tick state table of a run
#'
#' Reads `state.tsv`, by default its host rows. A mosquito's row counts the lineages it carries, so a
#' prevalence or MOI taken over every row mixes the two populations. Every run gives the same 12
#' columns, in the order and classes the model writes them ([state.tsv][malariafwd-outputs]
#' describes them),
#' whatever the run holds; a column the file lacks is `NA`.
#' @inheritParams load_events
#' @param entity the rows to give: `"human"` (the hosts), `"mosquito"`, or `"all"`, both in the
#'   order of the file.
#' @return a `data.table` with one row per host (or mosquito, or both) per day, from day 1 to the
#'   run's last day.
#' @examples
#' st <- load_state(example_run())
#' st[tick == 400]   # every host on the last day
#' @export
load_state <- function(run, entity = c("human", "mosquito", "all")) {
  entity <- match.arg(entity)
  d <- .read_output(run, "state.tsv", .STATE_SCHEMA)
  if (entity != "all") d <- d[entity_type == entity]
  d[]
}

#' Path to a run's tree-sequence file
#' @inheritParams load_events
#' @return the path to `plasmo.trees`.
#' @seealso [plasmo.trees][malariafwd-outputs]
#' @examples
#' trees_path(example_run())
#' file.exists(trees_path(example_run()))
#' @export
trees_path <- function(run) file.path(.run_dir(run), "plasmo.trees")

#' The parameter set a run was made from
#'
#' Reads `params.tsv` ([params.tsv][malariafwd-outputs]) into a [vivax_params()] object holding
#' every parameter, in order and of its type: `run_sim(load_params(run), seed = run$seed,
#' outdir = tempfile())` gives the same run again. Values are not checked here ([run_sim()] checks
#' them). The seed, the versions, the derived genome and the outcome rows are on the [load_run()]
#' object instead.
#'
#' A file from another version of the package can lack some parameters, which take their
#' defaults, or hold rows that are not parameters, which are left out. Either warns, once per
#' directory per session, with class `malariafwd_params_mismatch`.
#' @inheritParams load_events
#' @return a [vivax_params()] object.
#' @examples
#' p <- load_params(example_run())
#' p
#' @examplesIf slim_available()
#' run <- run_sim(p, seed = 3, outdir = tempfile("seed3"))   # the same settings, another seed
#' @export
load_params <- function(run) {
  if (.has_current_set(run)) return(run$params)
  dir <- .run_dir(run)
  r <- .read_params_tsv(file.path(dir, "params.tsv"))
  .warn_params_mismatch(r, dir)
  r$params
}
