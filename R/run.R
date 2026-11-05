## run.R -- launching the SLiM subprocess (single runs and parallel sweeps).
## Design (see DECISIONS.md D15): SLiM runs as a subprocess via processx (NO shell, so
## Eidos string quotes survive), the model writes its own output files, and reticulate+
## tskit are used only later for IBD. run_grid() parallelises via future/furrr.

#' Run one simulation
#'
#' Launches SLiM with the given parameters and writes its outputs into `outdir`
#' (`events.tsv`, `state.tsv`, `plasmo.trees`, `params.tsv`), plus a provenance snapshot
#' (`model.slim`, `cmdline.txt`, `run.log`). The launch is shell-free (processx), so it
#' handles the SLiM quoting/`mkdir` pitfalls for you.
#'
#' @param params a [vivax_params()] object, or a named list of overrides. Either is checked by
#'   [vivax_params()] again before the run, so an object edited after it was built is checked too.
#' @param seed integer RNG seed. SLiM starts from it (`-s`) and the model sets it as `SEED`,
#'   so the outputs and the tree-sequence provenance both record it.
#' @param outdir output directory (created if needed).
#' @param model path to the SLiM model (defaults to the bundled template).
#' @param verbose if `TRUE`, run the model in VERBOSE mode (full per-event narration).
#'   `NULL` (default) leaves the param as set (quiet).
#' @param timeout seconds before the run is killed (`Inf` = no limit).
#' @param echo if `TRUE`, stream SLiM output live to the console.
#' @return invisibly, a `vivax_run` object (`outdir`, `params`, `seed`).
#' @examples
#' \dontrun{
#' run_sim(vivax_params(N_HOSTS = 200, BITING_RATE = 0.09), seed = 1, outdir = "runs/demo")
#' }
#' @export
run_sim <- function(params = vivax_params(), seed = 1L, outdir, model = slim_model(),
                    verbose = NULL, timeout = Inf, echo = FALSE) {
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) || seed != round(seed) ||
      abs(seed) > .Machine$integer.max)
    stop("`seed` must be a single whole number (|seed| <= ", .Machine$integer.max, ").", call. = FALSE)
  ## Rebuilt every time, so that a set edited after vivax_params() (with $<- or modifyList(),
  ## which keep the class) is checked as well: a misspelt name or an impossible value is refused
  ## here rather than run.
  params <- do.call(vivax_params, as.list(unclass(params)))
  if (!slim_available())
    stop("SLiM not found. Install SLiM (>= 5.2) on PATH or set env var MALARIAFWD_SLIM.", call. = FALSE)
  .check_slim_version()
  if (missing(outdir)) stop("`outdir` is required.", call. = FALSE)
  if (!nzchar(model)) stop("model template not found (is the package installed?).", call. = FALSE)

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)  # SLiM won't create OUTDIR itself
  params$SEED   <- as.integer(seed)
  params$OUTDIR <- normalizePath(outdir, winslash = "/", mustWork = FALSE)
  if (!is.null(verbose)) params$VERBOSE <- isTRUE(verbose)

  pkg_version <- .eidos_quote(as.character(utils::packageVersion("malariafwd")))
  args <- c("-s", as.character(params$SEED), as_slim_args(params),
            "-d", paste0("MALARIAFWD_VERSION=", pkg_version), model)
  ## SLiM reads and writes UTF-8, and vivax_params() has made every string UTF-8. On Unix,
  ## processx converts arguments to the native encoding, which in a C locale (cron, say) turns
  ## UTF-8 text into "<U+00F4>", so they go as the bytes they are; the paths are native already.
  if (.Platform$OS.type == "unix") Encoding(args) <- "bytes"
  file.copy(model, file.path(outdir, "model.slim"), overwrite = TRUE)  # provenance
  writeLines(c(paste0("# ", format(Sys.time(), tz = "UTC")), args), file.path(outdir, "cmdline.txt"),
             useBytes = TRUE)

  res <- processx::run(slim_path(), args, error_on_status = FALSE, timeout = timeout,
                       echo = echo, spinner = FALSE, encoding = "UTF-8")
  writeLines(res$stdout, file.path(outdir, "run.log"), useBytes = TRUE)
  if (!identical(res$status, 0L)) {
    tail_err <- utils::tail(strsplit(paste0(res$stdout, res$stderr), "\n")[[1]], 6)
    stop("SLiM run failed (exit ", res$status, "). Log: ", file.path(outdir, "run.log"), "\n",
         paste(tail_err, collapse = "\n"), call. = FALSE)
  }
  written <- file.exists(file.path(outdir, c("params.tsv", "events.tsv")))
  if (!all(written))
    stop("SLiM exited without an error but did not write ",
         paste(c("params.tsv", "events.tsv")[!written], collapse = " or "), " to ", outdir,
         ". Log: ", file.path(outdir, "run.log"), call. = FALSE)
  invisible(load_run(outdir))  # resolved params from params.tsv is the authoritative record
}

#' Run a grid of parameter sets x seeds in parallel
#'
#' Runs every `(parameter set, seed)` combination, each in its own output directory,
#' parallelised with the current `future` plan. Set `workers` for a quick local
#' multisession; for a cluster, set `future::plan(...)` yourself before calling and
#' leave `workers = NULL` (see DECISIONS.md D15).
#'
#' @param param_sets a [vivax_params()] object, or a (named) list of them.
#' @param seeds integer vector of seeds to run for each set.
#' @param outdir_base parent directory; each run goes in `<base>/<set>_seed<seed>`.
#' @param model path to the SLiM model.
#' @param workers if set, run under `future::plan(multisession, workers)` for this call
#'   (the previous plan is restored afterward). If `NULL`, use the current plan.
#' @param verbose passed to [run_sim()].
#' @param ... passed to [run_sim()].
#' @return a data.frame with one row per run: `set`, `seed`, `outdir`, `ok`, and `error`
#'   (the failure message, or `NA` for a cell that ran). A failed cell does not abort the
#'   sweep -- one bad cell should not lose a grid -- but it warns, and the reason travels
#'   with the result rather than being discarded.
#' @export
run_grid <- function(param_sets, seeds = 1L, outdir_base, model = slim_model(),
                     workers = NULL, verbose = NULL, ...) {
  if (inherits(param_sets, "vivax_params")) param_sets <- list(param_sets)
  if (is.null(names(param_sets)) || any(!nzchar(names(param_sets))))
    names(param_sets) <- sprintf("set%02d", seq_along(param_sets))

  jobs <- expand.grid(set = names(param_sets), seed = as.integer(seeds),
                      stringsAsFactors = FALSE)
  jobs$outdir <- file.path(outdir_base, sprintf("%s_seed%d", jobs$set, jobs$seed))

  if (!is.null(workers)) {
    oldplan <- future::plan(future::multisession, workers = workers)
    on.exit(future::plan(oldplan), add = TRUE)
  }

  out <- furrr::future_pmap(
    list(jobs$set, jobs$seed, jobs$outdir),
    function(set, seed, outdir) {
      ## Keep the error TEXT. `ok = FALSE` alone is unactionable: a cell that failed because
      ## SLiM is not on PATH and one that failed on a bad parameter are indistinguishable, and
      ## a failed cell leaves no directory and no log to read afterwards. The sweep still does
      ## not abort - one bad cell should not lose a grid - but the reason now travels with it.
      err <- NA_character_
      ok <- tryCatch({
        run_sim(param_sets[[set]], seed = seed, outdir = outdir, model = model, verbose = verbose, ...)
        TRUE
      }, error = function(e) { err <<- conditionMessage(e); FALSE })
      data.frame(set = set, seed = seed, outdir = outdir, ok = ok, error = err,
                 stringsAsFactors = FALSE)
    },
    .options = furrr::furrr_options(seed = TRUE)
  )
  res <- do.call(rbind, out)
  ## Say it out loud. A grid is usually launched and left, and a silent half-empty result is
  ## how a sweep gets analysed with a third of its cells missing.
  if (any(!res$ok)) {
    n <- sum(!res$ok)
    warning(n, " of ", nrow(res), " grid cell(s) FAILED and produced no run. First error: ",
            res$error[!res$ok][1], "\nSee the `error` column for the rest.", call. = FALSE)
  }
  res
}

#' @export
print.vivax_run <- function(x, ...) {
  cat("<vivax_run> seed", x$seed, "\n  outdir:", x$outdir, "\n")
  invisible(x)
}
