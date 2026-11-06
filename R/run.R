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
#' @param timeout seconds before the run is stopped, with an error saying so (`Inf` = no limit).
#' @param echo if `TRUE`, stream SLiM output live to the console.
#' @return invisibly, a `vivax_run` object (`outdir`, `params`, `seed`, `extinct_tick`,
#'   `never_started`; see [load_run()]). A run that goes extinct, with no parasite left and no
#'   importation still to come, ends on that tick and gives a message; so does a founder-free run
#'   whose imports never established a parasite.
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
  params <- do.call(vivax_params, as.list(unclass(params)), quote = TRUE)
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

  ## SLiM's output is UTF-8. processx's own buffer translates it to the native encoding and then
  ## cuts it at the UTF-8 byte count, which in a C locale loses the end of the log, so the chunks
  ## are kept here as they arrive (and echoed here when asked).
  out <- character(); err <- character()
  res <- processx::run(slim_path(), args, error_on_status = FALSE, timeout = timeout, echo = FALSE,
                       spinner = FALSE, encoding = "UTF-8",
                       stdout_callback = function(x, proc) { out[[length(out) + 1L]] <<- x; if (echo) cat(x) },
                       stderr_callback = function(x, proc) { err[[length(err) + 1L]] <<- x; if (echo) cat(x, file = stderr()) })
  out <- paste(out, collapse = ""); err <- paste(err, collapse = "")
  writeLines(out, file.path(outdir, "run.log"), useBytes = TRUE)
  if (isTRUE(res$timeout))
    stop("SLiM run timed out after ", if (inherits(timeout, "difftime")) format(timeout) else paste(timeout, "s"),
         " and was stopped. Log: ", file.path(outdir, "run.log"), call. = FALSE)
  if (!identical(res$status, 0L)) {
    tail_err <- utils::tail(strsplit(paste0(out, err), "\n")[[1]], 6)
    stop("SLiM run failed (exit ", res$status, "). Log: ", file.path(outdir, "run.log"), "\n",
         paste(tail_err, collapse = "\n"), call. = FALSE)
  }
  written <- file.exists(file.path(outdir, c("params.tsv", "events.tsv")))
  if (!all(written))
    stop("SLiM exited without an error but did not write ",
         paste(c("params.tsv", "events.tsv")[!written], collapse = " or "), " to ", outdir,
         ". Log: ", file.path(outdir, "run.log"), call. = FALSE)
  run <- load_run(outdir)  # resolved params from params.tsv is the authoritative record
  if (!is.na(run$extinct_tick))
    message("The run went extinct on tick ", run$extinct_tick, " of ", .slim_num(params$SIM_DAYS),
            ": no parasite was left in any host or mosquito, hypnozoites included, and no ",
            "importation was still to come, so its outputs end there.")
  if (isTRUE(run$never_started))
    message("The run never started and ended on tick ",
            read_params_tsv(file.path(outdir, "params.tsv"))[["NEVER_STARTED"]], " of ",
            .slim_num(params$SIM_DAYS), ": no import established a parasite and none was still to ",
            "come, so its outputs hold no parasite.")
  invisible(run)
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
#' @return a data.frame with one row per run: `set`, `seed`, `outdir`, `ok`, `error` (the failure
#'   message, or `NA` for a cell that ran), `extinct_tick` (the tick on which a run went extinct
#'   and ended, `NA` for a run that did not go extinct) and `never_started` (`TRUE` for a
#'   founder-free run whose imports never established a parasite, `FALSE` for a run that started;
#'   see [load_run()]). Both outcome columns are `NA` for a failed cell. A cell that never started
#'   is `ok`: not starting is an outcome of the grid, which `never_started` records. Its run holds
#'   no parasite, though, so give only the cells with `ok & !never_started` to the readers that
#'   need one, such as [arg_summary()] and [arg_diversity()]. The cells that went extinct, and
#'   those that never started, are counted in one message each for the grid, in place of
#'   [run_sim()]'s message for each. A failed cell does not abort the sweep -- one bad cell should
#'   not lose a grid -- but it warns, and the reason travels with the result rather than being
#'   discarded. A warning given in the cells (such as a parameter set's own warnings, which
#'   [run_sim()] repeats for each run) is given once, after the FAILED summary, with the number of
#'   cells it came from. If the grid stops part-way with an error, that error's message lists the
#'   warnings of the cells whose results had reached this session; with workers, a cell's warnings
#'   travel with its chunk's result, so those of a chunk lost with its worker are lost. Under
#'   `options(warn = 2)` a warning in a cell fails that cell, in any plan and also inside
#'   `suppressWarnings()`, and the FAILED summary is then an error, so nothing is returned; the
#'   cells that ran are on disk.
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

  ## A set that warns when it is checked (biting above 1 a day, say) warns again in each of its
  ## cells, and R keeps only the first 50 warnings, which would crowd out the FAILED summary below.
  ## So each distinct warning is kept back (its first condition, so its class survives) and given
  ## once, after the summary: a caller's handler that catches the summary then finds nothing left
  ## to give. A grid stopped part-way by an error passes them on in that error's message, and an
  ## interrupted one gives them from on.exit(). Under warn = 2 a cell's warning is made its error
  ## in the cell itself, so that it fails the cell whatever the plan (workers do not inherit warn)
  ## and whatever the version of future.
  warn_level <- getOption("warn")
  msgs <- character(); conds <- list(); counts <- integer()
  keep <- function(w) {
    if (warn_level >= 2) return()
    m <- conditionMessage(w); i <- match(m, msgs)
    if (is.na(i)) { msgs <<- c(msgs, m); conds[[length(conds) + 1L]] <<- w; counts <<- c(counts, 1L) }
    else counts[i] <<- counts[i] + 1L
    invokeRestart("muffleWarning")
  }
  counted <- function() paste0(msgs, ifelse(counts > 1, paste0(" [", counts, " cells]"), ""))
  give <- function() {
    m <- counted(); cs <- conds
    msgs <<- character(); conds <<- list(); counts <<- integer()
    for (i in seq_along(m)) { w <- cs[[i]]; w$message <- m[i]; warning(w) }
  }
  on.exit(give(), add = TRUE)
  out <- tryCatch(withCallingHandlers(furrr::future_pmap(
    list(jobs$set, jobs$seed, jobs$outdir),
    function(set, seed, outdir) {
      ## Keep the error TEXT. `ok = FALSE` alone is unactionable: a cell that failed because
      ## SLiM is not on PATH and one that failed on a bad parameter are indistinguishable, and
      ## a failed cell leaves no directory and no log to read afterwards. The sweep still does
      ## not abort - one bad cell should not lose a grid - but the reason now travels with it.
      err <- NA_character_; extinct <- NA_integer_; never <- NA
      run <- function() run_sim(param_sets[[set]], seed = seed, outdir = outdir, model = model,
                                verbose = verbose, ...)
      ok <- tryCatch({
        ## the extinction message is summed up once for the grid, below
        r <- suppressMessages(if (warn_level < 2) run() else withCallingHandlers(run(), warning =
          function(w) stop("(converted from warning) ", conditionMessage(w), call. = FALSE)))
        extinct <- r$extinct_tick; never <- isTRUE(r$never_started)
        TRUE
      }, error = function(e) { err <<- conditionMessage(e); FALSE })
      data.frame(set = set, seed = seed, outdir = outdir, ok = ok, error = err,
                 extinct_tick = extinct, never_started = never, stringsAsFactors = FALSE)
    },
    .options = furrr::furrr_options(seed = TRUE)
  ), warning = keep), error = function(e) {
    if (length(msgs)) {
      e$message <- paste0(conditionMessage(e), "\nThe cells that had reported gave these warnings:\n",
                          paste0("  ", counted(), collapse = "\n"))
      msgs <<- character(); conds <<- list(); counts <<- integer()
    }
    stop(e)
  })
  res <- do.call(rbind, out)
  ## Say it out loud. A grid is usually launched and left, and a silent half-empty result is
  ## how a sweep gets analysed with a third of its cells missing. Under warn = 2 the summary is an
  ## error, and nothing is returned.
  held <- list(msgs, conds, counts); msgs <- character(); conds <- list(); counts <- integer()
  if (any(!res$ok)) {
    n <- sum(!res$ok)
    warning(n, " of ", nrow(res), " grid cell(s) FAILED and produced no run. First error: ",
            res$error[!res$ok][1], if (warn_level < 2) "\nSee the `error` column for the rest.",
            call. = FALSE)
  }
  msgs <- held[[1]]; conds <- held[[2]]; counts <- held[[3]]
  give()
  if (any(!is.na(res$extinct_tick)))
    message(sum(!is.na(res$extinct_tick)), " of ", nrow(res), " grid cell(s) went extinct; see the ",
            "`extinct_tick` column.")
  if (any(res$never_started %in% TRUE))
    message(sum(res$never_started %in% TRUE), " of ", nrow(res), " grid cell(s) never started: no ",
            "import established a parasite; see the `never_started` column.")
  res
}

#' @export
print.vivax_run <- function(x, ...) {
  cat("<vivax_run> seed", x$seed, "\n  outdir:", x$outdir, "\n")
  if (!is.null(x$extinct_tick) && !is.na(x$extinct_tick))
    cat("  went extinct on tick", x$extinct_tick, "of", x$params$SIM_DAYS %||% NA, "\n")
  if (isTRUE(x$never_started))
    cat("  never started: no import established a parasite\n")
  invisible(x)
}
