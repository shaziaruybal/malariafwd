## run.R -- launching the SLiM subprocess (single runs and parallel sweeps).
## Design (see DECISIONS.md D15): SLiM runs as a subprocess via processx (NO shell, so
## Eidos string quotes survive), the model writes its own output files, and reticulate+
## tskit are used only later for IBD. run_grid() parallelises via future/furrr.

#' Run one simulation
#'
#' Launches SLiM with the given parameters and writes the run into `outdir`: `events.tsv`,
#' `state.tsv`, `plasmo.trees` and `params.tsv`, with `model.slim`, `cmdline.txt`, `run.log` and
#' `DONE` (see *The run directory*). `params.tsv` records every parameter, floats to the 15
#' significant digits passed to SLiM. The launch is shell-free (processx), so it handles the SLiM
#' quoting pitfalls for you.
#'
#' @section The run directory:
#' The run is written into a new directory beside `outdir` and moved to `outdir` only when it has
#' finished, so `outdir` holds a whole run or what it held before, never a mix. `DONE`, written
#' last, records when the run finished, how long it took, its seed, a key of its parameters and
#' seed, the md5 of the model, the versions, and the size of every file; [load_run()] reports a
#' directory whose `DONE` is missing or whose files have changed since. SLiM runs in the run's
#' directory (`OUTDIR` is `'.'`) on the copy of the model it holds, `model.slim`, so the directory
#' is self-contained: `cmdline.txt` gives the arguments, and running `slim` with them in a
#' directory holding `model.slim` repeats the run. A custom model that reads other files must use
#' absolute paths. A run that fails, times out or is interrupted leaves its files and `run.log` in
#' `<outdir>.failed`; if R itself is killed, a `<outdir>.partial-*` directory is left beside
#' `outdir`, to delete by hand.
#'
#' @param params a [vivax_params()] object, or a named list of overrides. Either is checked by
#'   [vivax_params()] again before the run, so an object edited after it was built is checked too.
#' @param seed integer RNG seed. SLiM starts from it (`-s`) and the model sets it as `SEED`,
#'   so the outputs and the tree-sequence provenance both record it.
#' @param outdir the run's directory, created when the run has finished (its parent is created if
#'   needed). A directory that already holds a run is not replaced unless `overwrite = TRUE`, and a
#'   directory holding other files is never used.
#' @param model path to the SLiM model (defaults to the bundled template).
#' @param verbose if `TRUE`, run the model in VERBOSE mode (full per-event narration).
#'   `NULL` (default) leaves the param as set (quiet).
#' @param timeout seconds before the run is stopped, with an error saying so (`Inf` = no limit).
#' @param echo if `TRUE`, stream SLiM output live to the console.
#' @param overwrite if `TRUE`, a run already in `outdir` is replaced, once the new run has
#'   finished. Only run files are ever deleted.
#' @param skip_existing if `TRUE` and `outdir` already holds this run (the same parameters, seed
#'   and model), it is returned without running SLiM.
#' @return invisibly, a `vivax_run` object: the run directory, its parameter set, seed, outcome,
#'   versions, genome and `complete` (see [load_run()]). A run that goes extinct, with no parasite
#'   left and no importation still to come, ends on that tick and gives a message; so does a
#'   founder-free run whose imports never established a parasite.
#' @examplesIf slim_available()
#' run <- run_sim(vivax_params(SIM_DAYS = 100), seed = 1, outdir = tempfile("demo"))
#' run
#' summarise_run(run)
#' @export
run_sim <- function(params = vivax_params(), seed = 1L, outdir, model = slim_model(),
                    verbose = NULL, timeout = Inf, echo = FALSE, overwrite = FALSE,
                    skip_existing = FALSE) {
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) || seed != round(seed) ||
      abs(seed) > .Machine$integer.max)
    stop("`seed` must be a single whole number (|seed| <= ", .Machine$integer.max, ").", call. = FALSE)
  ## Rebuilt every time, so that a set edited after vivax_params() (with $<- or modifyList(),
  ## which keep the class) is checked as well: a misspelt name or an impossible value is refused
  ## here rather than run.
  if (inherits(params, "vivax_params")) .refuse_unknown(names(params), built = TRUE)
  params <- do.call(vivax_params, as.list(unclass(params)), quote = TRUE)
  if (missing(outdir)) stop("`outdir` is required.", call. = FALSE)
  if (!is.character(outdir) || length(outdir) != 1L || is.na(outdir) || !nzchar(outdir))
    stop("`outdir` must be one directory path.", call. = FALSE)
  outdir <- path.expand(outdir)
  if (basename(outdir) %in% c("", ".", ".."))
    stop("`outdir` must name a directory of its own, not \"", outdir, "\".", call. = FALSE)
  if (file.exists(outdir) && !dir.exists(outdir)) stop("`outdir` is a file: ", outdir, call. = FALSE)
  for (f in c("overwrite", "skip_existing"))
    if (!isTRUE(get(f)) && !isFALSE(get(f))) stop("`", f, "` must be TRUE or FALSE.", call. = FALSE)
  if (!is.character(model) || length(model) != 1L || is.na(model) || !nzchar(model))
    stop("model template not found (is the package installed?).", call. = FALSE)
  if (!file.exists(model)) stop("model file not found: ", model, call. = FALSE)
  model <- normalizePath(model, winslash = "/")
  if (!is.null(verbose)) params$VERBOSE <- isTRUE(verbose)

  ## what is at outdir already, and what to do with it
  key <- .run_key(params, seed)
  model_md5 <- unname(tools::md5sum(model))
  rs <- .run_state(outdir, key, model_md5)
  if (rs$state == "same" && skip_existing) {
    message(outdir, " already holds this run (same parameters, seed and model; finished ",
            sub("T", " ", sub("Z$", "", rs$done[["finished"]])), " UTC): not run again.")
    return(invisible(load_run(outdir)))
  }
  if (length(rs$foreign))
    .outdir_in_use(outdir, " holds files that are not part of a run (",
                   paste(utils::head(rs$foreign, 5), collapse = ", "), if (length(rs$foreign) > 5) ", ...",
                   "): choose another outdir or move them.")
  if (rs$state == "same" && !overwrite)
    .outdir_in_use(outdir, " already holds this run: skip_existing = TRUE returns it, overwrite = TRUE ",
                   "runs it again.")
  if (rs$state %in% c("other", "incomplete") && !overwrite)
    .outdir_in_use(outdir, " already holds a run of other parameters, seed or model (or one that did ",
                   "not finish): overwrite = TRUE replaces it once the new run has finished, or choose ",
                   "another outdir.")
  replace <- rs$state != "none"

  slim_v <- .require_slim()

  ## The run is written beside outdir, so that the move is a rename on one file system
  dir.create(dirname(outdir), recursive = TRUE, showWarnings = FALSE)
  outdir <- file.path(normalizePath(dirname(outdir), winslash = "/", mustWork = TRUE), basename(outdir))
  partial <- tempfile(paste0(basename(outdir), ".partial-"), tmpdir = dirname(outdir))
  dir.create(partial)
  failed <- paste0(outdir, ".failed")
  out <- character()
  finished <- FALSE; keep_partial <- FALSE
  on.exit(if (!finished && !keep_partial) {
    if (!file.exists(file.path(partial, "run.log")))
      writeLines(paste(out, collapse = ""), file.path(partial, "run.log"), useBytes = TRUE)
    .remove_run_dir(failed)
    if (!dir.exists(failed)) file.rename(partial, failed)
  }, add = TRUE)

  params$SEED   <- as.integer(seed)
  params$OUTDIR <- "."
  pkg_version <- .eidos_quote(as.character(utils::packageVersion("malariafwd")))
  args <- c("-s", as.character(params$SEED), as_slim_args(params),
            "-d", paste0("MALARIAFWD_VERSION=", pkg_version), "model.slim")
  ## SLiM reads and writes UTF-8, and vivax_params() has made every string UTF-8. On Unix,
  ## processx converts arguments to the native encoding, which in a C locale (cron, say) turns
  ## UTF-8 text into "<U+00F4>", so they go as the bytes they are; the paths are native already.
  if (.Platform$OS.type == "unix") Encoding(args) <- "bytes"
  file.copy(model, file.path(partial, "model.slim"))
  writeLines(c(paste0("# ", format(Sys.time(), tz = "UTC"), " UTC"),
               paste("# SLiM was started in the run's directory (OUTDIR is '.'), running its copy of the",
                     "model, model.slim; to repeat the run, run slim with these arguments from a directory",
                     "holding model.slim"), args),
             file.path(partial, "cmdline.txt"), useBytes = TRUE)

  ## SLiM's output is UTF-8. processx's own buffer translates it to the native encoding and then
  ## cuts it at the UTF-8 byte count, which in a C locale loses the end of the log, so the chunks
  ## are kept here as they arrive (and echoed here when asked). Its errors go to the same log, and
  ## supervise and cleanup_tree stop SLiM if R dies.
  t0 <- Sys.time()
  res <- processx::run(slim_path(), args, wd = partial, error_on_status = FALSE, timeout = timeout,
                       echo = FALSE, spinner = FALSE, encoding = "UTF-8", stderr_to_stdout = TRUE,
                       cleanup_tree = TRUE, supervise = TRUE,
                       stdout_callback = function(x, proc) { out[[length(out) + 1L]] <<- x; if (echo) cat(x) })
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  log_text <- paste(out, collapse = "")
  writeLines(log_text, file.path(partial, "run.log"), useBytes = TRUE)
  if (isTRUE(res$timeout))
    stop("SLiM run timed out after ", if (inherits(timeout, "difftime")) format(timeout) else paste(timeout, "s"),
         " and was stopped. Its files and log are in ", failed, call. = FALSE)
  if (!identical(res$status, 0L)) {
    tail_err <- utils::tail(strsplit(log_text, "\n")[[1]], 6)
    stop("SLiM run failed (exit ", res$status, "). Log: ", file.path(failed, "run.log"), "\n",
         paste(tail_err, collapse = "\n"), call. = FALSE)
  }
  outputs <- c("params.tsv", "events.tsv", "state.tsv", "plasmo.trees")
  written <- file.exists(file.path(partial, outputs))
  if (!all(written))
    stop("SLiM exited without an error but did not write ", paste(outputs[!written], collapse = " or "),
         ". Log: ", file.path(failed, "run.log"), call. = FALSE)
  .write_done(partial, seed = as.integer(seed), run_key = key, model_md5 = model_md5,
              elapsed_s = elapsed, slim_version = slim_v)

  ## Place the run. A run being replaced is moved aside first and deleted only once the new one is
  ## in place; anything this call may not replace that appeared during the run stops it, and the
  ## finished run stays where it is.
  now <- .run_state(outdir, key, model_md5)
  if (length(now$foreign) || (now$state != "none" && !replace)) {
    keep_partial <- TRUE
    .outdir_in_use(outdir, " appeared while this run was running; the finished run is in ", partial)
  }
  old <- NULL
  if (dir.exists(outdir)) {
    if (now$state == "none") unlink(outdir, recursive = TRUE)
    else {
      old <- tempfile(paste0(basename(outdir), ".old-"), tmpdir = dirname(outdir))
      if (!file.rename(outdir, old)) { keep_partial <- TRUE
        .outdir_in_use("the run in ", outdir, " could not be moved aside; the finished run is in ", partial) }
    }
  }
  if (!file.rename(partial, outdir)) {
    if (!is.null(old)) file.rename(old, outdir)
    keep_partial <- TRUE
    stop("the finished run could not be moved to ", outdir, "; it is in ", partial, call. = FALSE)
  }
  finished <- TRUE
  if (!is.null(old) && !isTRUE(.remove_run_dir(old)))
    warning("the run this one replaced could not be deleted: ", old, call. = FALSE)
  .remove_run_dir(failed)

  run <- load_run(outdir)  # resolved params from params.tsv is the authoritative record
  if (!is.na(run$extinct_tick))
    message("The run went extinct on tick ", run$extinct_tick, " of ", .slim_num(params$SIM_DAYS),
            ": no parasite was left in any host or mosquito, hypnozoites included, and no ",
            "importation was still to come, so its outputs end there.")
  if (isTRUE(run$never_started))
    message("The run never started and ended on tick ",
            .read_params_tsv(file.path(outdir, "params.tsv"))$record$NEVER_STARTED, " of ",
            .slim_num(params$SIM_DAYS), ": no import established a parasite and none was still to ",
            "come, so its outputs hold no parasite.")
  invisible(run)
}

#' @export
print.vivax_run <- function(x, ...) {
  p <- x$params
  typed <- inherits(p, "vivax_params")
  cat("<vivax_run> seed", format(x$seed), "\n  outdir:", x$outdir, "\n")
  if (typed)
    cat("  ", .slim_num(p$SIM_DAYS), " days, ", .slim_num(p$N_HOSTS), " hosts, ", .slim_num(p$N_MOSQ),
        " mosquitoes\n", sep = "")
  if (!is.null(x$extinct_tick) && !is.na(x$extinct_tick))
    cat("  went extinct on tick", x$extinct_tick, "of", format(p$SIM_DAYS), "\n")
  if (isTRUE(x$never_started))
    cat("  never started: no import established a parasite\n")
  if (typed) .print_non_default(p)
  if (!is.null(x$slim_version))
    cat("  SLiM ", x$slim_version, ", malariafwd ", x$malariafwd_version, "\n", sep = "")
  if (isFALSE(x$complete)) cat("  not marked complete: the directory", .run_complete(x$outdir)$reason, "\n")
  invisible(x)
}

## ---- the run directory ---------------------------------------------------------------------
## A run is written into a sibling directory and moved into place when it has finished, with DONE
## written last (D58). These are the files a run holds; the ignored ones are what file managers
## leave behind.
.run_files <- c("DONE", "params.tsv", "events.tsv", "state.tsv", "plasmo.trees", "model.slim",
                "cmdline.txt", "run.log")
.ignored_files <- c(".DS_Store", "Thumbs.db", "desktop.ini")

## What identifies a run: the md5 of the sorted NAME=value tokens SLiM is given, without OUTDIR
## and VERBOSE (narration only), with SEED. The tokens are what SLiM receives, so the key depends
## neither on params.tsv's precision nor on the order of the parameters.
.run_key <- function(params, seed) {
  a <- as_slim_args(params[setdiff(names(params), c("OUTDIR", "VERBOSE", "SEED"))])
  ## sorted and hashed as their UTF-8 bytes, so that the key is the same in every locale
  tok <- enc2utf8(c(a[c(FALSE, TRUE)], paste0("SEED=", as.integer(seed))))
  Encoding(tok) <- "bytes"
  tok <- sort(tok, method = "radix")
  f <- tempfile(); on.exit(unlink(f))
  writeBin(unlist(lapply(tok, function(t) c(charToRaw(t), as.raw(10L)))), f)
  unname(tools::md5sum(f))
}

## DONE: key<TAB>value rows, then one file:<name> row per file the run holds, with its size.
.write_done <- function(dir, seed, run_key, model_md5, elapsed_s, slim_version) {
  files <- sort(setdiff(list.files(dir, all.files = TRUE, no.. = TRUE), c("DONE", .ignored_files)),
                method = "radix")
  rows <- c(finished = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
            elapsed_s = sprintf("%.3f", elapsed_s), seed = as.character(seed), run_key = run_key,
            model_md5 = model_md5, malariafwd_version = as.character(utils::packageVersion("malariafwd")),
            slim_version = if (is.null(slim_version) || is.na(slim_version)) "NA" else slim_version,
            stats::setNames(format(file.size(file.path(dir, files)), scientific = FALSE, trim = TRUE),
                            paste0("file:", files)))
  writeLines(c("key\tvalue", paste0(names(rows), "\t", rows)), file.path(dir, "DONE"), useBytes = TRUE)
}

## A run's DONE as a named character vector, or NULL when it has none or it cannot be read.
.read_done <- function(dir) {
  f <- file.path(dir, "DONE")
  if (!file.exists(f)) return(NULL)
  l <- readLines(f, warn = FALSE, encoding = "UTF-8")
  if (!length(l) || l[1] != "key\tvalue") return(NULL)
  l <- l[-1][grepl("\t", l[-1], fixed = TRUE)]
  tab <- regexpr("\t", l, fixed = TRUE)
  stats::setNames(substring(l, tab + 1), substr(l, 1, tab - 1))
}

## Whether a directory holds a whole run: DONE, and every file it lists, of the size it gives.
.run_complete <- function(dir) {
  done <- .read_done(dir)
  if (is.null(done))
    return(list(complete = FALSE, reason = paste(
      "has no DONE file: run_sim() did not finish this run, or it was made another way (such as",
      "slim -d, or an earlier build), so its outputs may stop early or come from different runs.")))
  listed <- sub("^file:", "", grep("^file:", names(done), value = TRUE))
  path <- file.path(dir, listed)
  gone <- listed[!file.exists(path)]
  resized <- setdiff(listed[file.exists(path) & file.size(path) != suppressWarnings(as.numeric(done[paste0("file:", listed)]))], gone)
  if (length(gone) || length(resized))
    return(list(complete = FALSE, reason = paste0(
      "has changed since its run finished: ",
      paste(c(if (length(gone)) paste(paste(gone, collapse = ", "), if (length(gone) == 1L) "is" else "are", "missing"),
              if (length(resized)) paste(paste(resized, collapse = ", "), if (length(resized) == 1L) "has" else "have",
                                         "another size")), collapse = ", and "),
      ", so its outputs may not be the run's.")))
  list(complete = TRUE, reason = NA_character_)
}

## What is at `dir`, against a run of key `key` from a model of md5 `model_md5`: `state` is "none"
## (nothing, or nothing but ignored files and files that are not a run's), "same" (a whole run of
## this key and model), "other" (a whole run of another) or "incomplete" (run files without a whole
## run); `foreign` names the entries that are not a run's, which nothing here deletes.
.run_state <- function(dir, key, model_md5) {
  if (!dir.exists(dir)) return(list(state = "none", foreign = character(), done = NULL))
  ent <- setdiff(list.files(dir, all.files = TRUE, no.. = TRUE), .ignored_files)
  done <- .read_done(dir)
  listed <- if (is.null(done)) character() else sub("^file:", "", grep("^file:", names(done), value = TRUE))
  ours <- intersect(ent, c(.run_files, listed))
  foreign <- setdiff(ent, ours)
  state <- if (!length(ours)) "none"
           else if (!.run_complete(dir)$complete) "incomplete"
           else if (identical(done[["run_key"]], key) && identical(done[["model_md5"]], model_md5)) "same"
           else "other"
  list(state = state, foreign = foreign, done = done)
}

## An error about what is at `outdir`, with a class callers can catch.
.outdir_in_use <- function(...)
  stop(structure(class = c("malariafwd_outdir_in_use", "error", "condition"),
                 list(message = paste0(...), call = NULL)))

## Removes `dir` if it holds nothing but run files and ignored files; never anything else.
.remove_run_dir <- function(dir) {
  if (!dir.exists(dir)) return(invisible(TRUE))
  ent <- setdiff(list.files(dir, all.files = TRUE, no.. = TRUE), .ignored_files)
  done <- .read_done(dir)
  listed <- if (is.null(done)) character() else sub("^file:", "", grep("^file:", names(done), value = TRUE))
  if (length(setdiff(ent, c(.run_files, listed)))) return(invisible(FALSE))
  invisible(unlink(dir, recursive = TRUE) == 0L)
}
