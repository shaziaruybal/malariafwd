## grid.R -- parameter set x seed sweeps. run_grid() checks every set and every cell's directory
## before any cell starts, runs the cells in parallel, resumes from what is on disk, and records
## the grid in grid.tsv, which load_grid() reads (D59).

## A set's name is its cells' directory names, so it must be safe as one on every file system;
## blank names are filled by position.
.grid_names <- function(sets) {
  nm <- names(sets)
  if (is.null(nm)) nm <- rep("", length(sets))
  blank <- is.na(nm) | !nzchar(nm)
  nm[blank] <- sprintf("set%02d", which(blank))
  bad <- nm[grepl("[/\\\\:*?\"<>|\\x00-\\x1f\\x7f]", nm, perl = TRUE) | nm %in% c(".", "..") |
              grepl("^[ .]|[ .]$", nm)]
  if (length(bad))
    stop("run_grid(): a set's name names its cells' directories, so it must not hold / \\ : * ? \" < > | ",
         "or a control character, be . or .., or start or end with a space or a dot: ",
         paste0("\"", bad, "\"", collapse = ", "), ".", call. = FALSE)
  dup <- unique(nm[duplicated(tolower(nm)) | duplicated(tolower(nm), fromLast = TRUE)])
  if (length(dup))
    stop("run_grid(): set names must differ, also ignoring case (the file systems of macOS and Windows ",
         "ignore it): ", paste(dup, collapse = ", "), ".", call. = FALSE)
  nm
}

.check_seeds <- function(seeds) {
  if (!is.numeric(seeds) || !length(seeds) || anyNA(seeds) || any(!is.finite(seeds)) ||
      any(seeds != round(seeds)) || any(abs(seeds) > .Machine$integer.max))
    stop("`seeds` must be one or more whole numbers (|seed| <= ", .Machine$integer.max, "), none missing.",
         call. = FALSE)
  if (anyDuplicated(seeds))
    stop("`seeds` must not repeat a seed: ", paste(unique(seeds[duplicated(seeds)]), collapse = ", "), ".",
         call. = FALSE)
  as.integer(seeds)
}

## Each parameter's value as SLiM is given it, as text: numbers as as_slim_args() writes them,
## logicals T or F, strings without their Eidos quotes.
.slim_values <- function(p) {
  tok <- as_slim_args(p)[c(FALSE, TRUE)]
  nm <- sub("=.*$", "", tok); v <- sub("^[^=]*=", "", tok)
  q <- grepl("^'.*'$", v)
  v[q] <- gsub("\\\\(.)", "\\1", substr(v[q], 2L, nchar(v[q]) - 1L))
  stats::setNames(v, nm)
}

## The parameters whose values, as SLiM is given them, differ between the sets.
.varying_params <- function(sets) {
  v <- vapply(sets, .slim_values, character(length(.malariafwd_defaults)))
  if (is.null(dim(v))) return(character())
  rownames(v)[apply(v, 1L, function(x) length(unique(x)) > 1L)]
}

## A column of parameter `name` from its text, typed as the parameter's default.
.typed_column <- function(name, values, where) {
  if (!name %in% names(.malariafwd_defaults)) return(values)
  v <- .type_param_text(name, values, where)
  if (is.integer(.malariafwd_defaults[[name]]) && all(abs(v) <= .Machine$integer.max)) as.integer(v) else v
}

## The grid's record in grid.tsv: this call's rows replace those of the same set and seed; the
## others are kept. Written beside the file and renamed, so a reader never sees half of it.
.grid_columns <- c(set = "character", seed = "integer", dir = "character", ok = "logical",
                   skipped = "logical", error = "character", extinct_tick = "integer",
                   never_started = "logical", elapsed_s = "double", finished = "character",
                   slim_version = "character", malariafwd_version = "character",
                   run_key = "character", model_md5 = "character")

.write_grid_manifest <- function(base, rows) {
  f <- file.path(base, "grid.tsv")
  old <- if (file.exists(f)) tryCatch(.read_grid_manifest(base), error = function(e) NULL)
  if (!is.null(old)) {
    old <- old[!paste(old$set, old$seed) %in% paste(rows$set, rows$seed)]
    rows <- data.table::rbindlist(list(old, rows), use.names = TRUE, fill = TRUE)
  }
  tmp <- tempfile("grid.tsv.tmp-", tmpdir = base)
  data.table::fwrite(rows, tmp, sep = "\t", na = "NA", quote = TRUE)
  if (!file.rename(tmp, f)) { unlink(tmp); warning("grid.tsv could not be written in ", base, call. = FALSE) }
  invisible(f)
}

## grid.tsv as written, every column as text (the parameters are typed by the caller).
.read_grid_manifest <- function(base) {
  d <- data.table::fread(file.path(base, "grid.tsv"), sep = "\t", colClasses = "character",
                         na.strings = NULL, showProgress = FALSE)
  cell <- intersect(names(.grid_columns), names(d))
  for (nm in cell) {
    x <- d[[nm]]; x[x == "NA"] <- NA
    data.table::set(d, j = nm, value = switch(.grid_columns[[nm]], integer = as.integer(x),
                                              double = as.numeric(x), logical = as.logical(x), x))
  }
  d[]
}

#' Run a grid of parameter sets x seeds in parallel
#'
#' Runs every (parameter set, seed) cell, each in its own directory `<set>_seed<seed>` under
#' `outdir_base`, in parallel with the current `future` plan. Set `workers` for a quick local
#' multisession; for a cluster, set `future::plan(...)` yourself before calling and leave
#' `workers = NULL` (see DECISIONS.md D15).
#'
#' Every set is checked, and every cell's directory, before any cell starts: a set
#' [vivax_params()] refuses, or a cell whose directory [run_sim()] would refuse, stops the grid
#' with all of them listed. A cell that already holds its run (the same parameters, seed and
#' model) is not run again, so calling `run_grid()` again with the same sets and seeds resumes a
#' grid that stopped part-way, and adding seeds or sets runs only the new cells. The cells of an
#' edited set hold another run and stop the grid unless `overwrite = TRUE`.
#'
#' `grid.tsv` in `outdir_base`, written when the grid finishes (also when cells failed), records
#' each cell: its outcome, versions, the key of its run and every parameter. [load_grid()] reads
#' it. Cells run by other means, such as one [run_sim()] per job of a job array, are recorded by
#' calling `run_grid()` with the same sets and seeds afterwards, which runs nothing. Two grids
#' writing into one `outdir_base` at the same time race on `grid.tsv`; a later call rebuilds it.
#'
#' @param param_sets a [vivax_params()] object, or a list of them (or of lists of overrides). Its
#'   names name the cells' directories: they must differ, also ignoring case, and be usable as a
#'   directory name; blank names become `set01`, `set02`, ... by position.
#' @param seeds whole-number seeds to run for each set, none repeated.
#' @param outdir_base the grid's directory, created if needed.
#' @param model path to the SLiM model.
#' @param workers if set, run under `future::plan(multisession, workers)` for this call (the
#'   previous plan is restored afterward). Set it to the cores you can spare. If `NULL`, use the
#'   current plan.
#' @param verbose passed to [run_sim()].
#' @param skip_existing if `TRUE` (the default), a cell that already holds its run is not run again.
#' @param overwrite if `TRUE`, a cell that holds another run, or one that did not finish, is run
#'   again and its run replaced once the new one has finished. Files that are not part of a run are
#'   never replaced.
#' @param ... passed to [run_sim()].
#' @return a `data.table` with one row per cell, by set then seed: `set`, `seed`, `outdir`, `ok`
#'   (the cell holds its run, run now or before), `skipped` (it was not run again), `error` (the
#'   failure message, or `NA`), `extinct_tick` (the tick on which the run went extinct and ended,
#'   `NA` if it did not) and `never_started` (`TRUE` for a founder-free run whose imports never
#'   established a parasite; see [load_run()]), both `NA` for a failed cell, `elapsed_s`,
#'   `finished` (UTC), `slim_version` and `malariafwd_version`, from the run's `DONE`, and then one
#'   column per parameter whose value differs between the sets, of the parameter's type. A cell
#'   that never started is `ok`: not starting is an outcome of the grid, which `never_started`
#'   records. Its run holds no parasite, though, so give only the cells with
#'   `ok & !never_started` to the readers that need one, such as [arg_summary()]. The cells that
#'   went extinct, those that never started, and those not run again are counted in one message
#'   each. A failed cell does not abort the grid -- its directory is left as it was and its files
#'   and log are in `<cell>.failed` -- but it warns, and the reason travels with the result. A
#'   warning given in the cells (such as a parameter set's own warnings, which [run_sim()]
#'   repeats for each run) is given once, after the FAILED summary, with the number of cells it
#'   came from. If the grid stops part-way with an error, that error's message lists the warnings
#'   of the cells whose results had reached this session; with workers, a cell's warnings travel
#'   with its chunk's result, so those of a chunk lost with its worker are lost. Under
#'   `options(warn = 2)` a warning in a cell fails that cell, in any plan and also inside
#'   `suppressWarnings()`, and the FAILED summary is then an error, so nothing is returned; the
#'   cells that ran are on disk. The caller's random-number stream is left as it was.
#' @examples
#' \dontrun{
#' sets <- list(low = vivax_params(BITING_RATE = 0.05), high = vivax_params(BITING_RATE = 0.12))
#' g <- run_grid(sets, seeds = 1:3, outdir_base = file.path(tempdir(), "sweep"), workers = 2)
#' load_grid(file.path(tempdir(), "sweep"))
#' }
#' @export
run_grid <- function(param_sets, seeds = 1L, outdir_base, model = slim_model(),
                     workers = NULL, verbose = NULL, skip_existing = TRUE, overwrite = FALSE, ...) {
  if (inherits(param_sets, "vivax_params")) param_sets <- list(set01 = param_sets)
  if (!is.list(param_sets) || !length(param_sets) ||
      !all(vapply(param_sets, function(x) inherits(x, "vivax_params") || (is.list(x) && !is.object(x)),
                  logical(1))))
    stop("param_sets must be a vivax_params() set or a list of them.", call. = FALSE)
  names(param_sets) <- .grid_names(param_sets)
  seeds <- .check_seeds(seeds)
  if (missing(outdir_base) || !is.character(outdir_base) || length(outdir_base) != 1L ||
      is.na(outdir_base) || !nzchar(outdir_base))
    stop("`outdir_base` must be one directory path.", call. = FALSE)
  for (f in c("overwrite", "skip_existing"))
    if (!isTRUE(get(f)) && !isFALSE(get(f))) stop("`", f, "` must be TRUE or FALSE.", call. = FALSE)
  if (!is.character(model) || length(model) != 1L || is.na(model) || !nzchar(model))
    stop("model template not found (is the package installed?).", call. = FALSE)
  if (!file.exists(model)) stop("model file not found: ", model, call. = FALSE)
  model <- normalizePath(model, winslash = "/")
  model_md5 <- unname(tools::md5sum(model))

  ## Every set is checked before any cell starts; the warnings they give come from the cells.
  refused <- character()
  sets <- lapply(names(param_sets), function(nm) tryCatch(suppressWarnings({
    p <- param_sets[[nm]]
    if (inherits(p, "vivax_params")) .refuse_unknown(names(p), built = TRUE)
    do.call(vivax_params, as.list(unclass(p)), quote = TRUE)
  }), error = function(e) { refused[[nm]] <<- conditionMessage(e); NULL }))
  names(sets) <- names(param_sets)
  if (length(refused))
    stop("run_grid() did not start: ", length(refused), " of ", length(sets), " parameter sets ",
         if (length(refused) == 1L) "is" else "are", " refused:\n",
         paste0("  ", names(refused), ": ", refused, collapse = "\n"), call. = FALSE)

  jobs <- data.frame(set = rep(names(sets), each = length(seeds)), seed = rep(seeds, times = length(sets)),
                     stringsAsFactors = FALSE)
  jobs$dir <- sprintf("%s_seed%d", jobs$set, jobs$seed)
  jobs$key <- vapply(seq_len(nrow(jobs)), function(i) .run_key(sets[[jobs$set[i]]], jobs$seed[i]), "")
  ## and every cell's directory, against what run_sim() would do with it
  states <- lapply(seq_len(nrow(jobs)), function(i)
    .run_state(file.path(outdir_base, jobs$dir[i]), jobs$key[i], model_md5))
  st <- vapply(states, `[[`, "", "state")
  foreign <- vapply(states, function(x) length(x$foreign) > 0L, logical(1))
  jobs$skip <- st == "same" & skip_existing
  conflict <- !jobs$skip & (foreign | (st != "none" & !overwrite))
  if (any(conflict))
    stop("run_grid() did not start: ", sum(conflict), " of ", nrow(jobs), " cells already hold a different ",
         "run, one that did not finish, or other files: ", paste(utils::head(jobs$dir[conflict], 10), collapse = ", "),
         if (sum(conflict) > 10) ", ...", ". overwrite = TRUE replaces their runs once the new ones have ",
         "finished; files that are not part of a run are never replaced. Or use another outdir_base.",
         call. = FALSE)

  dir.create(outdir_base, recursive = TRUE, showWarnings = FALSE)
  base <- normalizePath(outdir_base, winslash = "/", mustWork = TRUE)
  jobs$outdir <- file.path(base, jobs$dir)
  todo <- which(!jobs$skip)

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
    ## run_sim() supervises SLiM so that it stops if R dies; processx then keeps two FIFOs open for
    ## the rest of the session, which future reports for the cell that opened them
    if (inherits(w, "ConnectionMisuseFutureWarning") && grepl("supervisor_std", conditionMessage(w), fixed = TRUE))
      invokeRestart("muffleWarning")
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

  out <- list()
  if (length(todo)) {
    .require_slim()
    if (!is.null(workers)) {
      oldplan <- future::plan(future::multisession, workers = workers)
      on.exit(future::plan(oldplan), add = TRUE)
    }
    out <- tryCatch(withCallingHandlers(furrr::future_pmap(
      list(jobs$set[todo], jobs$seed[todo], jobs$outdir[todo]),
      function(set, seed, outdir) {
        ## Keep the error TEXT: a cell that failed because SLiM is not on PATH and one that failed
        ## on a bad parameter must be told apart. The grid still does not abort - one bad cell
        ## should not lose a grid - and the cell's directory is left as it was, its files and log
        ## in <cell>.failed.
        err <- NA_character_; extinct <- NA_integer_; never <- NA
        run <- function() run_sim(sets[[set]], seed = seed, outdir = outdir, model = model,
                                  verbose = verbose, overwrite = overwrite, skip_existing = FALSE, ...)
        ok <- tryCatch({
          ## the extinction message is summed up once for the grid, below
          r <- suppressMessages(if (warn_level < 2) run() else withCallingHandlers(run(), warning =
            function(w) stop("(converted from warning) ", conditionMessage(w), call. = FALSE)))
          extinct <- r$extinct_tick; never <- isTRUE(r$never_started)
          TRUE
        }, error = function(e) { err <<- conditionMessage(e); FALSE })
        data.frame(set = set, seed = seed, outdir = outdir, ok = ok, skipped = FALSE, error = err,
                   extinct_tick = extinct, never_started = never, stringsAsFactors = FALSE)
      },
      .options = furrr::furrr_options(seed = FALSE)
    ), warning = keep), error = function(e) {
      if (length(msgs)) {
        e$message <- paste0(conditionMessage(e), "\nThe cells that had reported gave these warnings:\n",
                            paste0("  ", counted(), collapse = "\n"))
        msgs <<- character(); conds <<- list(); counts <<- integer()
      }
      stop(e)
    })
  }
  ## the cells that already held their run, from their directories
  for (i in which(jobs$skip)) {
    r <- load_run(jobs$outdir[i])
    out[[length(out) + 1L]] <- data.frame(set = jobs$set[i], seed = jobs$seed[i], outdir = jobs$outdir[i],
                                         ok = TRUE, skipped = TRUE, error = NA_character_,
                                         extinct_tick = r$extinct_tick, never_started = isTRUE(r$never_started),
                                         stringsAsFactors = FALSE)
  }
  res <- data.table::rbindlist(out)
  res <- res[match(jobs$outdir, res$outdir)]
  done <- lapply(res$outdir, .read_done)
  from_done <- function(k) vapply(done, function(d) if (is.null(d) || is.null(d[[k]])) NA_character_ else d[[k]], "")
  res[, `:=`(elapsed_s = suppressWarnings(as.numeric(from_done("elapsed_s"))), finished = from_done("finished"),
             slim_version = from_done("slim_version"), malariafwd_version = from_done("malariafwd_version"))]
  vary <- .varying_params(sets)
  for (nm in vary) data.table::set(res, j = nm, value = unlist(lapply(sets[res$set], `[[`, nm), use.names = FALSE))

  ## the record, also when cells failed or none ran
  vals <- lapply(sets, .slim_values)
  rec <- data.table::data.table(set = res$set, seed = res$seed, dir = jobs$dir, ok = res$ok, skipped = res$skipped,
                                error = res$error, extinct_tick = res$extinct_tick, never_started = res$never_started,
                                elapsed_s = res$elapsed_s, finished = res$finished, slim_version = res$slim_version,
                                malariafwd_version = res$malariafwd_version, run_key = jobs$key,
                                model_md5 = model_md5)
  pv <- do.call(rbind, lapply(res$set, function(s) vals[[s]][names(.malariafwd_defaults)]))
  rec <- cbind(rec, data.table::as.data.table(pv))
  .write_grid_manifest(base, rec)

  ## Say it out loud. A grid is usually launched and left, and a silent half-empty result is
  ## how a sweep gets analysed with a third of its cells missing. Under warn = 2 the summary is an
  ## error, and nothing is returned.
  held <- list(msgs, conds, counts); msgs <- character(); conds <- list(); counts <- integer()
  if (any(!res$ok)) {
    n <- sum(!res$ok)
    warning(n, " of ", nrow(res), " grid cell(s) FAILED: their directories were left as they were, and ",
            "each cell's files and log are in <cell>.failed. First error: ", res$error[!res$ok][1],
            if (warn_level < 2) "\nSee the `error` column for the rest.", call. = FALSE)
  }
  msgs <- held[[1]]; conds <- held[[2]]; counts <- held[[3]]
  give()
  if (any(res$skipped))
    message(sum(res$skipped), " of ", nrow(res), " grid cell(s) already held their run and were not run ",
            "again; see the `skipped` column.")
  if (any(!is.na(res$extinct_tick)))
    message(sum(!is.na(res$extinct_tick)), " of ", nrow(res), " grid cell(s) went extinct; see the ",
            "`extinct_tick` column.")
  if (any(res$never_started %in% TRUE))
    message(sum(res$never_started %in% TRUE), " of ", nrow(res), " grid cell(s) never started: no ",
            "import established a parasite; see the `never_started` column.")
  res[]
}

#' Load the record of a grid
#'
#' Reads the `grid.tsv` that [run_grid()] writes in its `outdir_base`: one row per cell, with
#' [run_grid()]'s columns. `outdir` is resolved against `outdir_base`, so a grid that was moved or
#' copied still loads. It warns, naming up to five, when a cell recorded as `ok` no longer holds
#' its whole run, and says so when the cells were run by different versions of SLiM or malariafwd.
#' @param outdir_base the grid's directory.
#' @param params `"varying"` (the default) for the parameters whose values differ between the
#'   cells, `"all"` for every parameter, each of its type.
#' @return a `data.table`, as [run_grid()] returns.
#' @examples
#' \dontrun{
#' run_grid(vivax_params(), seeds = 1:2, outdir_base = file.path(tempdir(), "sweep"))
#' load_grid(file.path(tempdir(), "sweep"))
#' }
#' @export
load_grid <- function(outdir_base, params = c("varying", "all")) {
  params <- match.arg(params)
  if (!is.character(outdir_base) || length(outdir_base) != 1L || is.na(outdir_base) || !nzchar(outdir_base))
    stop("`outdir_base` must be one directory path.", call. = FALSE)
  base <- normalizePath(outdir_base, winslash = "/", mustWork = FALSE)
  if (!file.exists(file.path(base, "grid.tsv")))
    stop("no grid.tsv in ", base, ": run_grid() writes it when it finishes. Call run_grid() again with the ",
         "same sets and seeds; cells that already hold their run are not run again.", call. = FALSE)
  d <- .read_grid_manifest(base)
  where <- paste("grid.tsv in", base)
  cell <- c("set", "seed", "outdir", "ok", "skipped", "error", "extinct_tick", "never_started", "elapsed_s",
            "finished", "slim_version", "malariafwd_version")
  data.table::set(d, j = "outdir", value = file.path(base, d$dir))
  pcols <- setdiff(names(d), c(names(.grid_columns), "outdir"))
  keep <- if (params == "all") pcols else pcols[vapply(pcols, function(nm) length(unique(d[[nm]])) > 1L, logical(1))]
  for (nm in keep) data.table::set(d, j = nm, value = .typed_column(nm, d[[nm]], where))
  out <- d[, c(cell, keep), with = FALSE]
  gone <- out$outdir[out$ok %in% TRUE & !vapply(out$outdir, function(x) dir.exists(x) && .run_complete(x)$complete, logical(1))]
  if (length(gone))
    warning(length(gone), " cell(s) recorded as ok no longer hold their whole run: ",
            paste(basename(utils::head(gone, 5)), collapse = ", "), if (length(gone) > 5) ", ...", ".", call. = FALSE)
  v <- unique(stats::na.omit(out$slim_version)); m <- unique(stats::na.omit(out$malariafwd_version))
  if (length(v) > 1L || length(m) > 1L)
    message("The grid's cells were run by ", paste(c(if (length(v) > 1L) paste("SLiM", paste(v, collapse = ", ")),
                                                      if (length(m) > 1L) paste("malariafwd", paste(m, collapse = ", "))),
                                                    collapse = " and "), ".")
  out[]
}
