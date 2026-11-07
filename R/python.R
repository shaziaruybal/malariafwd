## python.R -- the one check that the tree-sequence functions can reach Python with tskit, pandas
## and numpy (through reticulate), and check_setup(), which reports everything the package uses.

## The modules those functions import, and the lowest version of each that works (NA: any).
.py_floors <- c(tskit = "1.0", pandas = NA, numpy = NA)
.py_requirements <- c("tskit>=1.0", "pandas", "numpy")

## A check that found everything holds for the session, since a Python that has started cannot
## change; one that found a problem is not kept, since py_install() or a network can mend it.
.py_cache <- new.env(parent = emptyenv())

## Thin wrappers around reticulate, so that each branch can be tested without Python.
.has_reticulate <- function() requireNamespace("reticulate", quietly = TRUE)
.reticulate_version <- function() as.character(utils::packageVersion("reticulate"))
.py_started <- function() reticulate::py_available(initialize = FALSE)
.py_declare <- function() reticulate::py_require(.py_requirements)
.py_start <- function() reticulate::py_config()
.py_module_version <- function(m)
  tryCatch(as.character(reticulate::import(m)$`__version__`), error = function(e) NA_character_)
.py_cached <- function() .py_cache$status

.and_list <- function(x, last = "and")
  if (length(x) < 2L) x else paste(paste(x[-length(x)], collapse = ", "), last, x[length(x)])

## A version's leading numbers, so that "1.0.4.dev3" reads as 1.0.4.
.leading_version <- function(v) numeric_version(sub("^([0-9]+(\\.[0-9]+)*).*$", "\\1", v), strict = FALSE)

## What the Python side looks like: which reticulate and Python, how that Python was chosen
## (`chosen_by`, reticulate's own environment being `managed`), each module's version (NA when it
## does not import), the warnings reticulate gave while starting, and `problem`, a sentence when
## Python cannot be used at all. Never errors or warns.
.python_status <- function() {
  hit <- .py_cached()
  if (!is.null(hit)) return(hit)
  st <- list(reticulate = NA_character_, python = NA_character_, python_version = NA_character_,
             chosen_by = NA_character_, managed = FALSE,
             versions = stats::setNames(rep(NA_character_, length(.py_floors)), names(.py_floors)),
             warnings = character(), problem = NULL)
  if (!.has_reticulate()) {
    st$problem <- "the reticulate package, which is not installed: install.packages(\"reticulate\")."
    return(st)
  }
  st$reticulate <- .reticulate_version()
  if (.leading_version(st$reticulate) < "1.41.0") {
    st$problem <- paste0("reticulate >= 1.41.0, and reticulate ", st$reticulate,
                         " is installed: install.packages(\"reticulate\").")
    return(st)
  }
  w <- character()
  quiet <- function(expr) withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  ## Declaring the modules only records them: a Python the user chose is used as it is, and
  ## reticulate builds its own environment with them only when none was chosen.
  if (!isTRUE(tryCatch(.py_started(), error = function(e) FALSE)))
    try(quiet(.py_declare()), silent = TRUE)
  cfg <- tryCatch(quiet(.py_start()), error = identity)
  if (inherits(cfg, "error")) {
    chosen <- Sys.getenv("RETICULATE_PYTHON")
    st$problem <- paste0("reticulate, which could not start Python: ",
                         strsplit(conditionMessage(cfg), "\n", fixed = TRUE)[[1]][1], ". ",
                         if (nzchar(chosen))
                           paste0("Check RETICULATE_PYTHON (", chosen, "), then restart R.")
                         else paste("With no Python chosen, reticulate sets up its own, which needs",
                                    "internet the first time; or set RETICULATE_PYTHON to a Python that",
                                    "has the modules and restart R."))
    st$warnings <- w
    return(st)
  }
  st$python <- as.character(cfg$python)
  st$python_version <- as.character(cfg$version)
  forced <- cfg$forced
  st$chosen_by <- if (is.character(forced) && length(forced) == 1L && nzchar(forced)) forced else NA_character_
  st$managed <- isTRUE(cfg$ephemeral)
  read <- function() vapply(names(.py_floors), .py_module_version, character(1))
  st$versions <- read()
  ## reticulate adds what is declared to its own environment after it has started
  if (anyNA(st$versions) && st$managed) {
    try(quiet(.py_declare()), silent = TRUE)
    st$versions <- read()
  }
  st$warnings <- w
  if (is.null(.py_problem(st))) .py_cache$status <- st
  st
}

## The sentence that follows "<function> needs" when `modules` cannot be used, or NULL.
.py_problem <- function(st, modules = names(.py_floors)) {
  floors <- .py_floors[modules]
  need <- paste0("Python with ", .and_list(ifelse(is.na(floors), modules,
                                                  paste0(modules, " (>= ", floors, ")"))))
  if (!is.null(st$problem)) return(paste0(need, ", through ", st$problem, " check_setup() shows what is found."))
  v <- st$versions[modules]
  missing <- modules[is.na(v)]
  low <- modules[!is.na(v) & !is.na(floors) & .leading_version(v) < .leading_version(floors)]
  if (!length(missing) && !length(low)) return(NULL)
  how <- if (st$managed) "reticulate's own environment"
         else paste0("Python ", st$python_version, if (!is.na(st$chosen_by)) paste0(", chosen by ", st$chosen_by))
  using <- paste0("The Python reticulate is using, ", st$python, " (", how, "), ")
  what <- c(if (length(missing)) paste("cannot import", .and_list(missing, "or")),
            if (length(low)) paste0("has ", paste0(low, " ", v[low], collapse = " and "), "; malariafwd needs ",
                                    paste0(low, " >= ", floors[low], collapse = " and ")))
  fix <- if (st$managed)
    paste0("Adding them to reticulate's own environment failed",
           if (length(st$warnings)) paste0(" (", st$warnings[1], ")"),
           "; check the internet connection, or set RETICULATE_PYTHON to a Python that has them and ",
           "restart R.")
  else paste0("Install them into that Python (for example `", st$python, " -m pip install ",
              paste0("\"", .py_requirements, "\"", collapse = " "), "`, or reticulate::py_install() for a ",
              "virtualenv or conda environment), or set RETICULATE_PYTHON to a Python that has them and ",
              "restart R.",
              if (length(st$warnings)) paste0(" reticulate said: ", st$warnings[1]))
  paste0(need, ". ", using, paste(what, collapse = ", and "), ". ", fix, " check_setup() shows what is found.")
}

## Stops, naming `what`, unless Python and `modules` can be used; called by every function that
## runs Python, after its own argument checks and before it reads a file.
.require_python <- function(what, modules = names(.py_floors)) {
  st <- .python_status()
  p <- .py_problem(st, modules)
  if (!is.null(p))
    stop(structure(class = c("malariafwd_no_python", "error", "condition"),
                   list(message = paste(what, "needs", p), call = NULL)))
  invisible(st)
}

#' Can the tree-sequence functions reach Python with tskit?
#'
#' `TRUE` when Python with tskit (>= 1.0), pandas and numpy can be used through reticulate. It
#' never errors or warns, so it can guard examples, vignette chunks and tests. It may start
#' Python, and with no Python chosen, reticulate then sets up its own environment.
#'
#' @section Python:
#' The IBD and ARG functions read a run's tree sequence with Python's tskit, through the
#' reticulate package. Either of two routes works. *Automatic*: with no Python chosen, reticulate
#' sets up its own environment with tskit (>= 1.0), pandas and numpy the first time one of these
#' functions runs; this needs internet once. *Your own Python*: set `RETICULATE_PYTHON` to a
#' Python that has those modules (in `~/.Renviron` to keep it) before the first such call of a
#' session. A function that cannot reach them stops with an error of class
#' `malariafwd_no_python` that says what to do, and [check_setup()] reports what is found.
#' @return `TRUE` or `FALSE`.
#' @examples
#' tskit_available()
#' @export
tskit_available <- function()
  isTRUE(tryCatch(is.null(.py_problem(.python_status())), error = function(e) FALSE))

#' Check what malariafwd needs: R, SLiM, the model, and Python with tskit
#'
#' Reports each component the package uses, with its version and where it was found, and what to
#' do about any that is missing. It never errors. SLiM and the model are needed to run; reticulate,
#' Python, tskit, pandas and numpy for the IBD and ARG functions (see [tskit_available()]); ggplot2
#' and patchwork for the plots.
#' @param python if `FALSE`, Python is not started and its rows are not checked.
#' @return a data.frame of class `malariafwd_setup`, one row per component, with columns
#'   `component`, `ok` (`NA` when not checked), `version`, `detail` and `advice` (`NA` when there
#'   is nothing to do).
#' @examples
#' check_setup(python = FALSE)
#' @export
check_setup <- function(python = TRUE) {
  row <- function(component, ok, version = NA_character_, detail = NA_character_, advice = NA_character_)
    data.frame(component = component, ok = ok, version = as.character(version), detail = detail,
               advice = if (isTRUE(ok) || is.na(ok)) NA_character_ else advice, stringsAsFactors = FALSE)
  pkg_version <- function(p) if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else NA_character_
  rows <- list(row("malariafwd", TRUE, pkg_version("malariafwd"), dirname(system.file(package = "malariafwd"))),
               row("R", TRUE, paste(R.version$major, R.version$minor, sep = "."), R.home()))
  ## SLiM
  env <- Sys.getenv("MALARIAFWD_SLIM")
  where <- paste0(slim_path(), if (nzchar(env)) " (from MALARIAFWD_SLIM)")
  sv <- tryCatch(if (slim_available()) slim_version() else NA_character_, error = function(e) NA_character_)
  rows[[length(rows) + 1]] <- if (!slim_available())
    row("SLiM", FALSE, NA, where, paste0("install SLiM >= ", .slim_floor, " and put it on PATH, or set ",
        "MALARIAFWD_SLIM=/path/to/slim in ~/.Renviron (RStudio on macOS may not see the shell's PATH)",
        if (nzchar(env)) paste0("; MALARIAFWD_SLIM is ", env, ", which was not found")))
  else row("SLiM", .slim_version_ok(sv) && !is.na(sv), sv, where,
           paste0("needs SLiM >= ", .slim_floor, ", found ", if (is.na(sv)) "an unreadable version" else sv))
  m <- slim_model()
  rows[[length(rows) + 1]] <- row("model", nzchar(m) && file.exists(m), NA, m,
                                  "the bundled model was not found: reinstall malariafwd")
  ## reticulate, then Python and its modules
  rv <- pkg_version("reticulate")
  rows[[length(rows) + 1]] <- row("reticulate", !is.na(rv) && .leading_version(rv) >= "1.41.0", rv,
                                  "for the IBD and ARG functions", "install.packages(\"reticulate\")")
  if (!python || is.na(rv)) {
    for (cmp in c("Python", names(.py_floors)))
      rows[[length(rows) + 1]] <- row(cmp, NA, NA, if (!python) "not checked" else "needs reticulate")
  } else {
    st <- .python_status()
    started <- is.null(st$problem)
    how <- if (st$managed) paste0("reticulate's own environment, in ", tools::R_user_dir("reticulate", "cache"))
           else if (!is.na(st$chosen_by)) paste("chosen by", st$chosen_by) else NA_character_
    rows[[length(rows) + 1]] <- row("Python", started, st$python_version,
                                    if (started) paste0(st$python, if (!is.na(how)) paste0(" (", how, ")")) else NA,
                                    .py_problem(st))
    for (mod in names(.py_floors))
      rows[[length(rows) + 1]] <- row(mod, started && is.null(.py_problem(st, mod)), st$versions[[mod]], NA,
                                      .py_problem(st, mod))
  }
  for (p in c("ggplot2", "patchwork")) {
    v <- pkg_version(p)
    rows[[length(rows) + 1]] <- row(p, !is.na(v), v, "optional, for the plots", paste0("install.packages(\"", p, "\")"))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  structure(out, class = c("malariafwd_setup", "data.frame"))
}

#' @export
print.malariafwd_setup <- function(x, ...) {
  d <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
  for (i in seq_len(nrow(d))) {
    mark <- if (is.na(d$ok[i])) "[  ]" else if (d$ok[i]) "[ok]" else "[--]"
    cat(sprintf("  %s %-11s %-10s %s\n", mark, d$component[i], if (is.na(d$version[i])) "" else d$version[i],
                if (is.na(d$detail[i])) "" else d$detail[i]))
    if (!is.na(d$advice[i])) cat("       ", d$advice[i], "\n", sep = "")
  }
  optional <- d$component %in% c("ggplot2", "patchwork")
  bad <- sum(d$ok %in% FALSE & !optional); bad_opt <- sum(d$ok %in% FALSE & optional)
  cat(if (bad) paste0(bad, " problem", if (bad > 1) "s", ": see the advice above.")
      else if (anyNA(d$ok)) "Everything checked is in place."
      else "Everything the package needs is in place.",
      if (bad_opt) paste0(" ", bad_opt, " optional package", if (bad_opt > 1) "s are" else " is",
                          " missing, for the plots."), "\n", sep = "")
  invisible(x)
}
