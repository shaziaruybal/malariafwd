## slim.R -- locating and checking the SLiM binary + the bundled model template.

#' Path to the SLiM binary
#'
#' Uses the `MALARIAFWD_SLIM` environment variable if set, otherwise `"slim"` on PATH.
#' @return the resolved path (or the bare name if not found on PATH).
#' @export
slim_path <- function() {
  p <- Sys.getenv("MALARIAFWD_SLIM", unset = "slim")
  found <- Sys.which(p)
  if (nzchar(found)) unname(found) else p
}

#' Is a SLiM binary available?
#' @return `TRUE` if SLiM is found on PATH (or via `MALARIAFWD_SLIM`).
#' @export
slim_available <- function() nzchar(Sys.which(slim_path()))

#' SLiM version string
#' @return the version (e.g. `"5.2"`) or `NA` if SLiM is not found.
#' @export
slim_version <- function() {
  if (!slim_available()) return(NA_character_)
  out <- tryCatch(processx::run(slim_path(), "-v")$stdout, error = function(e) "")
  l1 <- strsplit(out, "\n")[[1]][1]
  if (grepl("SLiM version", l1)) sub(".*SLiM version ([0-9.]+).*", "\\1", l1) else NA_character_
}

#' Path to the bundled SLiM model template
#'
#' The single source of truth for the model; `run_sim()` snapshots it into each run's
#' output directory for provenance.
#' @return file path to `inst/slim/malariafwd_vivax.slim`.
#' @export
slim_model <- function() system.file("slim", "malariafwd_vivax.slim", package = "malariafwd")
