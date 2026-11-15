## The bundled example run: the settings and seed data-raw/example_run.R runs it with. The tests
## compare the shipped run with them and with the model, so a change to a default or to the model
## says that the run must be made again.
.example_run_spec <- list(
  seed = 2L,
  params = list(N_HOSTS = 60, N_MOSQ = 150, BITING_RATE = 0.1, SIM_DAYS = 400,
                N_FOUNDER_HOSTS = 16, N_FOUNDER_LINEAGES = 10,
                LAMBDA_OOCYSTS = 1, LAMBDA_SPOROZOITES = 1.5,
                IMPORT_RATE_PER_TICK = 0.05, IMPORT_START = 1,
                BLOODCURE_COVERAGE = 0.5, RADCURE_COVERAGE = 0.3, P_BLOODCURE_SURVIVAL = 0.1,
                P_RECRUD_GIVEN_FAIL = 0.5, P_DETECT_RELAPSE = 0.5))

#' A small run that ships with the package
#'
#' A run of the package's own model, small enough to bundle: 60 hosts and 150 mosquitoes for 400
#' days, with imports, case management, recrudescence and co-transmission switched on, so that
#' every kind of case, recurrence and transmission link the package names occurs in it. Every
#' function that reads a run takes it, so the examples, the vignettes and a first look at the
#' package need no SLiM;
#' the relatedness functions still need Python with tskit ([tskit_available()]).
#'
#' The first call in a session unpacks the run into a temporary directory, and later calls return
#' the same run, unpacking it again if a file in its directory is missing or has changed. Its
#' settings are `load_params(example_run())` and its seed `example_run()$seed`;
#' `data-raw/example_run.R`, in the package's source repository, makes it.
#' @return a `vivax_run` (see [load_run()]).
#' @seealso [load_run()], [run_sim()]
#' @examples
#' run <- example_run()
#' run
#' table(classify_cases(run)$origin)
#' @export
example_run <- function() {
  tarball <- system.file("extdata", "example_run.tar.xz", package = "malariafwd")
  if (!nzchar(tarball)) stop("the example run is not installed with the package.", call. = FALSE)
  ## one directory per archive, so a new archive is unpacked afresh; unpacked beside it and moved in
  ## whole, so that the directory, like one run_sim() leaves, holds the whole run or nothing
  dir <- file.path(tempdir(), paste0("malariafwd-example-run-", substr(unname(tools::md5sum(tarball)), 1, 8)))
  if (!isTRUE(.run_complete(dir)$complete)) {
    tmp <- tempfile("malariafwd-example-run-")
    utils::untar(tarball, exdir = tmp, tar = "internal")
    if (dir.exists(dir) && !isTRUE(.run_complete(dir)$complete)) unlink(dir, recursive = TRUE)
    if (!file.rename(tmp, dir)) unlink(tmp, recursive = TRUE)   # another process moved its copy in first
    if (!isTRUE(.run_complete(dir)$complete)) stop("could not unpack the example run into ", dir, call. = FALSE)
  }
  load_run(dir)
}
