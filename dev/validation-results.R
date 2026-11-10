## Runs the validation scripts of inst/validation on this source tree and writes
## inst/validation/RESULTS.md: when, on which commit and versions, and what each script found.
## Run it from the package's root before a release, with SLiM and Python with tskit available:
##
##   Rscript dev/validation-results.R
##
## The package is installed into a temporary library for the run, never into yours.

rbin <- function(x) file.path(R.home("bin"), x)
lib <- file.path(tempdir(), "lib"); dir.create(lib, showWarnings = FALSE)
if (system2(rbin("R"), c("CMD", "INSTALL", "--no-test-load", "-l", shQuote(lib), "."), stdout = FALSE, stderr = FALSE) != 0L)
  stop("R CMD INSTALL failed", call. = FALSE)
env <- paste0("R_LIBS=", lib)
rscript <- function(code, log) system2(rbin("Rscript"), c("-e", shQuote(code)), stdout = log, stderr = log, env = env)

## the versions, as the installed package sees them (no paths)
vf <- tempfile()
rscript(paste0("s <- malariafwd::check_setup(); v <- stats::setNames(s$version, s$component); ",
               "writeLines(paste(names(v), v, sep = '\\t'), '", vf, "')"), tempfile())
ver <- utils::read.delim(vf, header = FALSE, col.names = c("component", "version"), colClasses = "character")
v <- function(x) { y <- ver$version[ver$component == x]; if (length(y) && !is.na(y[1]) && nzchar(y[1])) y[1] else "not found" }

scripts <- c("audit_events_and_ibd.R", "check_parameters.R", "audit_arg_recombination.R", "audit_arg_junctions.R")
res <- do.call(rbind, lapply(scripts, function(s) {
  log <- tempfile(fileext = ".log"); t0 <- Sys.time()
  status <- rscript(sprintf("source(system.file('validation', '%s', package = 'malariafwd'))", s), log)
  out <- readLines(log, warn = FALSE)
  line <- utils::tail(grep("[0-9]+ passed, [0-9]+ failed", out, value = TRUE), 1)
  data.frame(script = s,
             result = if (length(line)) sub(".*?([0-9]+ passed, [0-9]+ failed).*", "\\1", line)
                      else if (status == 0L) "ran without error (a numeric audit: run it to see its numbers)" else "stopped with an error",
             ok = status == 0L && (!length(line) || grepl(" 0 failed", line)),
             minutes = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1))
}))

git <- function(...) suppressWarnings(system2("git", c(...), stdout = TRUE, stderr = FALSE))
## R CMD INSTALL took everything under R/, inst/ and man/, untracked files too, so any change there
## or to DESCRIPTION and NAMESPACE means the run was not of the commit alone
changed <- git("status", "--porcelain", "--", "R", "inst", "man", "DESCRIPTION", "NAMESPACE")
changed <- changed[!grepl("inst/validation/RESULTS.md$", changed)]
commit <- paste0(git("rev-parse", "--short", "HEAD"), if (length(changed)) " (with uncommitted changes)" else "")
desc <- read.dcf("DESCRIPTION", fields = "Version")[1, 1]
md <- c(
  "# Validation results", "",
  sprintf("The validation scripts in this directory, run on malariafwd %s at commit %s on %s UTC.", desc, commit,
          format(Sys.time(), "%Y-%m-%d %H:%M", tz = "UTC")),
  sprintf("%s.", if (all(res$ok)) "Every script passed" else "NOT every script passed"), "",
  "| | version |", "|---|---|",
  sprintf("| R | %s |", paste(R.version$major, R.version$minor, sep = ".")),
  sprintf("| SLiM | %s |", v("SLiM")),
  sprintf("| Python | %s |", v("Python")),
  sprintf("| tskit | %s |", v("tskit")),
  sprintf("| numpy | %s |", v("numpy")),
  sprintf("| platform | %s |", R.version$platform), "",
  "| script | result | minutes |", "|---|---|---|",
  sprintf("| `%s` | %s | %.1f |", res$script, res$result, res$minutes), "",
  "`dev/validation-results.R`, in the package's source repository, writes this file.")
writeLines(md, file.path("inst", "validation", "RESULTS.md"))
cat(md, sep = "\n")
if (!all(res$ok)) quit(status = 1)
