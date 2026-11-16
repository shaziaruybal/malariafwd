## Writes the generated calibration table of MODEL_DESCRIPTION.md section 13 from the numbers
## inst/validation/calibrate_nu.R measures, so nobody retypes them. Run from the package root:
##
##   Rscript dev/calibration-table.R                      # runs the calibration (about 6 min), then fills
##   Rscript dev/calibration-table.R <dir>                 # fills from a run's tables in <dir>
##
## <dir> holds calibrate_nu.tsv and calibrate_nu_targets.tsv, which calibrate_nu.R writes into
## MALARIAFWD_CALIBRATE_OUT. dev/validation-results.R fills the same block on a release run.
## The paragraph BELOW the block reads the table: check it still says what the numbers do.

MARK <- "calibration table"
MD <- "MODEL_DESCRIPTION.md"

## one markdown table of the per-nu means, and one sentence of what the targets say about them
calibration_md <- function(dir) {
  m <- utils::read.delim(file.path(dir, "calibrate_nu.tsv"), check.names = FALSE)
  tg <- utils::read.delim(file.path(dir, "calibrate_nu_targets.tsv"), check.names = FALSE)
  rng <- function(metric) {
    r <- tg[tg$metric == metric, ]
    sprintf("%g-%g", r$lo, r$hi)
  }
  shipped <- malariafwd::vivax_params()$NU_HYPNOZOITE_BATCH
  nu <- ifelse(m$nu == shipped, sprintf("**%g**", m$nu), sprintf("%g", m$nu))
  two <- function(a, b) sprintf("%.2f / %.2f", a, b)
  rows <- sprintf("| %s | %.2f | %s | %s | %.2f | %s | %.2f | %.3f |",
                  nu, m$prevalence,
                  two(m$relapses_per_py, m$relapses_per_py_unmasked),
                  two(m$molfob_per_py, m$molfob_per_py_unmasked),
                  m$moi_median,
                  two(m$relapse_frac, m$relapse_frac_unmasked),
                  m$eir_per_py, m$multi_bites)
  inside <- function(flag) {
    v <- m$nu[flag]
    if (!length(v)) "none of the values tried"
    else sprintf("\u03bd = %s", paste(sprintf("%g", v), collapse = ", "))
  }
  c(sprintf("| ν | prevalence | relapses per person-year (%s): all / masked left out | molFOB (%s): all / masked left out | median MOI (%s) | relapse fraction (%s): all / masked left out | EIR | bites with more than one genotype |",
            rng("relapses_per_py"), rng("molfob_per_py"), rng("moi_median"), rng("relapse_frac")),
    "|---|---|---|---|---|---|---|---|",
    rows,
    "",
    sprintf("Inside all four targets, counting every relapse: %s. With the masked relapses left out of both readings: %s.",
            inside(as.logical(m$in_all)), inside(as.logical(m$in_all_unmasked))))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args)) {
  dir <- args[1]
} else {
  dir <- file.path(tempdir(), "calibrate_out")
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  cat("running inst/validation/calibrate_nu.R; about 6 minutes\n")
  status <- system2(file.path(R.home("bin"), "Rscript"),
                    c("-e", shQuote("source('inst/validation/calibrate_nu.R')")),
                    env = paste0("MALARIAFWD_CALIBRATE_OUT=", dir))
  if (status != 0L) stop("calibrate_nu.R did not pass; the table is unchanged", call. = FALSE)
}
stopifnot(file.exists(file.path(dir, "calibrate_nu.tsv")))
suppressMessages(pkgload::load_all(".", quiet = TRUE))

txt <- readLines(MD, warn = FALSE, encoding = "UTF-8")
begin <- grep(paste0("^<!-- BEGIN ", MARK), txt)
end <- grep(paste0("^<!-- END ", MARK, " -->$"), txt)
if (length(begin) != 1L || length(end) != 1L || end < begin)
  stop(MD, " must hold one BEGIN and one END ", MARK, " marker", call. = FALSE)
txt <- c(txt[seq_len(begin)], "", calibration_md(dir), "", txt[end:length(txt)])
writeLines(enc2utf8(txt), MD, useBytes = TRUE)
cat("wrote section 13's calibration table into", MD, "\n")
