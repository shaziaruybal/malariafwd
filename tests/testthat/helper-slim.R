## Shared test helpers. testthat sources helper*.R before any test file, so anything defined
## here is visible to all of them -- unlike a definition inside one test file, which is not.

skip_if_no_slim <- function() skip_if_not(slim_available(), "SLiM not on PATH")

## A stand-in for SLiM that reports `version` and otherwise runs `body`, for the paths that need
## SLiM to fail, hang or be too old. It is an sh script, so the tests that use it skip on Windows.
fake_slim <- function(body = character(), version = "5.2") {
  fake <- tempfile("fakeslim")
  writeLines(c("#!/bin/sh", paste0('if [ "$1" = "-v" ]; then echo "SLiM version ', version, '"; exit 0; fi'), body),
             fake)
  Sys.chmod(fake, "755")
  fake
}

## For the tests that read a tree sequence with Python: they skip, rather than error, without it.
skip_if_no_tskit <- function() {
  skip_if_not_installed("reticulate", "1.41.0")
  skip_if_not(tskit_available(), "Python with tskit (>= 1.0) and numpy not available")
}
