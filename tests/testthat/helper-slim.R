## Shared test helpers. testthat sources helper*.R before any test file, so anything defined
## here is visible to all of them -- unlike a definition inside one test file, which is not.

skip_if_no_slim <- function() skip_if_not(slim_available(), "SLiM not on PATH")

## For the tests that read a tree sequence with Python: they skip, rather than error, without it.
skip_if_no_tskit <- function() {
  skip_if_not_installed("reticulate", "1.41.0")
  skip_if_not(tskit_available(), "Python with tskit (>= 1.0), pandas and numpy not available")
}
