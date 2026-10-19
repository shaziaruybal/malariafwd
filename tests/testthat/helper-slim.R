## Shared test helpers. testthat sources helper*.R before any test file, so anything defined
## here is visible to all of them -- unlike a definition inside one test file, which is not.

skip_if_no_slim <- function() skip_if_not(slim_available(), "SLiM not on PATH")
