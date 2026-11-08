## The one Python check (R/python.R) and check_setup(). The Python side is mocked through the
## package's thin wrappers, so most of these need neither Python nor SLiM.

## Mocks reticulate as seen by the check, until the calling test ends. `modules` gives each module's
## version (NA: it does not import); `start` is what py_config() gives, or does.
mock_py <- function(has = TRUE, version = "1.46.0", started = FALSE, declare = function() invisible(),
                    start = function() list(python = "/opt/py/bin/python", version = "3.12",
                                            forced = "RETICULATE_PYTHON", ephemeral = FALSE),
                    modules = c(tskit = "1.0.4", numpy = "2.0.0"),
                    env = parent.frame()) {
  local_mocked_bindings(.py_cached = function() NULL, .has_reticulate = function() has,
                        .reticulate_version = function() version, .py_started = function() started,
                        .py_declare = declare, .py_start = start,
                        .py_module_version = if (is.function(modules)) modules else function(m) as.character(unname(modules[m])),
                        .env = env)
}
problem <- function(modules = c("tskit", "numpy"))
  malariafwd:::.py_problem(malariafwd:::.python_status(), modules)

test_that("without reticulate the check says to install it, and tskit_available() is FALSE quietly", {
  mock_py(has = FALSE)
  expect_match(problem(), "through the reticulate package, which is not installed", fixed = TRUE)
  expect_match(problem(), "install.packages(\"reticulate\")", fixed = TRUE)
  expect_silent(expect_false(tskit_available()))
})

test_that("a reticulate older than 1.41.0 is named", {
  mock_py(version = "1.40.0")
  expect_match(problem(), "reticulate >= 1.41.0, and reticulate 1.40.0 is installed")
})

test_that("a Python that does not start is quoted, with what to do", {
  mock_py(start = function() stop("Python specified in RETICULATE_PYTHON (/x) does not exist\nmore"))
  old <- Sys.getenv("RETICULATE_PYTHON", NA)
  Sys.setenv(RETICULATE_PYTHON = "/x")
  on.exit(if (is.na(old)) Sys.unsetenv("RETICULATE_PYTHON") else Sys.setenv(RETICULATE_PYTHON = old), add = TRUE)
  p <- problem()
  expect_match(p, "could not start Python: Python specified in RETICULATE_PYTHON (/x) does not exist.", fixed = TRUE)
  expect_match(p, "Check RETICULATE_PYTHON (/x), then restart R.", fixed = TRUE)
  expect_false(grepl("more", p))
})

test_that("a chosen Python without the modules is named, with how it was chosen and how to fix it", {
  mock_py(modules = c(tskit = NA, numpy = NA))
  p <- problem()
  expect_match(p, "The Python reticulate is using, /opt/py/bin/python (Python 3.12, chosen by RETICULATE_PYTHON)", fixed = TRUE)
  expect_match(p, "cannot import tskit or numpy", fixed = TRUE)
  expect_match(p, "/opt/py/bin/python -m pip install \"tskit>=1.0\" \"numpy\"", fixed = TRUE)
  expect_match(p, "restart R")
})

test_that("a tskit below 1.0 is refused, naming the floor", {
  mock_py(modules = c(tskit = "0.5.6", numpy = "2.0.0"))
  expect_match(problem(), "has tskit 0.5.6; malariafwd needs tskit >= 1.0")
  ## a development version reads as its leading numbers
  mock_py(modules = c(tskit = "1.0.4.dev3", numpy = "2.0.0"))
  expect_null(problem())
})

test_that("reticulate's own environment gets the modules declared again, then read again", {
  n <- 0L
  mock_py(declare = function() { n <<- n + 1L; invisible() },
          start = function() list(python = "/cache/env/bin/python", version = "3.13", forced = "py_require()",
                                  ephemeral = TRUE),
          modules = function(m) if (n < 2L) NA_character_ else c(tskit = "1.0.4", numpy = "2.0.0")[[m]])
  st <- malariafwd:::.python_status()
  expect_identical(n, 2L)
  expect_true(st$managed)
  expect_null(malariafwd:::.py_problem(st))
})

test_that("reticulate's warnings while starting do not escape when every module imports", {
  mock_py(start = function() {
    warning("some Python packages are not installed")
    list(python = "/venv/bin/python", version = "3.12", forced = "VIRTUAL_ENV", ephemeral = FALSE)
  })
  expect_no_warning(st <- malariafwd:::.python_status())
  expect_null(malariafwd:::.py_problem(st))
})

test_that("with every module present the check passes, and a function asks only for what it uses", {
  mock_py()
  expect_true(tskit_available())
  expect_silent(malariafwd:::.require_python("ibd_pairs()"))
  mock_py(modules = c(tskit = "1.0.4", numpy = NA))
  expect_false(tskit_available())
  expect_silent(malariafwd:::.require_python("a_tskit_only()", "tskit"))
  expect_error(malariafwd:::.require_python("ibd_pairs()"), "^ibd_pairs\\(\\) needs Python with tskit \\(>= 1.0\\) and numpy\\. .*cannot import numpy")
})

test_that("each function that runs Python stops, naming itself, before it reads a file", {
  mock_py(has = FALSE)
  d <- file.path(tempdir(), "no_such_run")
  calls <- list("ibd_pairs()" = function() ibd_pairs(d), "ibd_segments()" = function() ibd_segments(d, 1, 2),
                "ibd_truth_pairs()" = function() ibd_truth_pairs(d), "arg_summary()" = function() arg_summary(d),
                "arg_tmrca()" = function() arg_tmrca(d, 1, 2))
  if (requireNamespace("ggplot2", quietly = TRUE))
    calls[["plot_ibd_pair()"]] <- function() plot_ibd_pair(d, 1, 2)
  for (fn in names(calls)) {
    e <- tryCatch(calls[[fn]](), error = identity)
    expect_s3_class(e, "malariafwd_no_python")
    expect_true(startsWith(conditionMessage(e), paste(fn, "needs Python")), label = fn)
  }
  ## a run with too few cases for the measure gives its empty row without Python
  r <- fake_run(ev_rows = "1\t1_1\tvh_inoculation_seed\t-1\t-1\t0\t5\tNA\tNA\t0\t0.0\thuman\t0\t0\t5\t-1\t-1\t0\t-1\t1\t-1")
  expect_identical(new_infection_diversity(r)$n, 0L)
})

test_that("check_setup() reports a missing SLiM and leaves Python unchecked when asked", {
  old <- Sys.getenv("MALARIAFWD_SLIM", NA)
  Sys.setenv(MALARIAFWD_SLIM = file.path(tempdir(), "no-slim-here"))
  on.exit(if (is.na(old)) Sys.unsetenv("MALARIAFWD_SLIM") else Sys.setenv(MALARIAFWD_SLIM = old), add = TRUE)
  s <- check_setup(python = FALSE)
  expect_s3_class(s, "malariafwd_setup")
  expect_identical(s$component, c("malariafwd", "R", "SLiM", "model", "reticulate", "Python", "tskit",
                                  "numpy", "ggplot2", "patchwork"))
  expect_false(s$ok[s$component == "SLiM"])
  expect_match(s$advice[s$component == "SLiM"], "MALARIAFWD_SLIM")
  expect_true(all(is.na(s$ok[s$component %in% c("Python", "tskit", "numpy")])))
  out <- capture.output(print(s))
  expect_length(out, nrow(s) + sum(!is.na(s$advice)) + 1L)
  expect_true(all(grepl("^[ -~]*$", out)))
  expect_match(out[length(out)], "problem")
})

test_that("check_setup() reports a SLiM that is too old", {
  skip_on_os("windows")
  fake <- tempfile("fakeslim51")
  writeLines(c("#!/bin/sh", 'if [ "$1" = "-v" ]; then echo "SLiM version 5.1, built Jan 1 2025"; exit 0; fi'), fake)
  Sys.chmod(fake, "755")
  old <- Sys.getenv("MALARIAFWD_SLIM", NA)
  Sys.setenv(MALARIAFWD_SLIM = fake)
  on.exit(if (is.na(old)) Sys.unsetenv("MALARIAFWD_SLIM") else Sys.setenv(MALARIAFWD_SLIM = old), add = TRUE)
  s <- check_setup(python = FALSE)
  expect_false(s$ok[s$component == "SLiM"])
  expect_identical(s$version[s$component == "SLiM"], "5.1")
  expect_match(s$advice[s$component == "SLiM"], "needs SLiM >= 5.2, found 5.1")
})

test_that("check_setup() finds this machine's SLiM, and its Python when tskit is there", {
  skip_if_no_slim()
  s <- check_setup(python = FALSE)
  expect_true(s$ok[s$component == "SLiM"])
  expect_identical(s$version[s$component == "SLiM"], slim_version())
  skip_if_not(tskit_available(), "Python with tskit is not available")
  s <- check_setup()
  expect_true(all(s$ok[s$component %in% c("reticulate", "Python", "tskit", "numpy")]))
  expect_true(numeric_version(sub("^([0-9.]+).*", "\\1", s$version[s$component == "tskit"])) >= "1.0")
})

test_that("skip_if_no_tskit() skips when Python with tskit is not there", {
  local_mocked_bindings(tskit_available = function() FALSE)
  expect_condition(skip_if_no_tskit(), class = "skip")
})
