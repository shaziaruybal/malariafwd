test_that("vivax_params defaults stay in sync with the SLiM model", {
  model <- default_model()
  skip_if(!nzchar(model) || !file.exists(model), "model template not found")
  src <- readLines(model, warn = FALSE)
  model_names <- unique(sub('.*defineDefault\\("([A-Z_]+)".*', "\\1",
                            grep('defineDefault\\("[A-Z_]+"', src, value = TRUE)))
  r_names <- names(malariafwd:::.malariafwd_defaults)
  # SEED and OUTDIR are run-level (set by run_sim), not part of vivax_params().
  expect_setequal(setdiff(model_names, c("SEED", "OUTDIR")), r_names)
})

test_that("validation catches bad input", {
  expect_error(vivax_params(GENOME_MODE = "xyz"))
  expect_error(vivax_params(N_HOSTS = 5, N_FOUNDERS = 10))
  expect_error(vivax_params(FOO = 1))
  expect_error(vivax_params(P_TRANS_VH = 2))
  expect_error(vivax_params(SIM_DAYS = -1))
})

test_that("as_slim_args quotes strings, formats logicals, passes numbers", {
  a <- as_slim_args(vivax_params(GENOME_MODE = "core", VERBOSE = TRUE, N_HOSTS = 200))
  expect_true("GENOME_MODE='core'" %in% a)
  expect_true("VERBOSE=T" %in% a)
  expect_true("RETAIN_FULL_HISTORY=T" %in% a)
  expect_true("N_HOSTS=200" %in% a)
  expect_true(all(a[seq(1, length(a), 2)] == "-d"))  # every other token is the -d flag
})

test_that("print.vivax_params reports non-defaults", {
  expect_output(print(vivax_params(N_HOSTS = 300)), "N_HOSTS")
})
