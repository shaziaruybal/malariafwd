test_that("a SLiM older than 5.2 is refused, and an unreadable version is let through", {
  expect_false(malariafwd:::.slim_version_ok("5.1"))
  expect_false(malariafwd:::.slim_version_ok("4.3"))
  expect_true(malariafwd:::.slim_version_ok("5.2"))
  expect_true(malariafwd:::.slim_version_ok("6.0"))
  expect_true(malariafwd:::.slim_version_ok(NA_character_))
})

test_that("params.tsv records the seed and the SLiM and package versions of the run", {
  skip_if_not(slim_available(), "SLiM not on PATH")
  run <- run_sim(vivax_params(N_HOSTS = 20, N_MOSQ = 40, SIM_DAYS = 5), seed = 5,
                 outdir = file.path(tempdir(), "mfwd_test_versions"))
  p <- load_params(run)
  expect_equal(p$ACTUAL_SEED, "5")
  expect_equal(p$SLIM_VERSION, slim_version())
  expect_equal(p$MALARIAFWD_VERSION, as.character(utils::packageVersion("malariafwd")))
})

test_that("the model refuses values its draws cannot take when SLiM is run directly", {
  skip_if_not(slim_available(), "SLiM not on PATH")
  ## validate_params() refuses these too, but the model can be run with -d, bypassing R
  slim <- function(..., days = 5) {
    out <- tempfile("mfwd_direct_"); dir.create(out)
    processx::run(slim_path(), c("-s", "1", "-d", paste0("SIM_DAYS=", days), "-d", "N_HOSTS=10", "-d", "N_MOSQ=20",
                                 "-d", paste0("OUTDIR='", out, "'"), ..., slim_model()),
                  error_on_status = FALSE)
  }
  bad <- c("TREATMENT_DELAY=0.0" = "must be >= 1 and finite",
           "TREATMENT_DELAY=0.5" = "must be >= 1 and finite",
           "MEAN_RECRUD_DELAY=0.5" = "must be >= 1 and finite",
           "MEAN_RECRUD_DELAY=INF" = "must be >= 1 and finite",
           "K_OOCYSTS=0.0001" = "K_OOCYSTS must be 0 \\(Poisson\\) or at least 0.001",
           "INTERVENTION_EFFICACY=1.5" = "INTERVENTION_EFFICACY must be in \\[0, 1\\]",
           "IMPORT_RATE_BOOST=-1.0" = "must be >= 0",
           "SEASONAL_AMP=-0.5" = "SEASONAL_AMP must be >= 0",
           ## infinite or huge values pass a sign check, then hang a draw or run as something else
           "K_OOCYSTS=INF" = "K_OOCYSTS must be 0 \\(Poisson\\) or at least 0.001 and finite",
           "LAMBDA_OOCYSTS=1001.0" = "must be at most 1000",
           "IMPORT_RATE_PER_TICK=INF" = "at most 1e6 imports a day",
           "IMPORT_RATE_PER_TICK=1.0e18" = "at most 1e6 imports a day",
           "IMPORT_RATE_BOOST=INF" = "must be >= 0 and finite",
           "BITING_RATE_BOOST=INF" = "must be >= 0 and finite",
           "SEASONAL_AMP=INF" = "SEASONAL_AMP must be >= 0 and finite",
           "NU_HYPNOZOITE_BATCH=2.0e6" = "at most 1e6 hypnozoites per bite",
           "ALPHA_ACTIVATE_BOOST=NAN" = "ALPHA_ACTIVATE_BOOST must be >= 0",
           "ALPHA_ACTIVATE_BOOST=INF" = "ALPHA_ACTIVATE_BOOST must be >= 0 and finite",
           "K_OOCYSTS=-1.0" = "K_OOCYSTS must be 0 \\(Poisson\\)")
  for (arg in names(bad)) {
    r <- slim("-d", arg)
    expect_true(r$status != 0, label = arg)
    expect_match(paste(r$stdout, r$stderr), bad[[arg]], label = arg)
  }
  ## a run long enough to reach day 99999, which the defaults use for 'never'
  r <- slim(days = 100000)
  expect_true(r$status != 0)
  expect_match(paste(r$stdout, r$stderr), "reaches day 99999")
  ## a founder-free run whose only import falls on day 0, before imports can arrive
  r <- slim("-d", "N_FOUNDERS=0", "-d", "IMPORT_TICK=0", "-d", "N_IMPORTS=5")
  expect_true(r$status != 0)
  expect_match(paste(r$stdout, r$stderr), "nothing can start this run")
  ## the integer spelling that -d X=7 produces is accepted at the start...
  expect_equal(slim("-d", "TREATMENT_DELAY=7", "-d", "MEAN_RECRUD_DELAY=4")$status, 0L)
  ## ...and each recrudescence shape draws its delay from it
  for (shape in c("geometric", "poisson", "fixed")) {
    r <- slim("-d", "MEAN_RECRUD_DELAY=4", "-d", paste0("RECRUD_DELAY_DIST='", shape, "'"),
              "-d", "IMPORT_TICK=3", "-d", "N_IMPORTS=10", "-d", "TREATMENT_DELAY=1",
              "-d", "BLOODCURE_COVERAGE=1.0", "-d", "P_BLOODCURE_SURVIVAL=1.0",
              "-d", "P_RECRUD_GIVEN_FAIL=1.0")
    expect_equal(r$status, 0L, label = paste0("a '", shape, "' delay given as an integer"))
  }
})
