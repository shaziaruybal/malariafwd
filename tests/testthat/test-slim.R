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
    processx::run(slim_path(), c("-s", "1", "-d", paste0("SIM_DAYS=", days), "-d", "N_HOSTS=20", "-d", "N_MOSQ=20",
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
           "K_OOCYSTS=-1.0" = "K_OOCYSTS must be 0 \\(Poisson\\)",
           ## at 0 the gate would refuse every inoculation, an empty host's included
           "MAX_MOI=0" = "MAX_MOI must be at least 1")
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
  ## ... or whose only window of imports ends on day 1
  r <- slim("-d", "N_FOUNDERS=0", "-d", "IMPORT_START=0", "-d", "IMPORT_END=1", "-d", "N_IMPORTS_PER_TICK=5")
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

test_that("a direct run refuses what validate_params() refuses, and takes either spelling of a number", {
  skip_if_not(slim_available(), "SLiM not on PATH")
  ## a 5-day, 20-host run; a setting given here replaces the default one
  slim <- function(...) {
    set <- c(SIM_DAYS = "5", N_HOSTS = "20", N_MOSQ = "20")
    a <- c(...); set[sub("=.*", "", a)] <- sub("^[^=]*=", "", a)
    out <- tempfile("mfwd_direct_"); dir.create(out)
    r <- processx::run(slim_path(), c("-s", "1", as.vector(rbind("-d", c(paste0(names(set), "=", set),
                                      paste0("OUTDIR='", out, "'")))), slim_model()), error_on_status = FALSE)
    r$outdir <- out; r
  }
  bad <- c("RADCURE_MODE='Person'" = "RADCURE_MODE must be 'person' or 'hypnozoite'",
           "P_TRANS_VH=1.5" = "P_TRANS_VH must be in \\[0, 1\\]",
           "ALPHA_ACTIVATE=NAN" = "ALPHA_ACTIVATE must be in \\[0, 1\\]",
           "N_MOSQ=0" = "N_MOSQ must be > 0",
           "BITING_RATE=INF" = "BITING_RATE must be > 0 and finite",
           "EIP_DAYS=-1" = "EIP_DAYS must be >= 0",
           "N_FOUNDERS=21" = "N_FOUNDER_HOSTS \\(21\\) must be at most N_HOSTS \\(20\\)",
           "N_HOSTS=20.5" = "N_HOSTS must be a whole number",
           "MOSQ_LIFESPAN=0.5" = "MOSQ_LIFESPAN must be >= 1",
           "MEAN_BLOOD_STAGE_DURATION=0.5" = "MEAN_BLOOD_STAGE_DURATION must be >= 1",
           "MU_HYPNOZOITE=0.0" = "MU_HYPNOZOITE must be in \\(0, 1\\]",
           "VERBOSE='F'" = "VERBOSE must be T or F",
           "SIM_ID=3" = "SIM_ID must be a string",
           "N_MOSQ=c(10, 20)" = "N_MOSQ must be a single value",
           "CHROMOSOMES='1,1'" = "CHROMOSOMES lists a chromosome twice",
           "CHROMOSOMES='3-1'" = "CHROMOSOMES has a range that runs high-low",
           "MDA_TICK=1" = "MDA_TICK must not be 1 in a run with founders")
  for (arg in names(bad)) {
    r <- slim(arg)
    expect_true(r$status != 0, label = arg)
    expect_match(paste(r$stdout, r$stderr), bad[[arg]], label = arg)
  }
  r <- slim("BLOODCURE_COVERAGE=0.5", "TREATMENT_START=10", "TREATMENT_END=5")
  expect_match(paste(r$stdout, r$stderr), "case-management window is closed")
  ## a number in either spelling gives the same run, through latency, an MDA round in the
  ## per-hypnozoite mode and a pulse of imports
  on <- c("K_LATENCY=1", "DELTA_LATENCY=0.5", "RADCURE_MODE='hypnozoite'", "MDA_TICK=3", "MDA_COVERAGE=1.0",
          "IMPORT_TICK=2")
  a <- slim(on, "ALPHA_ACTIVATE=1", "MU_HYPNOZOITE=1", "P_RAD_SURVIVAL=0", "N_IMPORTS=2.0", "N_HOSTS=20.0",
            "SIM_DAYS=5.0", "SEED=5.0")
  b <- slim(on, "ALPHA_ACTIVATE=1.0", "MU_HYPNOZOITE=1.0", "P_RAD_SURVIVAL=0.0", "N_IMPORTS=2", "N_HOSTS=20",
            "SIM_DAYS=5", "SEED=5")
  expect_equal(c(a$status, b$status), c(0L, 0L))
  for (f in c("events.tsv", "state.tsv"))
    expect_identical(unname(tools::md5sum(file.path(a$outdir, f))), unname(tools::md5sum(file.path(b$outdir, f))),
                     label = f)
})
