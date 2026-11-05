## Values that validate_params() accepts must run to the end: the model may not stop part-way
## with an Eidos error. Each case below stopped the run before the guard it tests.

edge_run <- function(tag, ...) {
  run_sim(vivax_params(N_HOSTS = 40, N_MOSQ = 100, CHROMOSOMES = "1", ...), seed = 3,
          outdir = file.path(tempdir(), paste0("mfwd_edge_", tag)))
}

test_that("SUBPATENT_INFECTIVITY = 0 runs when a host carries only dormant batches", {
  skip_if_no_slim()
  run <- edge_run("subpat0", SIM_DAYS = 60, SUBPATENT_INFECTIVITY = 0)
  expect_equal(max(load_state(run)$tick), 60)
})

test_that("the hypnozoite death probability stays a probability when alpha + mu = 1", {
  skip_if_no_slim()
  ## 0.1 / (1 - 0.9) is 1.0000000000000002 in floating point
  run <- edge_run("alphamu", SIM_DAYS = 60, ALPHA_ACTIVATE = 0.9, MU_HYPNOZOITE = 0.1)
  expect_equal(max(load_state(run)$tick), 60)
})

test_that("an import rate that reaches zero means no imports, not an error", {
  skip_if_no_slim()
  ## the border closes after the trigger...
  run <- edge_run("impboost0", SIM_DAYS = 30, IMPORT_START = 5, IMPORT_END = 30,
                  IMPORT_RATE_PER_TICK = 0.2, RESURGENCE_TRIGGER_TICK = 10, IMPORT_RATE_BOOST = 0)
  ev <- load_events(run)
  expect_equal(max(load_state(run)$tick), 30)
  expect_equal(nrow(ev[event_type == "vh_inoculation_import" & tick > 10]), 0L)
  ## ...and a seasonal amplitude above 1 floors the seasonal factor at 0 for part of the year
  run <- edge_run("impseason", SIM_DAYS = 140, IMPORT_START = 5, IMPORT_END = 140,
                  IMPORT_RATE_PER_TICK = 0.2, SEASONAL_AMP = 1.5)
  expect_equal(max(load_state(run)$tick), 140)
})

test_that("a small K_OOCYSTS runs, although its gamma draws can underflow to zero", {
  skip_if_no_slim()
  ## what uptake_calibrate(3, var = 30) returns
  run <- edge_run("kooc", SIM_DAYS = 10, LAMBDA_OOCYSTS = 0.009644, K_OOCYSTS = 0.001)
  expect_equal(max(load_state(run)$tick), 10)
})
