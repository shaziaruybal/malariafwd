## Values that validate_params() accepts must run to the end: the model may not stop part-way
## with an Eidos error. Without the guard it tests, each case below stops the run, and each
## checks that the run reached the state the guard handles, so a lucky seed cannot pass it.

edge_run <- function(tag, ...) {
  run_sim(vivax_params(N_HOSTS = 40, N_MOSQ = 100, CHROMOSOMES = "1", ...), seed = 3,
          outdir = tempfile(paste0("mfwd_edge_", tag)))
}

test_that("SUBPATENT_INFECTIVITY = 0 runs when a host carries only dormant batches", {
  skip_if_no_slim()
  run <- edge_run("subpat0", SIM_DAYS = 120, SUBPATENT_INFECTIVITY = 0)
  st <- load_state(run)
  expect_equal(max(st$tick), 120)
  ## a host with nothing but dormant batches is infected and not infectious
  dormant_only <- st[entity_type == "human" & moi_active == 0 & moi_subpatent == 0 & moi_dormant > 0]
  expect_gt(nrow(dormant_only), 0)
  expect_true(all(dormant_only$tag == 0))
})

test_that("the hypnozoite death probability stays a probability when alpha + mu = 1", {
  skip_if_no_slim()
  ## 0.1 / (1 - 0.9) is 1.0000000000000002 in floating point
  run <- edge_run("alphamu", SIM_DAYS = 120, ALPHA_ACTIVATE = 0.9, MU_HYPNOZOITE = 0.1)
  expect_equal(max(load_state(run)$tick), 120)
  ## each batch is gone after its first dormant day: what does not activate dies
  clear <- load_events(run)[event_type == "hypnozoite_clearance"]
  expect_gt(nrow(clear), 0)
  expect_equal(max(table(clear$pedigree_id)), 1)
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

test_that("local vector control scales biting, not imports, and gives no course", {
  skip_if_no_slim()
  ## full vector control from day 10: no mosquito bites from then on, while imports, infections
  ## acquired elsewhere, keep arriving
  run <- edge_run("vc1", SIM_DAYS = 60, IMPORT_START = 5, IMPORT_END = 60, IMPORT_RATE_PER_TICK = 1,
                  INTERVENTION_START = 10, INTERVENTION_EFFICACY = 1)
  ev <- load_events(run)
  expect_equal(max(load_state(run)$tick), 60)
  bites <- c("vh_inoculation", "vh_inoculation_no_take", "hv_uptake_clonal", "hv_uptake_recombinant")
  expect_gt(nrow(ev[event_type %in% bites & tick < 10]), 0L)
  expect_equal(nrow(ev[event_type %in% bites & tick >= 10]), 0L)
  expect_gt(nrow(ev[event_type == "vh_inoculation_import" & tick > 10]), 40L)
  ## the imports are detected cases, but treatment comes from a coverage, and there is none (D30)
  expect_equal(nrow(ev[event_type %in% c("bloodcure_given", "radcure_given")]), 0L)
})

test_that("a small K_OOCYSTS runs, although its gamma draws can underflow to zero", {
  skip_if_no_slim()
  ## a dispersed uptake at the K_OOCYSTS floor (a mean of 4.08 oocysts)
  run <- edge_run("kooc", SIM_DAYS = 10, LAMBDA_OOCYSTS = 0.009644, K_OOCYSTS = 0.001)
  expect_equal(max(load_state(run)$tick), 10)
})
