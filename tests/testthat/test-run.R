# These tests actually run SLiM, so they skip when it (or tskit) is unavailable.

test_that("run_sim produces the expected outputs", {
  skip_if_not(slim_available(), "SLiM not on PATH")
  od <- file.path(tempdir(), "mfwd_test_run")
  run <- run_sim(vivax_params(N_HOSTS = 120, N_MOSQ = 240, BITING_RATE = 0.12,
                              LAMBDA_SPOROZOITES = 2.0, SIM_DAYS = 365),
                 seed = 1, outdir = od)
  expect_s3_class(run, "vivax_run")
  expect_true(all(file.exists(file.path(od,
    c("events.tsv", "state.tsv", "plasmo.trees", "params.tsv", "model.slim", "run.log")))))
  expect_gt(nrow(load_events(run)), 0)
})

test_that("state.tsv has one row per host and per mosquito on every tick it covers", {
  skip_if_not(slim_available(), "SLiM not on PATH")
  run <- run_sim(vivax_params(N_HOSTS = 30, N_MOSQ = 60, SIM_DAYS = 40),
                 seed = 4, outdir = file.path(tempdir(), "mfwd_test_state"))
  st <- load_state(run)
  per_tick <- st[, .N, by = tick]
  expect_equal(per_tick$tick, seq_len(max(st$tick)))
  expect_true(all(per_tick$N == 30 + 60))
})

test_that("hypnozoite dynamics match alpha/(alpha+mu) [transmission-invariant]", {
  skip_if_not(slim_available(), "SLiM not on PATH")
  run <- run_sim(vivax_params(N_HOSTS = 120, N_MOSQ = 240, BITING_RATE = 0.12,
                              LAMBDA_SPOROZOITES = 2.0, SIM_DAYS = 365),
                 seed = 2, outdir = file.path(tempdir(), "mfwd_test_hyp"))
  ev <- load_events(run)
  act  <- ev[event_type == "hypnozoite_activation", sum(n_hypnozoites)]
  died <- ev[event_type == "hypnozoite_clearance",  sum(n_hypnozoites)]
  skip_if(act + died < 200, "too few resolved hypnozoites to test")
  expect_equal(act / (act + died), (1/170) / (1/170 + 1/442), tolerance = 0.06)
})

test_that("relapse-to-origin IBD is clonal (= 1)", {
  skip_if_not(slim_available(), "SLiM not on PATH")
  skip_if_not(requireNamespace("reticulate", quietly = TRUE), "reticulate not available")
  run <- run_sim(vivax_params(N_HOSTS = 120, N_MOSQ = 240, BITING_RATE = 0.12,
                              LAMBDA_SPOROZOITES = 2.0, SIM_DAYS = 365),
                 seed = 3, outdir = file.path(tempdir(), "mfwd_test_ibd"))
  ev <- load_events(run)
  cs <- classify_cases(ev)
  idx <- tryCatch(ibd_index(ibd_pairs(run, ids = cs$pedigree_id)),
                  error = function(e) skip(paste("tskit unavailable:", conditionMessage(e))))
  roi <- ibd_relapse_origin(ev, idx)
  skip_if(nrow(roi) < 10, "too few relapse-origin pairs")
  expect_equal(median(roi$ibd), 1, tolerance = 0.001)
})

test_that("run_sim() works with an apostrophe and a space in the output path and SIM_ID", {
  skip_if_not(slim_available(), "SLiM not on PATH")
  od <- file.path(tempdir(), "o'neil runs", "a")
  run <- run_sim(vivax_params(N_HOSTS = 20, N_MOSQ = 40, SIM_DAYS = 5, SIM_ID = "Bob's sweep"),
                 seed = 1, outdir = od)
  expect_true(file.exists(file.path(od, "events.tsv")))
  expect_equal(load_params(run)$SIM_ID, "Bob's sweep")
})

test_that("run_sim() stops when SLiM exits without writing the run's outputs", {
  skip_if_not(slim_available(), "SLiM not on PATH")
  tiny <- tempfile(fileext = ".slim")
  writeLines(c("initialize() { }", "1 early() { sim.simulationFinished(); }"), tiny)
  expect_error(run_sim(vivax_params(), seed = 1, outdir = tempfile(), model = tiny),
               "did not write")
})

test_that("blood-stage infections last a geometric time with mean MEAN_BLOOD_STAGE_DURATION days", {
  skip_if_not(slim_available(), "SLiM not on PATH")
  # founders only: they carry no hypnozoites, and a single mosquito almost never bites
  run <- run_sim(vivax_params(N_HOSTS = 1000, N_FOUNDERS = 1000, N_FOUNDER_HOSTS = 1000,
                              N_FOUNDER_LINEAGES = 1000, N_MOSQ = 1, BITING_RATE = 1e-6,
                              MEAN_BLOOD_STAGE_DURATION = 2, SIM_DAYS = 40),
                 seed = 6, outdir = file.path(tempdir(), "mfwd_test_clearance"))
  ev <- load_events(run)
  d <- merge(ev[event_type == "vh_inoculation_seed", .(pedigree_id, t0 = tick)],
             ev[event_type == "human_clearance", .(pedigree_id, t1 = tick)], by = "pedigree_id")
  expect_gt(nrow(d), 900)
  expect_lt(abs(mean(d$t1 - d$t0) - 2), 0.15)   # a daily probability of 1 - exp(-1/2) gives 2.54
})
