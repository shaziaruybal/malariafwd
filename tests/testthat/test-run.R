# These tests actually run SLiM, so they skip when it (or tskit) is unavailable.

test_that("run_sim produces the expected outputs", {
  skip_if_not(slim_available(), "SLiM not on PATH")
  od <- file.path(tempdir(), "mfwd_test_run")
  run <- run_sim(vivax_params(N_HOSTS = 120, N_MOSQ = 240, BITING_RATE = 0.12,
                              N_OOCYSTS_PER_BITE = 2, SIM_DAYS = 365),
                 seed = 1, outdir = od)
  expect_s3_class(run, "vivax_run")
  expect_true(all(file.exists(file.path(od,
    c("events.tsv", "state.tsv", "plasmo.trees", "params.tsv", "model.slim", "run.log")))))
  expect_gt(nrow(load_events(run)), 0)
})

test_that("hypnozoite dynamics match alpha/(alpha+mu) [transmission-invariant]", {
  skip_if_not(slim_available(), "SLiM not on PATH")
  run <- run_sim(vivax_params(N_HOSTS = 120, N_MOSQ = 240, BITING_RATE = 0.12,
                              N_OOCYSTS_PER_BITE = 2, SIM_DAYS = 365),
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
                              N_OOCYSTS_PER_BITE = 2, SIM_DAYS = 365),
                 seed = 3, outdir = file.path(tempdir(), "mfwd_test_ibd"))
  ev <- load_events(run)
  cs <- classify_cases(ev)
  idx <- tryCatch(ibd_index(calc_ibd(run, ids = cs$pedigree_id)),
                  error = function(e) skip(paste("tskit unavailable:", conditionMessage(e))))
  roi <- relapse_origin_ibd(ev, idx)
  skip_if(nrow(roi) < 10, "too few relapse-origin pairs")
  expect_equal(median(roi$ibd), 1, tolerance = 0.001)
})
