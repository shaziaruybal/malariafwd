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

test_that("a non-ASCII SIM_ID reaches SLiM unchanged, in a C locale too", {
  skip_if_not(slim_available(), "SLiM not on PATH")
  old <- Sys.getlocale("LC_CTYPE")
  on.exit(Sys.setlocale("LC_CTYPE", old), add = TRUE)
  skip_if(Sys.setlocale("LC_CTYPE", "C") == "", "cannot switch to the C locale")
  ## marked UTF-8, as read.csv(encoding = "UTF-8") or data.table::fread() give it
  id <- intToUtf8(c(0x43, 0xf4, 0x74, 0x65, 0x20, 0x65e5, 0x672c))
  od <- tempfile("mfwd_utf8_")
  run_sim(vivax_params(N_HOSTS = 20, N_MOSQ = 40, SIM_DAYS = 3, SIM_ID = id), seed = 1, outdir = od)
  for (f in c("params.tsv", "cmdline.txt")) {
    txt <- rawToChar(readBin(file.path(od, f), "raw", file.size(file.path(od, f))))
    expect_true(grepl(id, txt, fixed = TRUE, useBytes = TRUE), label = f)
  }
  ## and SLiM's own output reaches run.log whole, its closing banner included
  log <- rawToChar(readBin(file.path(od, "run.log"), "raw", file.size(file.path(od, "run.log"))))
  hourglass <- intToUtf8(c(0x231B, 0xFE0F))
  expect_true(grepl(paste0(hourglass, " SIMULATION FINISHED! ", hourglass), log, fixed = TRUE, useBytes = TRUE))
})

test_that("a set at the edge of the rounding allowance on the caps runs", {
  skip_if_not(slim_available(), "SLiM not on PATH")
  ## 333333.33333366667 x 3 is accepted at 1e6 x (1 + 1e-12); SLiM gets the 15-digit 333333.333333667
  p <- vivax_params(N_HOSTS = 20, N_MOSQ = 40, SIM_DAYS = 3, NU_HYPNOZOITE_BATCH = 333333.33333366667,
                    NU_HYPNOZOITE_BATCH_BOOST = 3, RESURGENCE_TRIGGER_TICK = 2L)
  expect_no_error(run_sim(p, seed = 1, outdir = tempfile("mfwd_edge_")))
})

test_that("run_grid() gives each warning from its cells once, after its FAILED summary", {
  skip_if_not(slim_available(), "SLiM not on PATH")
  ## R keeps only the first 50 warnings, so a set that warns in every cell would crowd it out
  p <- suppressWarnings(vivax_params(N_HOSTS = 20, N_MOSQ = 40, SIM_DAYS = 3, BITING_RATE = 2))
  q <- vivax_params(N_HOSTS = 20, N_MOSQ = 40, SIM_DAYS = 3)
  q$N_MOSQ <- 40.5                                   # refused in each of its cells
  w <- character()
  g <- withCallingHandlers(
    run_grid(list(warns = p, broken = q), seeds = 1:3, outdir_base = tempfile("mfwd_grid_")),
    warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  expect_identical(g$ok, g$set == "warns")
  expect_length(w, 2)
  expect_match(w[1], "3 of 6 grid cell\\(s\\) FAILED")
  expect_match(w[2], "bites at most once a day.* \\[3 cells\\]$")
})

test_that("under options(warn = 2) a warning fails its cell, as R makes it an error", {
  skip_if_not(slim_available(), "SLiM not on PATH")
  p <- suppressWarnings(vivax_params(N_HOSTS = 20, N_MOSQ = 40, SIM_DAYS = 3, BITING_RATE = 2))
  op <- options(warn = 2); on.exit(options(op), add = TRUE)
  expect_error(run_grid(p, seeds = 1:2, outdir_base = tempfile("mfwd_grid2_")),
               "2 of 2 grid cell\\(s\\) FAILED.*converted from warning")
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

test_that("no blood-stage infection is forced to clear at a fixed age", {
  skip_if_not(slim_available(), "SLiM not on PATH")
  run <- run_sim(vivax_params(N_HOSTS = 200, N_FOUNDERS = 200, N_FOUNDER_HOSTS = 200,
                              N_FOUNDER_LINEAGES = 200, N_MOSQ = 1, BITING_RATE = 1e-6,
                              MEAN_BLOOD_STAGE_DURATION = 1000, SIM_DAYS = 200),
                 seed = 7, outdir = file.path(tempdir(), "mfwd_test_nocap"))
  ev <- load_events(run)
  d <- merge(ev[event_type == "vh_inoculation_seed", .(pedigree_id, t0 = tick)],
             ev[event_type == "human_clearance", .(pedigree_id, t1 = tick)], by = "pedigree_id")
  expect_lt(sum(d$t1 - d$t0 == 181), 5)   # a 180-day cap cleared about 83% of founders on day 181
})

test_that("latent hypnozoites die at MU_HYPNOZOITE, so a temperate batch relapses less", {
  skip_if_not(slim_available(), "SLiM not on PATH")
  ## one pulse of imports, no onward transmission, long enough for every batch to resolve
  run <- run_sim(vivax_params(N_HOSTS = 300, N_MOSQ = 10, BITING_RATE = 1e-6, SIM_DAYS = 1500,
                              N_FOUNDERS = 0, IMPORT_TICK = 2, N_IMPORTS = 300,
                              K_LATENCY = 2L, DELTA_LATENCY = 1 / 100, CHROMOSOMES = "1"),
                 seed = 3, outdir = file.path(tempdir(), "mfwd_test_latency"))
  ev <- load_events(run)
  act  <- ev[event_type == "hypnozoite_activation", sum(n_hypnozoites)]
  died <- ev[event_type == "hypnozoite_clearance",  sum(n_hypnozoites)]
  a <- 1 / 170; mu <- 1 / 442; d <- 1 / 100
  ## each latent day a hypnozoite dies w.p. mu, then the batch moves on w.p. d (two stages)
  stage <- d * (1 - mu) / (1 - (1 - mu) * (1 - d))
  expect_gt(act + died, 600)
  expect_equal(act / (act + died), stage^2 * a / (a + mu), tolerance = 0.05 / 0.48)
})
