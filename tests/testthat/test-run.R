# These tests actually run SLiM, so they skip when it (or tskit) is unavailable.

test_that("run_sim produces the expected outputs", {
  skip_if_not(slim_available(), "SLiM not on PATH")
  od <- tempfile("mfwd_test_run")
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
                 seed = 4, outdir = tempfile("mfwd_test_state"))
  st <- load_state(run, entity = "all")
  per_tick <- st[, .N, by = tick]
  expect_equal(per_tick$tick, seq_len(max(st$tick)))
  expect_true(all(per_tick$N == 30 + 60))
  expect_equal(nrow(load_state(run)), 30 * max(st$tick))
  expect_equal(nrow(load_state(run, entity = "mosquito")), 60 * max(st$tick))
})

test_that("hypnozoite dynamics match alpha/(alpha+mu) [transmission-invariant]", {
  skip_if_not(slim_available(), "SLiM not on PATH")
  run <- run_sim(vivax_params(N_HOSTS = 120, N_MOSQ = 240, BITING_RATE = 0.12,
                              LAMBDA_SPOROZOITES = 2.0, SIM_DAYS = 365),
                 seed = 2, outdir = tempfile("mfwd_test_hyp"))
  ev <- load_events(run)
  act  <- ev[event_type == "hypnozoite_activation", sum(n_hypnozoites)]
  died <- ev[event_type == "hypnozoite_clearance",  sum(n_hypnozoites)]
  skip_if(act + died < 200, "too few resolved hypnozoites to test")
  expect_equal(act / (act + died), (1/170) / (1/170 + 1/442), tolerance = 0.06)
})

test_that("relapse-to-origin IBD is clonal (= 1)", {
  skip_if_no_slim(); skip_if_no_tskit()
  run <- run_sim(vivax_params(N_HOSTS = 120, N_MOSQ = 240, BITING_RATE = 0.12,
                              LAMBDA_SPOROZOITES = 2.0, SIM_DAYS = 365),
                 seed = 3, outdir = tempfile("mfwd_test_ibd"))
  ev <- load_events(run)
  cs <- classify_cases(ev)
  idx <- ibd_index(ibd_pairs(run, ids = cs$pedigree_id))
  roi <- ibd_relapse_origin(ev, idx)
  skip_if(nrow(roi) < 10, "too few relapse-origin pairs")
  expect_equal(median(roi$ibd), 1, tolerance = 0.001)
})

test_that("run_sim() works with an apostrophe and a space in the output path and SIM_ID", {
  skip_if_not(slim_available(), "SLiM not on PATH")
  od <- file.path(tempfile("o'neil runs "), "a")
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
  ## 333333.333333666 x 3 = 1e6 x (1 + 9.98e-13) is accepted at 1e6 x (1 + 1e-12); it is the
  ## 15-digit value SLiM gets
  p <- vivax_params(N_HOSTS = 20, N_MOSQ = 40, SIM_DAYS = 3, NU_HYPNOZOITE_BATCH = 333333.333333666,
                    NU_HYPNOZOITE_BATCH_BOOST = 3, RESURGENCE_TRIGGER_TICK = 2L)
  expect_no_error(run_sim(p, seed = 1, outdir = tempfile("mfwd_edge_")))
})

test_that("run_grid() gives each warning from its cells once, after its FAILED summary", {
  ## R keeps only the first 50 warnings, so a set that warns in every cell would crowd it out
  p <- suppressWarnings(vivax_params(N_HOSTS = 20, N_MOSQ = 40, SIM_DAYS = 3, BITING_RATE = 2))
  q <- vivax_params(N_HOSTS = 20, N_MOSQ = 40, SIM_DAYS = 3)
  q$N_MOSQ <- 40.5
  ## a set vivax_params() refuses stops the grid before any cell runs
  base <- tempfile("mfwd_grid_")
  expect_error(run_grid(list(warns = p, broken = q), seeds = 1:3, outdir_base = base),
               "run_grid() did not start: 1 of 2 parameter sets is refused:\n  broken: N_MOSQ", fixed = TRUE)
  expect_false(dir.exists(base))
  ## cells that warn, and one that fails: the FAILED summary, then each warning once
  local_mocked_bindings(run_sim = function(params, seed, ...) {
    warning("BITING_RATE = 2: a mosquito bites at most once a day")
    if (seed == 2) stop("cell error")
    list(extinct_tick = NA_integer_)
  }, .require_slim = function() "5.2")
  w <- character()
  g <- withCallingHandlers(run_grid(list(warns = p), seeds = 1:3, outdir_base = base),
                           warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  expect_identical(g$ok, c(TRUE, FALSE, TRUE))
  expect_length(w, 2)
  expect_match(w[1], "1 of 3 grid cell\\(s\\) FAILED")
  expect_match(w[2], "bites at most once a day \\[3 cells\\]$")
})

test_that("under options(warn = 2) a warning fails its cell, as R makes it an error", {
  skip_if_not(slim_available(), "SLiM not on PATH")
  p <- suppressWarnings(vivax_params(N_HOSTS = 20, N_MOSQ = 40, SIM_DAYS = 3, BITING_RATE = 2))
  op <- options(warn = 2); on.exit(options(op), add = TRUE)
  expect_error(run_grid(p, seeds = 1:2, outdir_base = tempfile("mfwd_grid2_")),
               "2 of 2 grid cell\\(s\\) FAILED.*converted from warning")
})

test_that("a caller catching run_grid()'s warnings gets the FAILED summary, and a stopped grid its error", {
  sets <- list(warns = vivax_params(SIM_ID = "warns"), broken = vivax_params(SIM_ID = "broken"))
  local_mocked_bindings(run_sim = function(params, seed, ...) {
    if (params$SIM_ID == "warns") warning("cell warning") else stop("cell error")
    list(extinct_tick = NA_integer_)
  }, .require_slim = function() "5.2")
  got <- tryCatch(run_grid(sets, seeds = 1:2, outdir_base = tempfile("mfwd_gridw_")),
                  warning = function(w) conditionMessage(w))
  expect_match(got, "^2 of 4 grid cell\\(s\\) FAILED")
  ## a value that cannot go in the result table stops the grid in its second cell: the error
  ## reaches the caller, carrying the warnings of the cells that ran, and no warning is given
  local_mocked_bindings(run_sim = function(params, seed, ...) {
    warning("cell warning")
    list(extinct_tick = if (seed == 2) sum else NA_integer_)
  }, .require_slim = function() "5.2")
  expect_no_warning(expect_error(run_grid(sets["warns"], seeds = 1:2, outdir_base = tempfile("mfwd_gridw_")),
                                 "cells that had reported gave these warnings:\n  cell warning [2 cells]",
                                 fixed = TRUE))
  ## under warn = 2 a warning fails its cell, also inside suppressWarnings()
  local_mocked_bindings(run_sim = function(params, seed, ...) {
    warning("cell warning")
    list(extinct_tick = NA_integer_)
  }, .require_slim = function() "5.2")
  op <- options(warn = 2); on.exit(options(op), add = TRUE)
  g <- suppressWarnings(run_grid(sets["warns"], seeds = 1:2, outdir_base = tempfile("mfwd_gridw_")))
  expect_false(any(g$ok))
  expect_match(g$error, "^\\(converted from warning\\) cell warning$")
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
  run <- run_sim(vivax_params(N_HOSTS = 1000, N_FOUNDER_HOSTS = 1000,
                              N_FOUNDER_LINEAGES = 1000, N_MOSQ = 1, BITING_RATE = 1e-6,
                              MEAN_BLOOD_STAGE_DURATION = 2, SIM_DAYS = 40),
                 seed = 6, outdir = tempfile("mfwd_test_clearance"))
  ev <- load_events(run)
  d <- merge(ev[event_type == "vh_inoculation_seed", .(pedigree_id, t0 = tick)],
             ev[event_type == "human_clearance", .(pedigree_id, t1 = tick)], by = "pedigree_id")
  expect_gt(nrow(d), 900)
  expect_lt(abs(mean(d$t1 - d$t0) - 2), 0.15)   # a daily probability of 1 - exp(-1/2) gives 2.54
})

test_that("no blood-stage infection is forced to clear at a fixed age", {
  skip_if_not(slim_available(), "SLiM not on PATH")
  run <- run_sim(vivax_params(N_HOSTS = 200, N_FOUNDER_HOSTS = 200,
                              N_FOUNDER_LINEAGES = 200, N_MOSQ = 1, BITING_RATE = 1e-6,
                              MEAN_BLOOD_STAGE_DURATION = 1000, SIM_DAYS = 200),
                 seed = 7, outdir = tempfile("mfwd_test_nocap"))
  ev <- load_events(run)
  d <- merge(ev[event_type == "vh_inoculation_seed", .(pedigree_id, t0 = tick)],
             ev[event_type == "human_clearance", .(pedigree_id, t1 = tick)], by = "pedigree_id")
  expect_lt(sum(d$t1 - d$t0 == 181), 5)   # a 180-day cap cleared about 83% of founders on day 181
})

test_that("latent hypnozoites die at MU_HYPNOZOITE, so a temperate batch relapses less", {
  skip_if_not(slim_available(), "SLiM not on PATH")
  ## one pulse of imports, no onward transmission, long enough for every batch to resolve
  run <- run_sim(vivax_params(N_HOSTS = 300, N_MOSQ = 10, BITING_RATE = 1e-6, SIM_DAYS = 1500,
                              N_FOUNDER_HOSTS = 0, IMPORT_TICK = 2, N_IMPORTS = 300,
                              K_LATENCY = 2L, DELTA_LATENCY = 1 / 100, CHROMOSOMES = "1"),
                 seed = 3, outdir = tempfile("mfwd_test_latency"))
  ev <- load_events(run)
  act  <- ev[event_type == "hypnozoite_activation", sum(n_hypnozoites)]
  died <- ev[event_type == "hypnozoite_clearance",  sum(n_hypnozoites)]
  a <- 1 / 170; mu <- 1 / 442; d <- 1 / 100
  ## each latent day a hypnozoite dies w.p. mu, then the batch moves on w.p. d (two stages)
  stage <- d * (1 - mu) / (1 - (1 - mu) * (1 - d))
  expect_gt(act + died, 600)
  expect_equal(act / (act + died), stage^2 * a / (a + mu), tolerance = 0.05 / 0.48)
})

test_that("clearance rows carry the is_imported of their lineage", {
  skip_if_not(slim_available(), "SLiM not on PATH")
  ## seeded only by imports, with no onward transmission, so every lineage traces to an import:
  ## batches that stay latent through the run, so that every hypnozoite death is a latent one, then
  ## batches without latency
  for (lat in list(list(K_LATENCY = 2L, DELTA_LATENCY = 1e-4), list(K_LATENCY = 0L))) {
    p <- do.call(vivax_params, c(list(N_HOSTS = 60, N_MOSQ = 10, BITING_RATE = 1e-6, SIM_DAYS = 200,
                                      N_FOUNDER_HOSTS = 0, IMPORT_TICK = 2, N_IMPORTS = 60,
                                      CHROMOSOMES = "1"), lat))
    ev <- load_events(run_sim(p, seed = 2, outdir = tempfile("mfwd_imported_")))
    expect_equal(ev[event_type == "vh_inoculation", .N], 0)
    if (lat$K_LATENCY > 0) expect_equal(ev[event_type == "hypnozoite_activation", .N], 0)
    for (type in c("hypnozoite_clearance", "human_clearance")) {
      rows <- ev[event_type == type]
      expect_gt(nrow(rows), 10, label = type)
      expect_true(all(rows$is_imported == 1), label = type)
    }
  }
})

test_that("a run that goes extinct says so and records the tick", {
  skip_if_not(slim_available(), "SLiM not on PATH")
  ## few bites and no hypnozoites: all of 60 seeds died out, by day 308 at the latest
  base <- list(N_HOSTS = 60, N_MOSQ = 60, BITING_RATE = 0.015, SIM_DAYS = 400, N_FOUNDER_HOSTS = 6,
               NU_HYPNOZOITE_BATCH = 0, CHROMOSOMES = "1")
  dies <- do.call(vivax_params, base)
  expect_message(run <- run_sim(dies, seed = 1, outdir = tempfile("mfwd_extinct_")), "went extinct on tick")
  expect_lt(run$extinct_tick, 400)
  expect_equal(max(load_state(run)$tick), run$extinct_tick)
  expect_equal(load_run(run$outdir)$extinct_tick, run$extinct_tick)
  full <- run_sim(vivax_params(N_HOSTS = 20, N_MOSQ = 40, SIM_DAYS = 5), seed = 1, outdir = tempfile("mfwd_full_"))
  expect_true(is.na(full$extinct_tick))
  ## EXTINCT_TICK is an outcome, not a parameter, so every run's parameters have the same names
  expect_null(load_params(run)$EXTINCT_TICK)
  expect_null(load_params(run$outdir)$EXTINCT_TICK)
  expect_identical(names(load_params(run)), names(load_params(full)))
  expect_message(g <- run_grid(dies, seeds = 1:2, outdir_base = tempfile("mfwd_extinct_grid_")),
                 "^2 of 2 grid cell\\(s\\) went extinct; see the `extinct_tick` column")
  expect_true(all(g$extinct_tick < 400))
  ## an import set outside the run never comes, so it does not keep a dead run going: a pulse after
  ## SIM_DAYS, a window after it, and a count whose window is left at the default days
  for (late in list(list(IMPORT_TICK = 500L, N_IMPORTS = 5L),
                    list(IMPORT_START = 450L, IMPORT_END = 600L, IMPORT_RATE_PER_TICK = 0.5),
                    list(N_IMPORTS_PER_TICK = 2L))) {
    r <- suppressMessages(run_sim(do.call(vivax_params, c(base, late)), seed = 1,
                                  outdir = tempfile("mfwd_extinct_late_")))
    expect_equal(r$extinct_tick, run$extinct_tick, label = names(late)[1])
  }
  ## the message gives the run's length as params.tsv writes it
  long <- c(utils::modifyList(base, list(SIM_DAYS = 1e5)),
            stats::setNames(as.list(rep(2e5, length(malariafwd:::.never_params))), malariafwd:::.never_params))
  expect_message(run_sim(do.call(vivax_params, long), seed = 1, outdir = tempfile("mfwd_extinct_long_")),
                 paste0("went extinct on tick ", run$extinct_tick, " of 100000: "), fixed = TRUE)
})

test_that("imports go to random hosts, and a host at MAX_MOI takes none", {
  skip_if_not(slim_available(), "SLiM not on PATH")
  ## 200 hosts, 100 of them founders, and a pulse of 100 imports on tick 3 with no transmission
  base <- list(N_HOSTS = 200, N_MOSQ = 10, BITING_RATE = 1e-6, SIM_DAYS = 5, N_FOUNDER_HOSTS = 100,
               N_FOUNDER_LINEAGES = 10, IMPORT_TICK = 3, N_IMPORTS = 100, CHROMOSOMES = "1")
  counts <- function(...) {
    run <- run_sim(do.call(vivax_params, utils::modifyList(base, list(...))), seed = 4,
                   outdir = tempfile("mfwd_import_"))
    ev <- load_events(run)
    founders <- unique(ev[event_type == "vh_inoculation_seed", human_id])
    ## the founder hosts that still carry their clone when the pulse lands: clearance runs in
    ## early() before importation, and nothing else changes a founder host's clones here
    holding <- setdiff(founders, ev[event_type == "human_clearance" & tick <= 3, human_id])
    imp <- ev[event_type == "vh_inoculation_import" & is_dormant == 0]
    c(into_founders = sum(imp$human_id %in% founders), into_holding = sum(imp$human_id %in% holding),
      refused = ev[event_type == "vh_inoculation_import_no_take", .N])
  }
  open <- counts()
  expect_gt(open[["into_founders"]], 30)
  expect_lt(open[["into_founders"]], 70)
  capped <- counts(MAX_MOI = 1L)
  expect_equal(capped[["into_holding"]], 0)
  expect_gt(capped[["refused"]], 30)
  ## a pulse larger than the population brings one import to each host, and warns
  expect_warning(over <- vivax_params(N_HOSTS = 30, N_MOSQ = 10, BITING_RATE = 1e-6, SIM_DAYS = 4,
                                      IMPORT_TICK = 2, N_IMPORTS = 50, CHROMOSOMES = "1"),
                 "at most one import a day")
  expect_warning(ev <- load_events(run_sim(over, seed = 1, outdir = tempfile("mfwd_import_"))),
                 "at most one import a day")
  imp <- ev[event_type == "vh_inoculation_import" & is_dormant == 0]
  expect_equal(nrow(imp), 30)
  expect_equal(length(unique(imp$human_id)), 30)
})

test_that("a bite on a host at MAX_MOI is refused whether or not the host is infectious", {
  skip_if_not(slim_available(), "SLiM not on PATH")
  ## every course fails into subpatency, and the suppressed clone stays for the rest of the run;
  ## at SUBPATENT_INFECTIVITY = 0 such a host is not infectious, but its clone still counts
  p <- vivax_params(N_HOSTS = 60, N_MOSQ = 200, BITING_RATE = 0.2, SIM_DAYS = 200, MAX_MOI = 1L,
                    SUBPATENT_INFECTIVITY = 0, BLOODCURE_COVERAGE = 1, P_BLOODCURE_SURVIVAL = 1,
                    P_RECRUD_GIVEN_FAIL = 1, RECRUD_DELAY_DIST = "fixed", MEAN_RECRUD_DELAY = 1e5,
                    NU_HYPNOZOITE_BATCH = 0, CHROMOSOMES = "1")
  ev <- load_events(run_sim(p, seed = 1, outdir = tempfile("mfwd_moi_subpatent_")))
  suppressed <- ev[event_type == "bloodcure_subpatent", .(since = min(tick)), by = human_id]
  bites <- ev[event_type %in% c("vh_inoculation", "vh_inoculation_no_take")]
  expect_gt(nrow(suppressed), 20)
  expect_gt(nrow(bites), 20)
  expect_equal(merge(bites, suppressed, by = "human_id")[tick > since, .N], 0)
})

test_that("a pulse of imports on the first day keeps the founders", {
  skip_if_not(slim_available(), "SLiM not on PATH")
  run <- run_sim(vivax_params(N_HOSTS = 40, N_MOSQ = 40, SIM_DAYS = 5, IMPORT_TICK = 1, N_IMPORTS = 5,
                              CHROMOSOMES = "1"), seed = 1, outdir = tempfile("mfwd_pulse1_"))
  ev <- load_events(run)
  expect_equal(ev[event_type == "vh_inoculation_seed", length(unique(human_id))], 12)
  expect_gt(ev[event_type == "vh_inoculation_import" & tick == 1, .N], 0)
})

test_that("into_patent_lineage marks a relapse into a patent clone of its own lineage", {
  skip_if_not(slim_available(), "SLiM not on PATH")
  ## treatment failures and an MDA round, so clones are also cleared, suppressed and recrudesce
  run <- run_sim(vivax_params(N_HOSTS = 60, N_MOSQ = 150, SIM_DAYS = 400, CHROMOSOMES = "1",
                              BLOODCURE_COVERAGE = 0.5, P_BLOODCURE_SURVIVAL = 0.3,
                              P_RECRUD_GIVEN_FAIL = 0.5, MDA_TICK = 200, MDA_COVERAGE = 0.5),
                 seed = 3, outdir = tempfile("mfwd_masked_"))
  ev <- load_events(run)
  et <- event_types()
  ## the dictionary names every type the run logged
  expect_true(all(unique(ev$event_type) %in% et$event_type))
  ## each blood-stage clone is patent from its onset row to its first row that ends an episode
  on <- ev[event_type %in% et$event_type[et$is_onset] & is_dormant == 0,
           .(pedigree_id, human_id, lineage_id, start = tick)]
  ends <- ev[event_type %in% et$event_type[et$ends_episode], .(end = min(tick)), by = pedigree_id]
  spans <- merge(on, ends, by = "pedigree_id", all.x = TRUE)
  act <- ev[event_type == "hypnozoite_activation"]
  ## an activation runs in reproduction(): it sees clones born on earlier days, and a clone that
  ## ends on its own day ends after it
  seen <- vapply(seq_len(nrow(act)), function(i) {
    s <- spans[human_id == act$human_id[i] & lineage_id == act$lineage_id[i]]
    as.integer(any(s$start < act$tick[i] & (is.na(s$end) | s$end >= act$tick[i])))
  }, integer(1))
  expect_gt(nrow(act), 50)
  expect_gt(sum(ev$event_type == "recrudescence"), 0)
  expect_gt(sum(ev$event_type == "mda_bloodstage_clear"), 0)
  expect_gt(sum(seen), 0)
  expect_lt(sum(seen), nrow(act))
  expect_equal(act$into_patent_lineage, seen)
  expect_true(all(ev[event_type != "hypnozoite_activation", into_patent_lineage] == -1))
})

test_that("run_sim() says when it stops a run at its timeout", {
  skip_on_os("windows")
  ## a stand-in for SLiM that reports its version and then hangs
  fake <- tempfile("fakeslim")
  writeLines(c("#!/bin/sh", 'if [ "$1" = "-v" ]; then echo "SLiM version 5.2"; exit 0; fi', "exec sleep 30"), fake)
  Sys.chmod(fake, "755")
  old <- Sys.getenv("MALARIAFWD_SLIM", NA)
  Sys.setenv(MALARIAFWD_SLIM = fake)
  on.exit(if (is.na(old)) Sys.unsetenv("MALARIAFWD_SLIM") else Sys.setenv(MALARIAFWD_SLIM = old), add = TRUE)
  expect_error(run_sim(vivax_params(), outdir = tempfile("mfwd_timeout_"), timeout = 1),
               "SLiM run timed out after 1 s and was stopped")
})

test_that("a founder-free run whose imports never take ends cleanly as never started", {
  skip_if_not(slim_available(), "SLiM not on PATH")
  ## an import rate so low that no import arrives in the window
  p <- vivax_params(N_HOSTS = 20, N_MOSQ = 20, SIM_DAYS = 60, N_FOUNDER_HOSTS = 0, IMPORT_START = 1,
                    IMPORT_END = 20, IMPORT_RATE_PER_TICK = 1e-6, CHROMOSOMES = "1")
  expect_message(run <- run_sim(p, seed = 1, outdir = tempfile("mfwd_nostart_")),
                 "The run never started and ended on tick 20 of 60: no import established a parasite")
  expect_true(run$never_started)
  expect_true(is.na(run$extinct_tick))
  expect_null(load_params(run)$NEVER_STARTED)
  expect_true(file.exists(trees_path(run)))
  expect_equal(max(load_state(run)$tick), 20)
  expect_output(print(run), "never started: no import established a parasite")
  expect_message(g <- run_grid(p, seeds = 1:2, outdir_base = tempfile("mfwd_nostart_grid_")),
                 "2 of 2 grid cell\\(s\\) never started")
  expect_identical(g$never_started, c(TRUE, TRUE))
  expect_true(all(is.na(g$extinct_tick)))
})

test_that("a founder-free run is marked never started whatever its import schedule, once it ends", {
  skip_if_not(slim_available(), "SLiM not on PATH")
  ## the tick a run ended on, from state.tsv, and its two outcomes
  outcome <- function(...) {
    p <- vivax_params(N_HOSTS = 20, N_MOSQ = 20, N_FOUNDER_HOSTS = 0, CHROMOSOMES = "1", ...)
    run <- suppressMessages(run_sim(p, seed = 1, outdir = tempfile("mfwd_nostart_")))
    list(end = max(load_state(run)$tick), never_started = run$never_started,
         extinct_tick = run$extinct_tick, events = load_events(run))
  }
  ## imports that can take in neither compartment
  no_take <- list(P_PRIMARY = 0, NU_HYPNOZOITE_BATCH = 0)
  ## a pulse on day 1 waits for tick 10; one on day 15 ends on the day after
  o <- do.call(outcome, c(no_take, SIM_DAYS = 60, IMPORT_TICK = 1, N_IMPORTS = 3))
  expect_identical(o[1:3], list(end = 10L, never_started = TRUE, extinct_tick = NA_integer_))
  expect_identical(o$events$event_type, rep("vh_inoculation_import_no_take", 3))
  o <- do.call(outcome, c(no_take, SIM_DAYS = 60, IMPORT_TICK = 15, N_IMPORTS = 3))
  expect_identical(o[1:2], list(end = 16L, never_started = TRUE))
  ## an import path that lasts to the run's last day, or a run shorter than 10 days, is marked
  ## on its last day
  o <- outcome(SIM_DAYS = 40, IMPORT_START = 1, IMPORT_RATE_PER_TICK = 1e-6)
  expect_identical(o[1:3], list(end = 40L, never_started = TRUE, extinct_tick = NA_integer_))
  o <- do.call(outcome, c(no_take, SIM_DAYS = 30, IMPORT_TICK = 30, N_IMPORTS = 1))
  expect_identical(o[1:2], list(end = 30L, never_started = TRUE))
  o <- do.call(outcome, c(no_take, SIM_DAYS = 8, IMPORT_TICK = 1, N_IMPORTS = 2))
  expect_identical(o[1:2], list(end = 8L, never_started = TRUE))
  ## an import that took, as hypnozoites only, and died before activating: the run started, and
  ## went extinct (at a batch mean of 1e6 an empty batch, which would not take, is all but
  ## impossible)
  o <- outcome(SIM_DAYS = 60, IMPORT_TICK = 1, N_IMPORTS = 1, P_PRIMARY = 0,
               NU_HYPNOZOITE_BATCH = 1e6, MU_HYPNOZOITE = 1, ALPHA_ACTIVATE = 0)
  expect_true("vh_inoculation_import" %in% o$events$event_type)
  expect_identical(o[1:3], list(end = 2L, never_started = FALSE, extinct_tick = 2L))
})

## A set whose floats need all 15 digits, a probability R stores as SLiM is given it, whole numbers
## beyond the integer range, and a UTF-8 SIM_ID with an apostrophe.
record_set <- function() vivax_params(N_HOSTS = 30, N_MOSQ = 60, SIM_DAYS = 60, ALPHA_ACTIVATE = 1/170,
  BITING_RATE = 0.1 + 0.02, SEASONAL_AMP = 1/3, SEASONAL_PHASE = 100/7, P_TRANS_VH = 1 + 2e-16,
  LAMBDA_SPOROZOITES = 2/3, LAMBDA_OOCYSTS = exp(1), K_OOCYSTS = 0.003, RECOMB_RATE = 1e-7/3,
  INTERVENTION_START = 3e9, INTERVENTION_END = 3e9, MDA_TICK = 30, MDA_COVERAGE = 0.123456789012345678,
  SIM_ID = intToUtf8(c(0x43, 0xf4, 0x74, 0x65, 0x27, 0x73)), GENOME_MODE = "core")

test_that("load_params() gives the set a run was made from, and run_sim() of it repeats the run", {
  skip_if_no_slim()
  p <- record_set()
  r1 <- suppressMessages(run_sim(p, seed = 7, outdir = tempfile("mfwd_rec1_")))
  expect_identical(load_params(r1), p)
  expect_identical(load_params(r1$outdir), p)
  r2 <- suppressMessages(run_sim(load_params(r1), seed = r1$seed, outdir = tempfile("mfwd_rec2_")))
  files <- c("events.tsv", "state.tsv", "params.tsv")
  expect_identical(unname(tools::md5sum(file.path(r2$outdir, files))),
                   unname(tools::md5sum(file.path(r1$outdir, files))))
})

test_that("load_run() carries the run's record, typed", {
  skip_if_no_slim()
  r <- suppressMessages(run_sim(record_set(), seed = 7, outdir = tempfile("mfwd_rec3_")))
  expect_identical(r$seed, 7L)
  expect_identical(r$slim_version, slim_version())
  expect_identical(r$malariafwd_version, as.character(utils::packageVersion("malariafwd")))
  expect_identical(r$n_chromosomes, 14L)
  expect_identical(r$genome_length, as.integer(sum(chrom_map(r)$length)))
  expect_identical(r$recomb_within, 3.33333333333333e-08)
  expect_identical(r$extinct_tick, NA_integer_)
  expect_false(r$never_started)
  expect_identical(load_run(r), r)
  expect_output(print(r), "60 days, 30 hosts, 60 mosquitoes")
  expect_output(print(r), "SLiM [0-9.]+, malariafwd")
  ## the test helper writes params.tsv as the model does
  expect_false(any(startsWith(readLines(file.path(r$outdir, "params.tsv")), "N_FOUNDERS\t")))
  rec <- malariafwd:::.read_params_tsv(file.path(r$outdir, "params.tsv"))$record
  expect_setequal(params_tsv_lines(load_params(r), rec), readLines(file.path(r$outdir, "params.tsv"), encoding = "UTF-8"))
})
