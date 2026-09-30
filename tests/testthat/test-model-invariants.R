## Model invariants that guard against regressions of the WS-G bugs. These need SLiM,
## so they skip when it is unavailable (e.g. on CRAN / CI without the engine).

skip_if_no_slim <- function() skip_if_not(slim_available(), "SLiM not on PATH")

test_that("hypnozoite batches never draw more than they hold, and ids are unique", {
  skip_if_no_slim()
  out <- file.path(tempdir(), "mfw-invariants")
  unlink(out, recursive = TRUE)
  p <- vivax_params(N_HOSTS = 60, N_MOSQ = 120, BITING_RATE = 0.12, SIM_DAYS = 300,
                    N_FOUNDER_LINEAGES = 8)
  run <- run_sim(p, seed = 1, outdir = out)
  ev <- as.data.table(load_events(run))

  ## (WS-G bug 1) a batch's activations + deaths must never exceed the size it was
  ## established with. The old PENDING_HYP_DEATHS queue was scoped local, so retired
  ## batches stayed alive and kept activating (64% of batches over-drew).
  est <- ev[grepl("^vh_inoculation", event_type) & is_dormant == 1,
            .(batch = pedigree_id, held = n_hypnozoites)]
  act <- ev[event_type == "hypnozoite_activation", .(batch = parent1_pedigree_id, n = n_hypnozoites)]
  die <- ev[grepl("hypnozoite_clear", event_type), .(batch = pedigree_id, n = n_hypnozoites)]
  drawn <- rbind(act, die)[, .(drew = sum(n)), by = batch]
  chk <- merge(est, drawn, by = "batch", all.x = TRUE)[is.na(drew), drew := 0]
  expect_equal(sum(chk$drew > chk$held), 0)

  ## (WS-G bugs 2 & 3) event ids unique; inoculation ids do not restart each tick
  expect_equal(uniqueN(ev$event_id), nrow(ev))
  inoc <- ev[inoculation_id > 0, inoculation_id]
  expect_gt(max(inoc), 0.5 * uniqueN(inoc))
})

test_that("classify_cases labels imported relapses", {
  skip_if_no_slim()
  out <- file.path(tempdir(), "mfw-imported")
  unlink(out, recursive = TRUE)
  p <- vivax_params(N_HOSTS = 60, N_MOSQ = 120, BITING_RATE = 0.06, SIM_DAYS = 400,
                    IMPORT_START = 1L, IMPORT_END = 99999L, IMPORT_RATE_PER_TICK = 0.15)
  run <- run_sim(p, seed = 1, outdir = out)
  cs <- classify_cases(load_events(run))
  expect_true("imported" %in% names(cs))
  ## a relapse of an imported genotype must be labelled imported (was hard-coded FALSE)
  expect_gt(nrow(cs[origin == "relapse" & imported == TRUE]), 0)
})

test_that("deprecated N_OOCYSTS_PER_BITE is rejected", {
  expect_error(vivax_params(N_OOCYSTS_PER_BITE = 2), "deprecated")
})

test_that("zero-truncated count helpers round-trip", {
  ## lambda is the PRE-truncation rate: the realised mean is higher
  expect_gt(uptake_moments(2)[["mean"]], 2)
  expect_equal(unname(uptake_moments(2)[["mean"]]), 2 / (1 - exp(-2)), tolerance = 1e-8)
  ## NB adds dispersion at the same lambda
  expect_gt(uptake_moments(2, k = 1)[["var"]], uptake_moments(2, k = 0)[["var"]])
  ## inverse recovers the targets
  p <- uptake_calibrate(2.5)
  expect_equal(unname(p$achieved[["mean"]]), 2.5, tolerance = 1e-6)
  expect_equal(p$k, 0)
  q <- uptake_calibrate(2.5, var = 6)
  expect_equal(unname(q$achieved[["mean"]]), 2.5, tolerance = 1e-3)
  expect_equal(unname(q$achieved[["var"]]),  6.0, tolerance = 1e-3)
  expect_gt(q$k, 0)
  expect_error(uptake_calibrate(0.8), "mean > 1")
})

test_that("ARG summaries behave on the forest genealogy", {
  skip_if_no_slim()
  skip_if_not(requireNamespace("reticulate", quietly = TRUE), "reticulate not available")
  out <- file.path(tempdir(), "mfw-arg")
  unlink(out, recursive = TRUE)
  run <- run_sim(vivax_params(N_HOSTS = 60, N_MOSQ = 120, BITING_RATE = 0.12,
                              SIM_DAYS = 300, N_FOUNDER_LINEAGES = 8), seed = 1, outdir = out)
  ev <- as.data.table(load_events(run))
  s <- arg_summary(run)
  expect_gt(s$num_trees, 1)          # recombination happened
  expect_gte(s$roots_min, 1)
  ## the genealogy is a FOREST: independent founder lineages never coalesce
  expect_gt(s$roots_median, 1)

  ## arg_tmrca is the relatedness measure that works on a forest
  rel <- ev[event_type == "hypnozoite_activation" & parent1_pedigree_id > 0][1:5]
  if (nrow(rel) >= 1 && !any(is.na(rel$pedigree_id))) {
    tm <- arg_tmrca(run, rel$pedigree_id, rel$parent1_pedigree_id)
    expect_true(all(tm$coalesced > 0.99, na.rm = TRUE))   # a relapse is clonal to its origin
  }
  fnd <- unique(ev[event_type == "vh_inoculation_seed"], by = "lineage_id")$pedigree_id
  if (length(fnd) >= 4) {
    tf <- arg_tmrca(run, fnd[1:2], fnd[3:4])
    expect_true(all(tf$coalesced < 0.5, na.rm = TRUE))    # distinct founders are independent roots
  }
})
