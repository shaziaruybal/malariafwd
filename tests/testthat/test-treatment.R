## Case management. These need SLiM (skip_if_no_slim() lives in helper-slim.R).

shared <- function(...) vivax_params(SIM_DAYS = 120L, N_HOSTS = 40L, N_MOSQ = 80L,
                                     BITING_RATE = 0.12, ...)

run_counts <- function(p, tag) {
  out <- file.path(tempdir(), paste0("mfw-tx-", tag))
  unlink(out, recursive = TRUE)
  r <- run_sim(p, seed = 2, outdir = out)
  ev <- data.table::as.data.table(load_events(r))
  list(bloodcure = ev[event_type == "case_bloodstage_clear", .N],
       radcure   = ev[event_type == "radcure_hypnozoite_clear", .N])
}

test_that("treatment is driven by coverage alone, not by the vector-control window", {
  skip_if_no_slim()
  ## D19. Every case-based treatment used to be gated on INTERVENTION_START/END - the same
  ## window that scales biting and importation - so setting a coverage without also opening
  ## that window produced zero treatments and no warning. That is how the issue was found.
  expect_equal(run_counts(shared(), "none")$bloodcure, 0)

  blood <- run_counts(shared(BLOODCURE_COVERAGE = 1), "blood")
  expect_gt(blood$bloodcure, 0)

  rad <- run_counts(shared(RADCURE_COVERAGE = 1), "rad")
  expect_gt(rad$radcure, 0)
})

test_that("the vector-control window does not switch treatment on", {
  skip_if_no_slim()
  ## ...and conversely: opening the intervention window must change transmission without
  ## producing a single treatment, or the two are still entangled.
  vc <- run_counts(shared(INTERVENTION_START = 1L, INTERVENTION_END = 120L,
                          INTERVENTION_EFFICACY = 0.9), "vc")
  expect_equal(vc$bloodcure, 0)
  expect_equal(vc$radcure, 0)
})

test_that("the case-management window bounds when treatment happens", {
  skip_if_no_slim()
  full <- run_counts(shared(BLOODCURE_COVERAGE = 1), "win-full")
  late <- run_counts(shared(BLOODCURE_COVERAGE = 1,
                            TREATMENT_START = 90L, TREATMENT_END = 120L), "win-late")
  expect_gt(full$bloodcure, late$bloodcure)
  expect_gt(late$bloodcure, 0)
})

test_that("a coverage against a closed treatment window is refused", {
  ## The failure mode D19 fixed was silence, so re-creating it in the NEW parameter has to
  ## be an error rather than a quietly empty result.
  expect_error(vivax_params(BLOODCURE_COVERAGE = 1, TREATMENT_START = 50L, TREATMENT_END = 50L),
               "case-management window is closed")
  expect_error(vivax_params(RADCURE_COVERAGE = 1, SIM_DAYS = 100L, TREATMENT_START = 500L),
               "case-management window is closed")
  ## the defaults are "always on", so a coverage on its own is fine
  expect_no_error(vivax_params(BLOODCURE_COVERAGE = 1))
})
