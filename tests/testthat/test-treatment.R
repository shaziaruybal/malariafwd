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

  ## radical cure is NESTED inside blood-stage treatment, so it needs both coverages
  rad <- run_counts(shared(BLOODCURE_COVERAGE = 1, RADCURE_COVERAGE = 1), "rad")
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

test_that("the case-management window bounds WHEN treatment happens", {
  skip_if_no_slim()
  ## Deliberately an assertion about TIMING, not about counts. Since treatment became
  ## host-level it suppresses transmission hard enough that an always-open window can
  ## produce FEWER courses than a late one - it prevents the cases it would have treated.
  out <- file.path(tempdir(), "mfw-tx-window")
  unlink(out, recursive = TRUE)
  r <- run_sim(shared(BLOODCURE_COVERAGE = 1, TREATMENT_START = 90L, TREATMENT_END = 120L),
               seed = 2, outdir = out)
  ev <- data.table::as.data.table(load_events(r))
  given <- ev[event_type == "bloodcure_given", tick]
  expect_gt(length(given), 0)
  ## Nothing before the window opens: no case can be DECIDED earlier.
  expect_true(all(given >= 90))
  ## Courses decided just before it shuts are delivered TREATMENT_DELAY days later, so a
  ## short tail past TREATMENT_END is correct - the window gates prescription, not delivery.
  expect_true(all(given < 120 + 5 * 7))
})

test_that("a coverage against a closed treatment window is refused", {
  ## The failure mode D19 fixed was silence, so re-creating it in the NEW parameter has to
  ## be an error rather than a quietly empty result.
  expect_error(vivax_params(BLOODCURE_COVERAGE = 1, TREATMENT_START = 50L, TREATMENT_END = 50L),
               "case-management window is closed")
  ## suppressWarnings: this also trips the nested-radical-cure warning, which is correct
  ## but not what this expectation is about
  suppressWarnings(
    expect_error(vivax_params(RADCURE_COVERAGE = 1, SIM_DAYS = 100L, TREATMENT_START = 500L),
                 "case-management window is closed"))
  ## the defaults are "always on", so a coverage on its own is fine
  expect_no_error(vivax_params(BLOODCURE_COVERAGE = 1))
})
