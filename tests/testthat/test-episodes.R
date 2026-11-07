## Episodes and relapses: the masked relapses, into a patent clone of their own lineage.

test_that("classify_cases() marks the masked relapse, and leaves it out when asked", {
  fx <- episodes_fixture()
  cs <- classify_cases(fx)
  expect_identical(nrow(cs), 12L)
  expect_identical(cs$pedigree_id[cs$into_patent_lineage], 104L)
  expect_false(anyNA(cs$into_patent_lineage))
  cm <- classify_cases(fx, exclude_masked = TRUE)
  expect_identical(nrow(cm), 11L)
  expect_false(104L %in% cm$pedigree_id)
})

test_that("relapse_fraction() leaves the masked relapse out of both readings when asked", {
  rf <- relapse_fraction(episodes_fixture(), exclude_masked = TRUE)
  expect_identical(rf$n_relapse, 3L)
  expect_equal(rf$relapse_fraction, 3 / 10)
  expect_equal(rf$relapse_frac_episodes, 3 / 11)
  expect_identical(relapse_fraction(episodes_fixture())$n_relapse, 4L)
})

test_that("force_of_infection() leaves the masked relapse out when asked", {
  r <- episodes_fixture_run()
  all <- force_of_infection(r); unm <- force_of_infection(r, exclude_masked = TRUE)
  rel <- function(fo) sum(fo[grepl("^relapse", origin)]$n)
  expect_identical(c(rel(all), rel(unm)), c(4L, 3L))
  expect_identical(all[!grepl("^relapse", origin)], unm[!grepl("^relapse", origin)])
})

test_that("exclude_masked must be TRUE or FALSE", {
  r <- episodes_fixture_run(); fx <- episodes_fixture()
  for (v in list(NA, "yes", c(TRUE, FALSE))) {
    expect_error(classify_cases(fx, exclude_masked = v), "`exclude_masked` must be TRUE or FALSE")
    expect_error(relapse_fraction(fx, exclude_masked = v), "`exclude_masked` must be TRUE or FALSE")
    expect_error(force_of_infection(r, exclude_masked = v), "`exclude_masked` must be TRUE or FALSE")
  }
})

test_that("a run without the flag gives NA on its relapses, and cannot leave them out", {
  fx <- episodes_fixture()
  for (ev in list(data.table::copy(fx)[event_type == "hypnozoite_activation", into_patent_lineage := NA_integer_],
                  data.table::copy(fx)[, into_patent_lineage := NULL])) {
    cs <- classify_cases(ev)
    expect_true(all(is.na(cs[origin == "relapse"]$into_patent_lineage)))
    expect_false(anyNA(cs[origin != "relapse"]$into_patent_lineage))
    expect_error(classify_cases(ev, exclude_masked = TRUE), "4 of these events' activation rows lack it")
  }
})

test_that("on a run, exclude_masked removes exactly the flagged activations", {
  skip_if_no_slim()
  r <- suppressMessages(run_sim(vivax_params(N_HOSTS = 60, N_MOSQ = 150, SIM_DAYS = 400, CHROMOSOMES = "1",
                                             BLOODCURE_COVERAGE = 0.5, P_BLOODCURE_SURVIVAL = 0.3,
                                             P_RECRUD_GIVEN_FAIL = 0.5, MDA_TICK = 200, MDA_COVERAGE = 0.5),
                                seed = 3, outdir = tempfile("mfwd_masked_")))
  ev <- load_events(r)
  flagged <- ev[event_type == "hypnozoite_activation" & into_patent_lineage == 1L, .N]
  expect_gt(flagged, 0L)
  expect_identical(nrow(classify_cases(ev)) - nrow(classify_cases(ev, exclude_masked = TRUE)), flagged)
})
