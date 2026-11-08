## Relapse metrics: relapse_fraction()'s two readings and relapses per person-year.

test_that("relapse_fraction() gives relapses over every case and over the acquisitions", {
  rf <- relapse_fraction(episodes_fixture())
  expect_identical(names(rf), c("n_relapse", "n_new_infection", "n_recrudescence",
                                "relapse_frac_episodes", "relapse_frac_acquisitions"))
  expect_identical(c(rf$n_relapse, rf$n_new_infection, rf$n_recrudescence), c(4L, 7L, 1L))
  expect_identical(c(rf$relapse_frac_episodes, rf$relapse_frac_acquisitions), c(4 / 12, 4 / 11))
})

test_that("relapse_fraction() leaves the masked relapse out of both readings when asked", {
  rf <- relapse_fraction(episodes_fixture(), exclude_masked = TRUE)
  expect_identical(rf$n_relapse, 3L)
  expect_identical(c(rf$relapse_frac_episodes, rf$relapse_frac_acquisitions), c(3 / 11, 3 / 10))
})

test_that("a run with no case gives zero counts and no fraction, without a warning", {
  expect_no_warning(rf <- relapse_fraction(fake_run()))
  expect_identical(c(rf$n_relapse, rf$n_new_infection, rf$n_recrudescence), c(0L, 0L, 0L))
  expect_identical(c(rf$relapse_frac_episodes, rf$relapse_frac_acquisitions), c(NA_real_, NA_real_))
})

test_that("relapse_fraction() reads a run as it reads its events", {
  expect_identical(relapse_fraction(episodes_fixture_run()), relapse_fraction(episodes_fixture()))
  expect_identical(relapse_fraction(episodes_fixture_run(), exclude_masked = TRUE),
                   relapse_fraction(episodes_fixture(), exclude_masked = TRUE))
})

test_that("relapses_per_person_year() counts every relapse over the person-time of hosts ever infected", {
  rp <- relapses_per_person_year(episodes_fixture_run())
  expect_identical(names(rp), c("n_relapse", "n_hosts", "person_years", "relapses_per_py"))
  ## host 5 from its batch-only bite on day 7, host 2 from its bite on day 25 (not its seeding),
  ## every host to day 100, not the last event's day 80
  expect_identical(c(rp$n_relapse, rp$n_hosts), c(4L, 5L))
  expect_equal(rp$person_years, 445 / 365)
  expect_equal(rp$relapses_per_py, 4 / (445 / 365))
  m <- malariafwd:::.relapses_per_person_year(episodes_fixture(), 100L, exclude_masked = TRUE)
  expect_identical(m$n_relapse, 3L)
  expect_equal(m$person_years, 445 / 365)
})

test_that("a window restricts both the relapses and the person-time", {
  r <- episodes_fixture_run()
  w <- function(...) { x <- relapses_per_person_year(r, ...); c(x$n_relapse, x$person_years * 365) }
  expect_equal(w(from = 40), c(2, 305))
  expect_equal(w(to = 30), c(2, 95))
  expect_equal(w(from = 100, to = 100), c(0, 5))
})

test_that("relapses_per_person_year() refuses a window outside the run, and an events table", {
  r <- episodes_fixture_run()
  expect_error(relapses_per_person_year(r, to = 101), "`to` = 101 is after the run's last day, SIM_DAYS = 100")
  expect_error(relapses_per_person_year(r, from = 50, to = 40), "`from` = 50 is after `to` = 40")
  expect_error(relapses_per_person_year(r, from = -1), "`from` must be a single whole number")
  expect_error(relapses_per_person_year(episodes_fixture()), "needs the run for its length")
})

test_that("a run with no infection gives no person-time, without a warning", {
  r <- episodes_fixture_run(events = episodes_fixture()[0])
  expect_no_warning(rp <- relapses_per_person_year(r))
  expect_identical(c(rp$n_relapse, rp$n_hosts), c(0L, 0L))
  expect_identical(c(rp$person_years, rp$relapses_per_py), c(0, NA_real_))
})

test_that("on a run every relapse is counted, and person-time runs to SIM_DAYS after an extinction", {
  skip_if_no_slim()
  r <- suppressMessages(run_sim(vivax_params(N_HOSTS = 100, N_MOSQ = 200, BITING_RATE = 0.11, SIM_DAYS = 548,
                                             P_PRIMARY = 0.5, CHROMOSOMES = "1"),
                                seed = 1, outdir = tempfile("mfwd_rpy_")))
  ev <- load_events(r)
  expect_identical(relapses_per_person_year(r)$n_relapse, ev[event_type == "hypnozoite_activation", .N])
  w <- w_run(2); evw <- load_events(w)
  expect_lt(max(load_state(w)$tick), 400L)
  t0 <- evw[event_type %in% c("vh_inoculation", "vh_inoculation_import"), .(t0 = min(tick)), by = human_id]$t0
  expect_equal(relapses_per_person_year(w)$person_years, sum(400 - t0) / 365)
})
