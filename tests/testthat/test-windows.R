## Every run-length window is SIM_DAYS long, so the days after an extinction count as days with no
## parasite.

test_that("summarise_run() counts the days after an extinction as days with no infected host", {
  skip_if_no_slim()
  for (s in 1:6) {
    r <- w_run(s)
    expect_lt(r$extinct_tick, 200L)
    m <- summarise_run(r)
    expect_identical(c(m$prevalence_frac, m$hypnozoites_per_host), c(0, 0), label = paste("seed", s))
    expect_true(is.na(m$coi_mean) && is.na(m$polyclonal_frac), label = paste("seed", s))
  }
})

test_that("force_of_infection() gives each window its days within SIM_DAYS", {
  skip_if_no_slim()
  hit <- FALSE
  for (s in 1:6) {
    r <- w_run(s)
    fo <- force_of_infection(r)[win > 0]
    expect_identical(fo$days, pmin(30L, 400L - fo$win + 1L), label = paste("seed", s))
    hit <- hit || any(fo$win == (r$extinct_tick %/% 30) * 30 & fo$n > 0)
  }
  expect_true(hit)
})

test_that("new_infection_diversity() takes its burn-in from SIM_DAYS", {
  skip_if_no_slim()
  d <- new_infection_diversity(w_run(2))
  expect_identical(d$n, 0L)
  expect_true(is.na(d$mean_ibd))
})

test_that("plot_diagnostics() runs to SIM_DAYS, at zero after an extinction", {
  skip_if_no_slim()
  skip_if_not_installed("ggplot2"); skip_if_not_installed("patchwork")
  r <- w_run(2)
  p <- suppressWarnings(plot_diagnostics(r))
  prev <- p[[1]]$data
  expect_identical(range(prev$w), c(0, 360))
  expect_true(all(prev$y[prev$w > (r$extinct_tick %/% 30) * 30] == 0))
})

test_that("a run that reaches SIM_DAYS has its prevalence over the second half of its days", {
  skip_if_no_slim()
  r <- suppressMessages(run_sim(vivax_params(N_HOSTS = 30, N_MOSQ = 60, SIM_DAYS = 40), seed = 4,
                                outdir = tempfile("mfwd_full_")))
  st <- load_state(r)
  m <- summarise_run(r)
  expect_equal(m$prevalence_frac, mean(st[tick >= 21]$moi_active > 0))
})

test_that(".sim_days() reads SIM_DAYS, and refuses an events table", {
  expect_identical(malariafwd:::.sim_days(fake_run(params = list(SIM_DAYS = 100))), 100L)
  r <- fake_run(drop = "SIM_DAYS", extra = "SIM_DAYS\t400.0")
  expect_identical(malariafwd:::.sim_days(r), 400L)
  expect_error(malariafwd:::.sim_days(load_events(r), "summarise_run"), "summarise_run\\(\\) needs the run for its length")
})

## a case row of events.tsv: a local inoculation of host `h` on day `t`
case_row <- function(t, h, ped) paste(t, paste0(t, "_", h), "vh_inoculation", 1, 3, h, ped, ped, ped, 0, "0.0", "human", 0, 0,
                                     ped, -1, -1, 0, 2, 1, -1, sep = "\t")

test_that("force_of_infection() has a row for every window and origin, each over the days it covers", {
  r <- fake_run(params = list(N_HOSTS = 10, SIM_DAYS = 200, N_FOUNDER_HOSTS = 1),
                ev_rows = c(case_row(5, 1, 11), case_row(6, 2, 12), case_row(95, 3, 13)))
  fo <- force_of_infection(r)
  expect_identical(nrow(fo), 35L)
  expect_identical(vapply(fo, function(x) class(x)[1], ""),
                   c(win = "integer", origin = "factor", n = "integer", days = "integer", rate = "numeric",
                     is_acquisition = "logical"))
  expect_identical(unique(fo$win), seq(0L, 180L, by = 30L))
  expect_identical(fo[, sum(n), by = win]$V1, c(2L, 0L, 0L, 1L, 0L, 0L, 0L))
  expect_identical(fo[origin == "imported", days], c(29L, 30L, 30L, 30L, 30L, 30L, 21L))
  expect_identical(fo[, sum(days), by = origin]$V1, rep(200L, 5))
  expect_equal(fo$rate, fo$n / 10 / (fo$days / 365))
  ## a run with no case: every window, at zero
  z <- force_of_infection(fake_run(params = list(N_HOSTS = 10, SIM_DAYS = 200, N_FOUNDER_HOSTS = 1)))
  expect_identical(nrow(z), 35L)
  expect_true(all(z$n == 0L))
  expect_identical(levels(z$origin), levels(fo$origin))
  for (b in list(0, 2.5, c(30, 60), NA, "30")) expect_error(force_of_infection(r, bin = b), "`bin` must be a whole number")
})

test_that("plot_force_of_infection() draws both series for a run with no acquisition", {
  skip_if_not_installed("ggplot2")
  z <- force_of_infection(fake_run(params = list(N_HOSTS = 10, SIM_DAYS = 200, N_FOUNDER_HOSTS = 1)))
  expect_no_warning(b <- ggplot2::ggplot_build(plot_force_of_infection(z, smooth = NA)))
  expect_identical(nrow(b$data[[1]]), 14L)
})

test_that("an extinct run's windows after its extinction are zeros", {
  skip_if_no_slim()
  r <- w_run(2)
  fo <- force_of_infection(r)
  expect_identical(nrow(fo), (400L %/% 30L + 1L) * 5L)
  expect_true(all(fo[win > r$extinct_tick]$n == 0L))
})

test_that("new_infection_diversity() starts on summarise_run()'s first day", {
  r <- episodes_fixture_run(state = fixture_state())
  ## burn_frac 0.45 of 100 days: day 45, with acquisition 404, is the last burn-in day
  expect_identical(summarise_run(r, burn_frac = 0.45)$from_tick, 46L)
  expect_identical(new_infection_diversity(r, burn_frac = 0.45)$n, 1L)
  expect_identical(new_infection_diversity(r, burn_frac = 0.44)$n, 2L)
})

test_that("force_of_infection() counts every case, or with detected_only the reported ones", {
  r <- episodes_fixture_run()
  cs <- classify_cases(episodes_fixture())
  for (det in c(FALSE, TRUE)) {
    fo <- force_of_infection(r, split_imported = FALSE, detected_only = det)
    want <- (if (det) cs[detected == 1L] else cs)[, .(n = .N), by = .(win = as.integer((tick %/% 30) * 30), origin)]
    got <- fo[n > 0, .(win, origin = as.character(origin), n)]
    data.table::setorder(want, win, origin); data.table::setorder(got, win, origin)
    expect_identical(as.data.frame(got), as.data.frame(want), label = paste("detected_only =", det))
  }
  every <- force_of_infection(r, split_imported = FALSE)
  shown <- force_of_infection(r, split_imported = FALSE, detected_only = TRUE)
  ## days 30-59 hold relapses 502, 403 and 105; 105 was not detected
  expect_identical(every[win == 30L & origin == "relapse", n], 3L)
  expect_identical(shown[win == 30L & origin == "relapse", n], 2L)
  expect_identical(shown[, .(win, origin, days, is_acquisition)], every[, .(win, origin, days, is_acquisition)])
  expect_error(force_of_infection(r, detected_only = NA), "`detected_only` must be TRUE or FALSE.", fixed = TRUE)
})

test_that("plot_force_of_infection() leaves out the trend of a series too short for its span", {
  skip_if_not_installed("ggplot2")
  run <- example_run()
  ## 7 windows of 60 days: loess at span 0.4 would fit 2 or 3 points and warn
  short <- plot_force_of_infection(force_of_infection(run, bin = 60))
  expect_no_warning(b <- ggplot2::ggplot_build(short))
  expect_length(b$data, 1L)                      # the raw series only
  ## 14 windows of 30 days carry the trend
  long <- ggplot2::ggplot_build(plot_force_of_infection(run))
  expect_length(long$data, 2L)
  expect_length(ggplot2::ggplot_build(plot_force_of_infection(run, smooth = NA))$data, 1L)
})
