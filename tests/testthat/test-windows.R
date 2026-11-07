## Every run-length window is SIM_DAYS long, so the days after an extinction count as days with no
## parasite.

test_that("summarise_run() counts the days after an extinction as days with no infected host", {
  skip_if_no_slim()
  for (s in 1:6) {
    r <- w_run(s)
    expect_lt(r$extinct_tick, 200L)
    m <- suppressWarnings(summarise_run(r))
    expect_identical(c(m$prev, m$mean_res), c(0, 0), label = paste("seed", s))
    expect_true(is.na(m$COI) && is.na(m$poly), label = paste("seed", s))
  }
})

test_that("force_of_infection() gives each window its days within SIM_DAYS", {
  skip_if_no_slim()
  hit <- FALSE
  for (s in 1:6) {
    r <- w_run(s)
    fo <- force_of_infection(r)[win > 0]
    expect_identical(fo$days, as.numeric(pmin(30L, 400L - fo$win + 1L)), label = paste("seed", s))
    hit <- hit || any(fo$win == (r$extinct_tick %/% 30) * 30)
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

test_that("a run that reaches SIM_DAYS gives the same summary as before", {
  skip_if_no_slim()
  r <- suppressMessages(run_sim(vivax_params(N_HOSTS = 30, N_MOSQ = 60, SIM_DAYS = 40), seed = 4,
                                outdir = tempfile("mfwd_full_")))
  st <- load_state(r)
  m <- suppressWarnings(summarise_run(r))
  expect_identical(m$prev, round(mean(st[tick >= 20]$moi_active > 0), 3))
})

test_that(".sim_days() reads SIM_DAYS, and refuses an events table", {
  expect_identical(malariafwd:::.sim_days(fake_run(params = list(SIM_DAYS = 100))), 100L)
  r <- fake_run(drop = "SIM_DAYS", extra = "SIM_DAYS\t400.0")
  expect_identical(malariafwd:::.sim_days(r), 400L)
  expect_error(malariafwd:::.sim_days(load_events(r), "summarise_run"), "summarise_run\\(\\) needs the run for its length")
})
