## summarise_run(): one window of days for every column, in typed columns named by their units.

sm_cols <- c("from_tick", "to_tick", "n_days", "exclude_masked", "extinct", "extinct_tick", "never_started",
             "prevalence_frac", "coi_mean", "coi_median", "polyclonal_frac", "hypnozoites_per_host",
             "eir_per_py", "foi_per_py", "molfob_clones_per_py", "molfob_episodes_per_py",
             "relapses_per_py_at_risk", "relapse_frac_episodes", "relapse_frac_acquisitions", "imported_frac",
             "n_imported", "n_local_reinfection", "n_relapse", "n_relapse_masked", "n_recrudescence",
             "n_bloodstage_courses", "person_years", "person_years_at_risk")
sm_class <- c(rep("integer", 3), "logical", "logical", "integer", "logical", rep("numeric", 13),
              rep("integer", 6), rep("numeric", 2))
classes <- function(x) unname(vapply(x, function(v) class(v)[1], ""))
sm_main <- function() episodes_fixture_run(state = fixture_state())

test_that("summarise_run() gives one row of named, typed columns, the same for every run", {
  rows <- list(summarise_run(sm_main()), summarise_run(extinct_fixture_run()), summarise_run(never_fixture_run()))
  for (m in rows) {
    expect_identical(names(m), sm_cols)
    expect_identical(classes(m), sm_class)
    expect_identical(nrow(m), 1L)
  }
  expect_identical(classes(data.table::rbindlist(rows)), sm_class)
})

test_that("summarise_run() over the whole run", {
  m <- summarise_run(sm_main(), burn_frac = 0)
  py <- 600 / 365
  expect_identical(c(m$from_tick, m$to_tick, m$n_days), c(1L, 100L, 100L))
  expect_identical(c(m$exclude_masked, m$extinct, m$never_started), c(FALSE, FALSE, FALSE))
  expect_identical(m$extinct_tick, NA_integer_)
  expect_equal(c(m$prevalence_frac, m$coi_mean, m$coi_median, m$polyclonal_frac, m$hypnozoites_per_host),
               c(277 / 600, 347 / 277, 1, 70 / 277, 198 / 600))
  expect_equal(c(m$eir_per_py, m$foi_per_py, m$molfob_clones_per_py, m$molfob_episodes_per_py), c(9, 7, 11, 12) / py)
  expect_equal(m$relapses_per_py_at_risk, 4 / (445 / 365))
  expect_equal(c(m$relapse_frac_episodes, m$relapse_frac_acquisitions, m$imported_frac), c(4 / 12, 4 / 11, 1 / 7))
  expect_identical(c(m$n_imported, m$n_local_reinfection, m$n_relapse, m$n_relapse_masked, m$n_recrudescence,
                     m$n_bloodstage_courses), c(1L, 6L, 4L, 1L, 1L, 2L))
  expect_equal(c(m$person_years, m$person_years_at_risk), c(600, 445) / 365)
})

test_that("the default window is the second half, and every column follows the window", {
  r <- sm_main()
  m <- summarise_run(r)
  py <- 300 / 365
  expect_identical(c(m$from_tick, m$to_tick, m$n_days), c(51L, 100L, 50L))
  expect_equal(c(m$prevalence_frac, m$coi_mean, m$polyclonal_frac, m$hypnozoites_per_host),
               c(150 / 300, 189 / 150, 39 / 150, 40 / 300))
  expect_equal(c(m$eir_per_py, m$foi_per_py, m$molfob_clones_per_py, m$molfob_episodes_per_py), c(2, 1, 1, 1) / py)
  expect_identical(m$relapses_per_py_at_risk, 0)
  expect_equal(m$person_years_at_risk, 250 / 365)
  expect_identical(c(m$relapse_frac_episodes, m$relapse_frac_acquisitions, m$imported_frac), c(0, 0, 1))
  m <- summarise_run(r, from = 41)
  py <- 360 / 365
  expect_equal(c(m$prevalence_frac, m$coi_mean, m$polyclonal_frac, m$hypnozoites_per_host),
               c(167 / 360, 212 / 167, 45 / 167, 53 / 360))
  expect_equal(c(m$eir_per_py, m$foi_per_py, m$molfob_clones_per_py), c(4, 2, 4) / py)
  expect_equal(m$relapses_per_py_at_risk, 2 / (300 / 365))
  expect_equal(c(m$relapse_frac_episodes, m$relapse_frac_acquisitions, m$imported_frac), c(2 / 4, 2 / 4, 1 / 2))
  m <- summarise_run(r, from = 1, to = 30)
  py <- 180 / 365
  expect_equal(m$prevalence_frac, 80 / 180)
  expect_equal(c(m$eir_per_py, m$molfob_clones_per_py, m$molfob_episodes_per_py), c(5, 7, 8) / py)
  expect_equal(m$relapses_per_py_at_risk, 2 / (95 / 365))
  expect_identical(m$n_bloodstage_courses, 1L)
  expect_identical(summarise_run(r, from = 0), summarise_run(r, from = 1))
})

test_that("exclude_masked leaves the masked relapse out of every relapse column, and only those", {
  r <- sm_main()
  a <- summarise_run(r, burn_frac = 0); m <- summarise_run(r, burn_frac = 0, exclude_masked = TRUE)
  expect_identical(c(m$n_relapse, m$n_relapse_masked, a$n_relapse_masked), c(3L, 1L, 1L))
  expect_equal(c(m$molfob_clones_per_py, m$molfob_episodes_per_py), c(10, 11) * 365 / 600)
  expect_equal(m$relapses_per_py_at_risk, 3 / (445 / 365))
  expect_equal(c(m$relapse_frac_episodes, m$relapse_frac_acquisitions), c(3 / 11, 3 / 10))
  expect_true(m$exclude_masked)
  same <- c("prevalence_frac", "coi_mean", "coi_median", "polyclonal_frac", "hypnozoites_per_host", "eir_per_py",
            "foi_per_py", "imported_frac", "person_years", "person_years_at_risk")
  expect_identical(m[, same, with = FALSE], a[, same, with = FALSE])
  ## an event log without the flag, and a state table without the bite count
  nf <- episodes_fixture_run(events = episodes_fixture()[, into_patent_lineage := NULL], state = fixture_state())
  x <- summarise_run(nf, burn_frac = 0)
  expect_identical(x$n_relapse_masked, NA_integer_)
  expect_equal(x$relapse_frac_acquisitions, 4 / 11)
  expect_error(summarise_run(nf, burn_frac = 0, exclude_masked = TRUE), "into_patent_lineage")
  nb <- episodes_fixture_run(state = fixture_state()[, infectious_bites := NULL])
  expect_identical(summarise_run(nb, burn_frac = 0)$eir_per_py, NA_real_)
})

test_that("the days after an extinction count as days with no parasite", {
  r <- extinct_fixture_run()
  m <- summarise_run(r, burn_frac = 0)
  expect_equal(c(m$prevalence_frac, m$coi_mean, m$polyclonal_frac), c(4 / 40, 1, 0))
  expect_equal(c(m$eir_per_py, m$foi_per_py), c(1, 1) * 365 / 40)
  expect_identical(m$relapses_per_py_at_risk, 0)
  expect_equal(m$person_years_at_risk, 18 / 365)
  expect_identical(c(m$relapse_frac_episodes, m$relapse_frac_acquisitions), c(0, 0))
  m <- summarise_run(r)
  expect_identical(c(m$from_tick, m$to_tick), c(11L, 20L))
  expect_identical(c(m$prevalence_frac, m$hypnozoites_per_host, m$eir_per_py, m$molfob_clones_per_py), rep(0, 4))
  expect_identical(c(m$coi_mean, m$coi_median, m$polyclonal_frac), rep(NA_real_, 3))
  expect_identical(c(m$relapse_frac_episodes, m$relapse_frac_acquisitions, m$imported_frac), rep(NA_real_, 3))
  expect_identical(m$relapses_per_py_at_risk, 0)
  expect_equal(m$person_years_at_risk, 10 / 365)
  expect_true(m$extinct); expect_false(m$never_started)
  expect_identical(m$extinct_tick, 6L)
})

test_that("a run that never started gives zeros and NA, without a warning", {
  expect_no_warning(m <- summarise_run(never_fixture_run()))
  expect_true(m$never_started); expect_false(m$extinct)
  expect_identical(m$extinct_tick, NA_integer_)
  expect_identical(c(m$prevalence_frac, m$hypnozoites_per_host, m$eir_per_py, m$foi_per_py, m$molfob_clones_per_py),
                   rep(0, 5))
  expect_identical(c(m$coi_mean, m$relapse_frac_acquisitions, m$imported_frac, m$relapses_per_py_at_risk),
                   rep(NA_real_, 4))
  expect_identical(m$person_years_at_risk, 0)
  expect_identical(c(m$n_imported, m$n_local_reinfection, m$n_relapse, m$n_relapse_masked, m$n_recrudescence,
                     m$n_bloodstage_courses), rep(0L, 6))
})

test_that("a run that stopped without a record is refused past its last day", {
  r <- extinct_fixture_run(record = FALSE)
  expect_error(summarise_run(r), "stops on day 6 of 20")
  expect_error(summarise_run(r), "`to = 6`", fixed = TRUE)
  expect_equal(summarise_run(r, from = 1, to = 6)$prevalence_frac, 4 / 12)
  st <- fixture_state()[!(tick == 60L & entity_type == "human" & entity_id == 3L)]
  expect_error(summarise_run(episodes_fixture_run(state = st)), "one row per host")
})

test_that("the window's arguments are checked", {
  r <- sm_main()
  expect_error(summarise_run(r, burn_frac = 0.5, from = 10), "not both")
  for (b in list(1, -0.1, NA, "0.5")) expect_error(summarise_run(r, burn_frac = b), "`burn_frac` must be")
  for (f in list(1.5, -1, c(1, 2))) expect_error(summarise_run(r, from = f), "`from` must be")
  expect_error(summarise_run(r, to = 101), "SIM_DAYS = 100")
  expect_error(summarise_run(r, from = 60, to = 50), "has no day")
  expect_error(summarise_run(r, to = 40), "has no day")
  expect_error(summarise_run(r, exclude_masked = NA), "`exclude_masked` must be TRUE or FALSE")
})

test_that("summarise_run() takes a run or its directory", {
  r <- sm_main()
  expect_identical(summarise_run(r), summarise_run(r$outdir))
  expect_error(summarise_run(data.table::data.table()), "`run` must be a vivax_run")
})

test_that("on a run, every column is its definition over the window", {
  skip_if_no_slim()
  r <- treated_run(); ev <- load_events(r); st <- load_state(r); nh <- 60
  for (b in c(0, 0.5)) {
    m <- summarise_run(r, burn_frac = b)
    f <- m$from_tick; t <- m$to_tick
    expect_identical(c(f, t), c(as.integer(round(400 * b) + 1), 400L))
    h <- st[tick >= f & tick <= t]; hd <- nh * (t - f + 1)
    expect_equal(c(m$prevalence_frac, m$hypnozoites_per_host), c(sum(h$moi_active > 0), sum(h$hyp_reservoir_size)) / hd)
    expect_equal(m$coi_mean, mean(h[moi_active > 0]$moi_active))
    expect_equal(m$eir_per_py, sum(h$infectious_bites) / (hd / 365))
    evw <- ev[tick >= f & tick <= t]; cs <- classify_cases(evw)
    expect_identical(c(m$n_imported, m$n_local_reinfection, m$n_relapse, m$n_recrudescence),
                     c(cs[origin == "imported", .N], cs[origin == "local reinfection", .N],
                       cs[origin == "relapse", .N], cs[origin == "recrudescence", .N]))
    expect_identical(m$n_relapse_masked, sum(cs[origin == "relapse"]$into_patent_lineage))
    expect_identical(m$n_bloodstage_courses, evw[event_type %in% c("bloodcure_given", "mda_bloodcure_given"), .N])
    rf <- relapse_fraction(evw)
    expect_identical(c(m$relapse_frac_episodes, m$relapse_frac_acquisitions),
                     c(rf$relapse_frac_episodes, rf$relapse_frac_acquisitions))
    rp <- relapses_per_person_year(r, from = f, to = t)
    expect_identical(c(m$relapses_per_py_at_risk, m$person_years_at_risk), c(rp$relapses_per_py, rp$person_years))
  }
  a <- summarise_run(r, burn_frac = 0); b <- summarise_run(r)
  for (col in c("relapse_frac_acquisitions", "relapses_per_py_at_risk", "imported_frac", "n_bloodstage_courses"))
    expect_false(isTRUE(all.equal(a[[col]], b[[col]])), label = col)
})

test_that("FOI and molFOB are force_of_infection()'s counts over the whole run", {
  skip_if_no_slim()
  r <- treated_run()
  for (x in c(FALSE, TRUE)) {
    m <- summarise_run(r, burn_frac = 0, exclude_masked = x)
    fo <- force_of_infection(r, exclude_masked = x)
    py <- m$person_years
    expect_equal(m$foi_per_py * py, sum(fo[origin %in% c("imported", "local reinfection")]$n))
    expect_equal(m$molfob_clones_per_py * py, sum(fo[is_acquisition == TRUE]$n))
    expect_equal(m$molfob_episodes_per_py * py, sum(fo$n))
    expect_equal((m$molfob_episodes_per_py - m$molfob_clones_per_py) * py, m$n_recrudescence)
  }
  expect_gt(m$n_recrudescence, 0L)
})

test_that("runs that went extinct or never started give the same typed row", {
  skip_if_no_slim()
  w <- w_run(2)
  m <- summarise_run(w, from = w$extinct_tick + 1)
  expect_identical(c(m$prevalence_frac, m$eir_per_py, m$molfob_clones_per_py), c(0, 0, 0))
  expect_identical(m$coi_mean, NA_real_)
  a <- summarise_run(w, burn_frac = 0)
  expect_equal(a$prevalence_frac, sum(load_state(w)$moi_active > 0) / (60 * 400))
  expect_identical(a$extinct_tick, w$extinct_tick)
  expect_no_warning(n <- summarise_run(never_run()))
  expect_true(n$never_started)
  expect_identical(classes(a), classes(n))
})

test_that("plot_diagnostics() gives summarise_run()'s numbers for the days after the burn-in", {
  skip_if_not_installed("ggplot2"); skip_if_not_installed("patchwork")
  r <- sm_main()
  sub <- function(p) p$patches$annotation$subtitle
  p <- plot_diagnostics(r)
  expect_match(sub(p), paste("days 51-100: prevalence 50% | COI 1.26 | EIR 2.4/py | FOI 1.2/py | molFOB 1.2/py",
                             "| relapse frac (episodes) 0.00"), fixed = TRUE)
  expect_match(sub(p), "burn-in 50 d (dashed)", fixed = TRUE)
  m <- summarise_run(r, burn_frac = 0)
  line <- sprintf(paste("days %d-%d: prevalence %.0f%% | COI %.2f | EIR %.1f/py | FOI %.1f/py | molFOB %.1f/py",
                        "| relapse frac (episodes) %.2f"), m$from_tick, m$to_tick, 100 * m$prevalence_frac,
                  m$coi_mean, m$eir_per_py, m$foi_per_py, m$molfob_clones_per_py, m$relapse_frac_episodes)
  expect_match(sub(plot_diagnostics(r, burn_frac = 0)), line, fixed = TRUE)
  ## the day-22 recrudescence is a new episode, not a new clone
  expect_true(any(p[[2]]$data[origin == "recrudescence"]$y > 0))
  expect_false(any(p[[3]]$data[origin == "recrudescence"]$y > 0))
  expect_identical(p[[5]]$labels$y, "Mean hypnozoites/host")
  for (x in list(extinct_fixture_run(), never_fixture_run())) {
    s <- sub(plot_diagnostics(x, bin = 5))
    expect_match(s, "days 11-20: prevalence 0% | COI NA", fixed = TRUE)
    expect_match(s, "(episodes) NA", fixed = TRUE)
  }
  expect_error(plot_diagnostics(extinct_fixture_run(record = FALSE), bin = 5),
               "plot_diagnostics() needs a finished run: state.tsv", fixed = TRUE)
})

test_that("on a run, plot_diagnostics() names its window and that window's numbers", {
  skip_if_no_slim(); skip_if_not_installed("ggplot2"); skip_if_not_installed("patchwork")
  r <- treated_run()
  m <- summarise_run(r, burn_frac = 0.9)
  s <- plot_diagnostics(r, burn_frac = 0.9)$patches$annotation$subtitle
  expect_match(s, sprintf("days 361-400: prevalence %.0f%% | COI %.2f", 100 * m$prevalence_frac, m$coi_mean), fixed = TRUE)
  expect_match(s, sprintf("relapse frac (episodes) %.2f", m$relapse_frac_episodes), fixed = TRUE)
})

test_that("the dashboard's force panel is force_of_infection()'s rate, its short first window included", {
  skip_if_not_installed("ggplot2"); skip_if_not_installed("patchwork")
  r <- sm_main()
  d <- plot_diagnostics(r)[[3]]$data
  fo <- force_of_infection(r, split_imported = FALSE)[is_acquisition == TRUE & win %in% d$w]
  expect_true(any(fo$win == 0L & fo$days < 30L & fo$rate > 0))
  expect_equal(d[fo, on = .(w = win, origin), y], fo$rate)
})

test_that("the dashboard draws the run's last whole window and leaves out a part-window after it", {
  skip_if_not_installed("ggplot2"); skip_if_not_installed("patchwork")
  r <- extinct_fixture_run()   # 20 days
  last_window <- function(bin) max(plot_diagnostics(r, bin = bin)[[1]]$data$w)
  expect_equal(last_window(7), 14)   # days 14 to 20 are a whole window
  expect_equal(last_window(5), 15)   # days 15 to 19; day 20 alone is left out
})
