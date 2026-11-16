## Relapse metrics: relapse_fraction() and each batch's time to relapse.

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

test_that("time_to_relapse() times each batch from its establishment to its first relapse", {
  tt <- time_to_relapse(episodes_fixture_run())
  expect_identical(names(tt), c("batch_pedigree_id", "human_id", "lineage_id", "imported", "batch_tick",
                                "n_hypnozoites", "time", "event", "end_reason"))
  expect_identical(tt$batch_pedigree_id, c(501L, 402L, 103L, 203L, 504L))
  expect_identical(tt$human_id, c(5L, 4L, 1L, 2L, 5L))
  expect_identical(tt$lineage_id, c("7", "5", "1", "3", "8"))
  expect_identical(tt$imported, c(FALSE, FALSE, FALSE, FALSE, TRUE))
  expect_identical(tt$batch_tick, c(7L, 8L, 10L, 25L, 61L))
  expect_identical(tt$n_hypnozoites, c(2L, 1L, 2L, 3L, 1L))
  ## 203 is emptied by radical cure on day 33; 504 still holds its hypnozoite on day 100
  expect_identical(tt$time, c(23L, 37L, 10L, 8L, 39L))
  expect_identical(tt$event, c(1L, 1L, 1L, 0L, 0L))
  expect_identical(tt$end_reason, c("relapse", "relapse", "relapse", "emptied", "run_end"))
})

test_that("an emptied batch can be followed to the run's end, and a masked relapse passed over", {
  r <- episodes_fixture_run()
  f <- time_to_relapse(r, emptied = "follow")
  expect_identical(f[batch_pedigree_id == 203L, .(time, event, end_reason)],
                   data.table::data.table(time = 75L, event = 0L, end_reason = "emptied"))
  expect_identical(f[batch_pedigree_id != 203L], time_to_relapse(r)[batch_pedigree_id != 203L])
  m <- time_to_relapse(r, exclude_masked = TRUE)
  expect_identical(m[batch_pedigree_id == 103L]$time, 40L)   # activation 105 on day 50
  expect_identical(m[batch_pedigree_id != 103L], time_to_relapse(r)[batch_pedigree_id != 103L])
  expect_error(time_to_relapse(r, emptied = "drop"), "'arg' should be one of")
  expect_error(time_to_relapse(r, exclude_masked = NA), "`exclude_masked` must be TRUE or FALSE")
})

test_that("plot_relapse_km() starts at 1 on day 0 and steps at each relapse", {
  skip_if_not_installed("ggplot2")
  p <- plot_relapse_km(time_to_relapse(episodes_fixture_run()))
  expect_identical(as.numeric(p$data$time), c(0, 8, 10, 23, 37, 39))
  expect_equal(p$data$surv, c(1, 1, 0.75, 0.5, 0.25, 0.25))
})

test_that("a run with no batch gives no row, and the plot says why it cannot draw", {
  r <- episodes_fixture_run(events = episodes_fixture()[0])
  expect_no_warning(tt <- time_to_relapse(r))
  expect_identical(nrow(tt), 0L)
  expect_identical(names(tt), names(time_to_relapse(episodes_fixture_run())))
  expect_identical(vapply(tt, function(v) class(v)[1], ""),
                   vapply(time_to_relapse(episodes_fixture_run()), function(v) class(v)[1], ""))
  expect_error(time_to_relapse(episodes_fixture()), "needs the run for its length")
  skip_if_not_installed("ggplot2")
  expect_error(plot_relapse_km(tt), "no hypnozoite batch was established in this run")
  expect_error(plot_relapse_km(data.frame(x = 1)), "it has no `time` or `event`.", fixed = TRUE)
})

test_that("on an import cohort every batch is timed from its establishment, and radical cure empties batches", {
  skip_if_no_slim()
  cohort <- function(rc) suppressMessages(run_sim(vivax_params(
    N_HOSTS = 1000, N_MOSQ = 1, BITING_RATE = 1e-6, SIM_DAYS = 400, N_FOUNDER_HOSTS = 0, IMPORT_TICK = 5,
    N_IMPORTS = 800, BLOODCURE_COVERAGE = 1, P_DETECT_RELAPSE = 0, RADCURE_COVERAGE = rc, CHROMOSOMES = "1"),
    seed = 1, outdir = tempfile("mfwd_ttr_")))
  r <- cohort(1); ev <- load_events(r); tt <- time_to_relapse(r)
  est <- ev[event_type == "vh_inoculation_import" & is_dormant == 1,
            .(batch_pedigree_id = as.integer(pedigree_id), tick = as.integer(tick))]
  expect_gt(nrow(est), 100)
  expect_setequal(tt$batch_pedigree_id, est$batch_pedigree_id)
  expect_identical(tt[est, on = "batch_pedigree_id"]$batch_tick, est$tick)
  ## a batch radical cure reached before any of its hypnozoites activated
  act <- ev[event_type == "hypnozoite_activation", .(t_act = min(tick)), by = .(batch_pedigree_id = parent1_pedigree_id)]
  rad <- ev[event_type == "radcure_hypnozoite_clear", .(t_rad = min(tick)), by = .(batch_pedigree_id = pedigree_id)]
  pre <- act[rad, on = "batch_pedigree_id"][is.na(t_act) | t_rad < t_act]
  expect_gt(nrow(pre), 0)
  got <- tt[pre, on = "batch_pedigree_id"]
  expect_true(all(got$end_reason == "emptied"))
  expect_identical(got$time, as.integer(got$t_rad - got$batch_tick))
  ev0 <- load_events(cohort(0))
  expect_false("radcure_hypnozoite_clear" %in% ev0$event_type)
})

test_that("a batch that empties before relapsing in a run that goes extinct is censored where it emptied", {
  skip_if_no_slim()
  r <- suppressMessages(run_sim(vivax_params(N_HOSTS = 20, N_MOSQ = 20, N_FOUNDER_HOSTS = 0, CHROMOSOMES = "1",
                                             SIM_DAYS = 60, IMPORT_TICK = 1, N_IMPORTS = 1, P_PRIMARY = 0,
                                             NU_HYPNOZOITE_BATCH = 1e6, MU_HYPNOZOITE = 1, ALPHA_ACTIVATE = 0),
                                seed = 1, outdir = tempfile("mfwd_ttr_ext_")))
  expect_identical(r$extinct_tick, 2L)
  tt <- time_to_relapse(r)
  expect_identical(tt[, .(batch_tick, time, event, end_reason)],
                   data.table::data.table(batch_tick = 1L, time = 1L, event = 0L, end_reason = "emptied"))
})

test_that("every liver-stage row that names a genome takes hypnozoites out of its batch", {
  skip_if_no_slim()
  ## one genotype per bite: at this seed the run then logs every way a hypnozoite leaves its batch
  r <- suppressMessages(run_sim(vivax_params(N_HOSTS = 100, N_MOSQ = 200, BITING_RATE = 0.11, SIM_DAYS = 200,
                                             BLOODCURE_COVERAGE = 1, RADCURE_COVERAGE = 1, RADCURE_MODE = "hypnozoite",
                                             MDA_TICK = 100, CHROMOSOMES = "1",
                                             LAMBDA_OOCYSTS = 0, LAMBDA_SPOROZOITES = 0),
                                seed = 1, outdir = tempfile("mfwd_ttr_liver_")))
  ev <- load_events(r)
  liver <- event_types()$event_type[event_types()$stage == "liver"]
  named <- unique(ev[event_type %in% liver & pedigree_id != -1]$event_type)
  expect_true(all(c("hypnozoite_clearance", "radcure_hypnozoite_clear", "mda_hypnozoite_clear") %in% named))
  expect_setequal(named, malariafwd:::.hypnozoite_exit_events())
})

test_that("dormancy_durations() gives each relapse its host, batch and both days", {
  d <- dormancy_durations(episodes_fixture())
  expect_identical(names(d), c("human_id", "batch_pedigree_id", "pedigree_id", "batch_tick", "tick", "duration"))
  expect_identical(d$human_id, c(5L, 4L, 1L, 1L))
  expect_identical(d$batch_pedigree_id, c(501L, 402L, 103L, 103L))
  expect_identical(d$pedigree_id, c(502L, 403L, 104L, 105L))
  expect_identical(d$batch_tick, c(7L, 8L, 10L, 10L))
  expect_identical(d$tick, c(30L, 45L, 20L, 50L))
  expect_identical(d$duration, c(23L, 37L, 10L, 40L))
  e <- dormancy_durations(make_events())
  expect_identical(nrow(e), 0L)
  expect_true(all(vapply(e, is.integer, logical(1))))
})
