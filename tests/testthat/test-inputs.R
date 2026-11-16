## Every reader of a run's tables takes the run itself, and refuses the other kind of table by name.

test_that("the events and state readers take a run as well as its table", {
  r <- episodes_fixture_run(state = fixture_state())
  ev <- load_events(r); st <- load_state(r)
  expect_identical(classify_travel(r), classify_travel(ev))
  expect_identical(dormancy_durations(r), dormancy_durations(ev))
  expect_identical(transmission_links(transmission_fixture_run()), transmission_links(transmission_fixture()))
  expect_identical(host_state_prevalence(r), host_state_prevalence(st))
  expect_identical(coi_distribution(r), coi_distribution(st))
})

test_that("a reader given the other kind of table says which it got", {
  r <- episodes_fixture_run(state = fixture_state())
  ev <- load_events(r); st <- load_state(r)
  for (f in list(transmission_links, dormancy_durations, classify_travel, episodes))
    expect_error(f(st), "must be an events table from load_events(), a vivax_run or a run directory; this is a state table (from load_state()).", fixed = TRUE)
  for (f in list(host_state_prevalence, coi_distribution))
    expect_error(f(ev), "`state` must be a state table from load_state(), a vivax_run or a run directory; this is an events table (from load_events()).", fixed = TRUE)
  expect_error(summarise_run(ev), "`run` must be a vivax_run (from run_sim() or load_run()) or the path of a run directory, got an events table (from load_events()).", fixed = TRUE)
  expect_error(summarise_run(st), "got a state table (from load_state()).", fixed = TRUE)
  expect_error(summarise_run(42), "got an object of class numeric.", fixed = TRUE)
})

test_that("every plot takes a run, or the table its function returns", {
  skip_if_not_installed("ggplot2")
  r <- episodes_fixture_run(state = fixture_state())
  st <- load_state(r)
  expect_s3_class(plot_relapse_km(r), "ggplot")
  expect_identical(plot_relapse_km(r)$data, plot_relapse_km(time_to_relapse(r))$data)
  expect_error(plot_relapse_km(st), "`x` must be a run, or the table time_to_relapse() returns; it has no `time` or `event`.",
               fixed = TRUE)
  expect_s3_class(plot_force_of_infection(r), "ggplot")
  expect_error(plot_force_of_infection(st), "the table force_of_infection() returns; it has no `win` or `origin` or `rate` or `is_acquisition`.",
               fixed = TRUE)
  for (x in list(r, st, host_state_prevalence(st)))
    expect_identical(plot_host_state_prevalence(x)$data, plot_host_state_prevalence(host_state_prevalence(st))$data)
  expect_error(plot_relatedness_ladder(st), "the table ibd_truth_pairs() returns; it has no `cls` or `ibd`.", fixed = TRUE)
})

test_that("the readers that need the run or its tree sequence say which table they got", {
  r <- episodes_fixture_run(state = fixture_state())
  expect_error(time_to_relapse(load_state(r)),
               "needs the run for its length, SIM_DAYS: pass the vivax_run or its directory, got a state table (from load_state()).",
               fixed = TRUE)
  expect_error(time_to_relapse(load_events(r)), "got an events table (from load_events()).", fixed = TRUE)
  expect_error(malariafwd:::.trees_file(load_state(r)),
               "`run` must be a vivax_run (from run_sim() or load_run()) or the path of a run directory, got a state table (from load_state()).",
               fixed = TRUE)
  skip_if_no_tskit()
  expect_error(ibd_pairs(load_events(r), ids = 1:3), "got an events table (from load_events()).", fixed = TRUE)
  expect_error(arg_summary(load_state(r)), "got a state table (from load_state()).", fixed = TRUE)
})
