## The exported functions, listed: one added, removed or renamed changes this list on purpose.

test_that("the package exports these functions", {
  ns <- readLines(system.file("NAMESPACE", package = "malariafwd"))
  exports <- sub("^export\\((.*)\\)$", "\\1", grep("^export\\(", ns, value = TRUE))
  expect_setequal(exports, c(
  "arg_divergence", "arg_diversity", "arg_summary", "arg_tmrca", "as_slim_args", "check_setup",
  "chrom_map", "classify_cases", "classify_who", "coi_distribution", "dormancy_durations",
  "event_types", "force_of_infection", "host_state_prevalence", "ibd_lookup", "ibd_pairs",
  "ibd_relapse_origin", "ibd_segments", "ibd_trajectory", "ibd_truth_pairs", "incidence_by_origin",
  "load_events", "load_grid", "load_params", "load_run", "load_state", "model_version", "moi_vs_relapse",
  "new_infection_diversity", "param_info", "pick_host", "plot_diagnostics", "plot_force_of_infection",
  "plot_host_state_prevalence", "plot_host_trajectory", "plot_ibd_pair", "plot_moi_vs_relapse",
  "plot_relapse_km", "plot_relatedness_ladder", "rc_by_origin", "relapse_fraction",
  "relapses_per_person_year", "reservoir_trajectory", "run_grid", "run_sim", "slim_available",
  "slim_model", "slim_path", "slim_version", "summarise_run", "theme_malariafwd",
  "time_to_relapse", "transmission_chains", "transmission_links", "trees_path", "tskit_available",
  "uptake_calibrate", "uptake_moments", "vivax_params"))
})
