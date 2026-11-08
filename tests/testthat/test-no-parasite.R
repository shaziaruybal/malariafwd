## Every reader takes a run with no parasite: its usual columns and classes, no error and no warning.

## (a) never started, no event; (b) never started, imports that did not take; (c) founders that
## infected no one; (d) extinct. Made once for this file.
np_runs <- local({
  runs <- NULL
  function() {
    if (is.null(runs)) {
      mk <- function(tag, ...) suppressMessages(run_sim(vivax_params(CHROMOSOMES = "1", ...), seed = 1,
                                                        outdir = tempfile(paste0("mfwd_np_", tag, "_"))))
      runs <<- list(
        a = never_run(),
        b = mk("b", N_HOSTS = 20, N_MOSQ = 20, SIM_DAYS = 60, N_FOUNDER_HOSTS = 0, IMPORT_TICK = 1, N_IMPORTS = 3,
               P_PRIMARY = 0, NU_HYPNOZOITE_BATCH = 0),
        c = mk("c", N_HOSTS = 20, N_MOSQ = 1, BITING_RATE = 1e-6, SIM_DAYS = 60, N_FOUNDER_HOSTS = 3),
        d = w_run(2))
    }
    runs
  }
})

quiet <- function(expr) expect_no_warning(expect_no_error(expr))

test_that("every reader takes a run with no parasite", {
  skip_if_no_slim()
  rs <- np_runs()
  expect_true(rs$a$never_started && rs$b$never_started)
  expect_true(all(load_events(rs$b)$event_type == "vh_inoculation_import_no_take"))
  for (k in names(rs)) {
    r <- rs[[k]]
    ev <- load_events(r); st <- load_state(r)
    lo <- data.frame(chrom = 1, pos = c(1000, 2000), plaf = 0.3)
    readers <- list(
      load_run = function() load_run(r$outdir), load_params = function() load_params(r),
      classify_cases = function() classify_cases(ev), classify_who = function() classify_who(ev),
      classify_chain = function() classify_chain(r), episodes = function() episodes(r),
      recurrences = function() recurrences(r), relapse_fraction = function() relapse_fraction(r),
      relapses_per_person_year = function() relapses_per_person_year(r),
      dormancy_durations = function() dormancy_durations(ev), reservoir_trajectory = function() reservoir_trajectory(st),
      coi_distribution = function() coi_distribution(st), host_state_prevalence = function() host_state_prevalence(st),
      incidence_by_origin = function() incidence_by_origin(ev), force_of_infection = function() force_of_infection(r),
      transmission_links = function() transmission_links(ev), rc_episode = function() rc_by_origin(r),
      rc_acquisition = function() rc_by_origin(r, unit = "acquisition"), time_to_relapse = function() time_to_relapse(r),
      moi_vs_relapse = function() moi_vs_relapse(st, ev), pick_host = function() pick_host(ev),
      summarise_run = function() summarise_run(r), chrom_map = function() chrom_map(r))
    for (nm in names(readers)) quiet(readers[[nm]]())
    if (k != "d") expect_identical(nrow(classify_cases(ev)), 0L, label = k)
    expect_identical(nrow(summarise_run(r)), 1L, label = k)
    if (k %in% c("a", "b")) expect_identical(nrow(moi_vs_relapse(st, ev)), 0L, label = k)
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      quiet(ggplot2::ggplot_build(plot_force_of_infection(force_of_infection(r), smooth = NA)))
      quiet(ggplot2::ggplot_build(plot_host_state_prevalence(host_state_prevalence(st))))
      if (k %in% c("a", "b", "c")) expect_error(plot_relapse_km(time_to_relapse(r)), "no hypnozoite batch was established")
      if (requireNamespace("patchwork", quietly = TRUE)) quiet(print(plot_diagnostics(r)))
    }
  }
})

test_that("the tree-sequence readers take a run with no parasite", {
  skip_if_no_slim(); skip_if_no_tskit()
  rs <- np_runs()
  for (k in names(rs)) {
    r <- rs[[k]]
    lo <- data.frame(chrom = 1, pos = c(1000, 2000), plaf = 0.3)
    quiet(ibd_pairs(r)); quiet(ibd_segments(r, integer(), integer())); quiet(ibd_truth_pairs(r))
    quiet(new_infection_diversity(r)); quiet(arg_summary(r)); quiet(arg_tmrca(r, integer(), integer()))
    quiet(simulate_genotypes(r, integer(), lo))
    expect_identical(new_infection_diversity(r)$n, 0L, label = k)
    if (k %in% c("a", "b", "c")) expect_identical(nrow(ibd_truth_pairs(r)), 0L, label = k)
  }
  expect_identical(arg_summary(rs$a)$num_samples, 0L)
  expect_identical(arg_summary(rs$b)$num_samples, 0L)
})
