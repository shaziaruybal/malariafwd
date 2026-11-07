## The run loaders. Most of these write a small run directory by hand (helper-runs.R), so they
## need no SLiM.

classes <- function(d) vapply(d, function(x) class(x)[1], "")
ev_classes <- c(tick = "integer", event_id = "character", event_type = "character",
                inoculation_id = "integer", mosq_id = "integer", human_id = "integer",
                lineage_id = "integer", parent1_lineage_id = "integer", parent2_lineage_id = "integer",
                is_recombinant = "integer", p_outcross = "numeric", recipient = "character",
                is_dormant = "integer", is_imported = "integer", pedigree_id = "integer",
                parent1_pedigree_id = "integer", parent2_pedigree_id = "integer",
                n_hypnozoites = "integer", oocyst_id = "integer", detected = "integer",
                into_patent_lineage = "integer")
st_classes <- c(tick = "integer", entity_type = "character", entity_id = "integer", tag = "integer",
                moi_active = "integer", moi_subpatent = "integer", moi_dormant = "integer",
                hyp_reservoir_size = "integer", lineage_ids = "character", n_transmissible = "integer",
                mosq_age_days = "integer", infectious_bites = "integer")

test_that("an empty event log has every column, in its class", {
  ev <- load_events(fake_run())
  expect_equal(nrow(ev), 0L)
  expect_identical(classes(ev), ev_classes)
})

test_that("the classes do not depend on what a run holds", {
  ## course rows only, with no lineage, and hosts that never carry two lineages, or any
  r <- fake_run(ev_rows = "5\t5_1\tbloodcure_given\t-1\t-1\t3\tNA\tNA\tNA\t0\t0.0\thuman\t-1\t0\t-1\t-1\t-1\t1\t-1\t-1\t-1",
                st_rows = c("1\thuman\t0\t0\t0\t0\t0\t0\tNA\tNA\tNA\t0", "1\tmosquito\t2\t0\t0\t0\t0\t0\tNA\t0\t1\tNA"))
  expect_identical(classes(load_events(r)), ev_classes)
  expect_identical(classes(load_state(r)), st_classes)
  r <- fake_run(st_rows = "1\thuman\t0\t1\t1\t0\t0\t0\t7\tNA\tNA\t0")
  expect_identical(load_state(r)$lineage_ids, "7")
})

test_that("a file without a column gets it, empty, and an extra column follows", {
  r <- fake_run(ev_head = c(ev_cols[-21], "extra"), ev_rows = paste(c(rep("1", 20), "x"), collapse = "\t"),
                st_head = sub("lineage_ids", "brood_ids", st_cols[-c(6, 12)]),
                st_rows = "1\thuman\t0\t1\t1\t0\t0\t7|9\tNA\tNA")
  ev <- load_events(r); st <- load_state(r)
  expect_identical(names(ev), c(ev_cols, "extra"))
  expect_identical(ev$into_patent_lineage, NA_integer_)
  expect_identical(names(st), st_cols)
  expect_identical(classes(st), st_classes)
  expect_identical(st$lineage_ids, "7|9")
  expect_identical(c(st$moi_subpatent, st$infectious_bites), c(NA_integer_, NA_integer_))
})

test_that("a value that does not fit its column stops, naming it", {
  r <- fake_run(ev_rows = c(paste(rep("1", 21), collapse = "\t"),
                            paste(c("2", "a", "b", "1", "1", "x", rep("1", 15)), collapse = "\t")))
  expect_error(load_events(r), "events.tsv in .*: column human_id holds values that are not integers \\(line 3")
})

test_that("a missing or empty file, or something that is not a run, is named", {
  d <- tempfile("mfwd_none_"); dir.create(d)
  expect_error(load_events(d), "no events.tsv in")
  r <- fake_run(); writeLines(character(), file.path(r, "state.tsv"))
  expect_error(load_state(r), "state.tsv in .* is empty")
  expect_error(load_events(data.table::data.table()), "`run` must be a vivax_run")
})

test_that("normal, extinct and never-started runs load with the same columns and classes", {
  skip_if_no_slim()
  base <- list(N_HOSTS = 20, N_MOSQ = 40, SIM_DAYS = 40, CHROMOSOMES = "1")
  run <- function(...) suppressMessages(run_sim(do.call(vivax_params, utils::modifyList(base, list(...))),
                                                seed = 1, outdir = tempfile("mfwd_load_")))
  runs <- list(run(), run(NU_HYPNOZOITE_BATCH = 0, BITING_RATE = 0.001, MEAN_BLOOD_STAGE_DURATION = 1),
               run(N_FOUNDER_HOSTS = 0, IMPORT_START = 1, IMPORT_END = 20, IMPORT_RATE_PER_TICK = 1e-6))
  expect_false(is.na(runs[[2]]$extinct_tick)); expect_true(runs[[3]]$never_started)
  for (r in runs) {
    expect_identical(classes(load_events(r)), ev_classes)
    expect_identical(classes(load_state(r)), st_classes)
  }
  ## the files' own headers are the loaders' columns
  expect_identical(names(data.table::fread(file.path(runs[[1]]$outdir, "events.tsv"), nrows = 0)), ev_cols)
  expect_identical(names(data.table::fread(file.path(runs[[1]]$outdir, "state.tsv"), nrows = 0)), st_cols)
  ## a run with no event joins as any other
  expect_equal(nrow(moi_vs_relapse(load_state(runs[[3]]), load_events(runs[[3]]))), 0L)
})

test_that("load_params() reads each parameter into its type, and text as it is", {
  r <- fake_run(drop = c("VERBOSE", "RETAIN_FULL_HISTORY", "N_HOSTS", "N_MOSQ", "INTERVENTION_START",
                         "ALPHA_ACTIVATE", "SIM_ID"),
                extra = c("VERBOSE\tT", "RETAIN_FULL_HISTORY\tFALSE", "N_HOSTS\t20.0", "N_MOSQ\t20",
                          "INTERVENTION_START\t3000000000", "ALPHA_ACTIVATE\t0.0058823529411764",
                          "SIM_ID\tNA"))
  p <- load_params(r)
  expect_s3_class(p, "vivax_params")
  expect_identical(names(p), names(vivax_params()))
  expect_identical(p$VERBOSE, TRUE)
  expect_identical(p$RETAIN_FULL_HISTORY, FALSE)
  expect_identical(c(p$N_HOSTS, p$N_MOSQ), c(20L, 20L))
  expect_identical(p$INTERVENTION_START, 3e9)
  expect_identical(p$ALPHA_ACTIVATE, 0.0058823529411764)
  expect_identical(p$SIM_ID, "NA")
  expect_identical(load_params(fake_run(drop = "SIM_ID", extra = "SIM_ID\t x\ty"))$SIM_ID, " x\ty")
  ## every parameter of the set as written, and nothing else
  expect_identical(load_params(fake_run(params = list(N_HOSTS = 30, BITING_RATE = 0.12))),
                   vivax_params(N_HOSTS = 30, BITING_RATE = 0.12))
})

test_that("a params.tsv that cannot be read stops, naming the line or the parameter", {
  expect_error(load_params(fake_run(drop = "N_HOSTS", extra = "N_HOSTS\tabc")),
               "params.tsv in .*: N_HOSTS is \"abc\", not a number")
  expect_error(load_params(fake_run(drop = "VERBOSE", extra = "VERBOSE\tyes")), "VERBOSE is \"yes\", not T or F")
  expect_error(load_params(fake_run(extra = "N_HOSTS\t20")), "N_HOSTS is given twice")
  expect_error(load_params(fake_run(extra = c("", "a line"))), "line [0-9]+ has no tab")
  r <- fake_run(); writeLines("x\ty", file.path(r, "params.tsv"))
  expect_error(load_params(r), "does not start with the line param<TAB>value")
  r <- fake_run(); unlink(file.path(r, "params.tsv"))
  expect_error(load_params(r), "no params.tsv in")
})

test_that("a params.tsv of another version warns once, naming what it lacks and what it leaves out", {
  r <- fake_run(drop = "MDA_COVERAGE", extra = "MAX_BROOD_AGE\t180")
  w <- list()
  p <- withCallingHandlers(load_params(r),
                           warning = function(x) { w[[length(w) + 1]] <<- x; invokeRestart("muffleWarning") })
  expect_length(w, 1L)
  expect_s3_class(w[[1]], "malariafwd_params_mismatch")
  expect_s3_class(w[[1]], "malariafwd_run_warning")
  expect_match(conditionMessage(w[[1]]), "lacks MDA_COVERAGE \\(its default is used\\)")
  expect_match(conditionMessage(w[[1]]), "MAX_BROOD_AGE is not a parameter \\(left out\\)")
  expect_identical(p$MDA_COVERAGE, vivax_params()$MDA_COVERAGE)
  expect_null(p$MAX_BROOD_AGE)
  ## once per directory per session, whichever reader asks
  expect_no_warning(load_params(r))
  expect_no_warning(load_run(r))
  ## and never from a run object, or from the package's own reads
  r2 <- fake_run(drop = "MDA_COVERAGE")
  expect_no_warning(malariafwd:::.params_of(r2))
  run <- suppressWarnings(load_run(r2))
  expect_no_warning(load_params(run))
  ## the founder counts a file lacks follow one another, as in vivax_params()
  p <- suppressWarnings(load_params(fake_run(params = list(N_FOUNDERS = 7),
                                             drop = c("N_FOUNDER_HOSTS", "N_FOUNDER_LINEAGES"))))
  expect_identical(c(p$N_FOUNDER_HOSTS, p$N_FOUNDER_LINEAGES), c(7L, 7L))
})

test_that("load_run() carries the record rows, typed, and refuses a directory that is not a run", {
  r <- load_run(fake_run(record = list(ACTUAL_SEED = 3e9, SEED = 4L, SLIM_VERSION = "5.2",
                                       MALARIAFWD_VERSION = "unknown", GENOME_LENGTH = 1021664L,
                                       N_CHROMOSOMES = 1L, RECOMB_WITHIN = 7.4e-07, EXTINCT_TICK = 40L)))
  expect_identical(r$seed, 3e9)
  expect_identical(r$extinct_tick, 40L)
  expect_false(r$never_started)
  expect_identical(r$malariafwd_version, "unknown")
  expect_identical(c(r$genome_length, r$n_chromosomes), c(1021664L, 1L))
  expect_identical(r$recomb_within, 7.4e-07)
  expect_null(load_params(r)$EXTINCT_TICK)
  ## SEED when there is no ACTUAL_SEED, and NA for what the file does not record
  r <- load_run(fake_run(record = list(SEED = 4L, NEVER_STARTED = 10L)))
  expect_identical(r$seed, 4L)
  expect_true(r$never_started)
  expect_identical(r$slim_version, NA_character_)
  expect_identical(r$genome_length, NA_integer_)
  d <- tempfile("mfwd_empty_"); dir.create(d)
  expect_error(load_run(d), "is not a malariafwd run: it has no params.tsv or events.tsv")
  r <- fake_run(); unlink(file.path(r, "events.tsv"))
  expect_error(load_run(r), "is not a malariafwd run: it has no events.tsv")
})

test_that("load_state() gives the host rows unless asked for the mosquitoes, or both", {
  r <- fake_run(st_rows = c("1\thuman\t0\t1\t1\t0\t0\t0\t7\tNA\tNA\t0", "1\tmosquito\t2\t1\t1\t0\t0\t0\t7\t0\t3\tNA",
                            "2\thuman\t0\t1\t1\t0\t0\t0\t7\tNA\tNA\t1", "2\tmosquito\t2\t1\t1\t0\t0\t0\t7\t1\t4\tNA"))
  expect_identical(load_state(r)$entity_type, c("human", "human"))
  expect_identical(load_state(r)$infectious_bites, c(0L, 1L))
  expect_identical(load_state(r, entity = "mosquito")$mosq_age_days, c(3L, 4L))
  expect_identical(load_state(r, entity = "all")$entity_type, c("human", "mosquito", "human", "mosquito"))
  for (e in c("human", "mosquito", "all")) {
    expect_identical(classes(load_state(r, entity = e)), st_classes)
    expect_identical(classes(load_state(fake_run(), entity = e)), st_classes)
  }
  expect_error(load_state(r, entity = "hosts"), "should be one of")
})
