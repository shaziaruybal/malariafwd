## The run loaders. Most of these write a small run directory by hand, so they need no SLiM.

ev_cols <- c("tick", "event_id", "event_type", "inoculation_id", "mosq_id", "human_id", "lineage_id",
             "parent1_lineage_id", "parent2_lineage_id", "is_recombinant", "p_outcross", "recipient",
             "is_dormant", "is_imported", "pedigree_id", "parent1_pedigree_id", "parent2_pedigree_id",
             "n_hypnozoites", "oocyst_id", "detected", "into_patent_lineage")
st_cols <- c("tick", "entity_type", "entity_id", "tag", "moi_active", "moi_subpatent", "moi_dormant",
             "hyp_reservoir_size", "lineage_ids", "n_transmissible", "mosq_age_days", "infectious_bites")

## a run directory with the given rows (each a tab-separated line) under the given headers
fake_run <- function(ev_rows = character(), st_rows = character(), ev_head = ev_cols, st_head = st_cols) {
  d <- tempfile("mfwd_fake_"); dir.create(d)
  writeLines(c(paste(ev_head, collapse = "\t"), ev_rows), file.path(d, "events.tsv"))
  writeLines(c(paste(st_head, collapse = "\t"), st_rows), file.path(d, "state.tsv"))
  writeLines(c("param\tvalue", "N_HOSTS\t2"), file.path(d, "params.tsv"))
  d
}
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
