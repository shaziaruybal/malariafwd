## Episodes and relapses: the masked relapses, into a patent clone of their own lineage.

test_that("classify_cases() marks the masked relapse, and leaves it out when asked", {
  fx <- episodes_fixture()
  cs <- classify_cases(fx)
  expect_identical(nrow(cs), 12L)
  expect_identical(cs$pedigree_id[cs$into_patent_lineage], 104L)
  expect_false(anyNA(cs$into_patent_lineage))
  cm <- classify_cases(fx, exclude_masked = TRUE)
  expect_identical(nrow(cm), 11L)
  expect_false(104L %in% cm$pedigree_id)
})

test_that("force_of_infection() leaves the masked relapse out when asked", {
  r <- episodes_fixture_run()
  all <- force_of_infection(r); unm <- force_of_infection(r, exclude_masked = TRUE)
  rel <- function(fo) sum(fo[grepl("^relapse", origin)]$n)
  expect_identical(c(rel(all), rel(unm)), c(4L, 3L))
  expect_identical(all[!grepl("^relapse", origin)], unm[!grepl("^relapse", origin)])
})

test_that("exclude_masked must be TRUE or FALSE", {
  r <- episodes_fixture_run(); fx <- episodes_fixture()
  for (v in list(NA, "yes", c(TRUE, FALSE))) {
    expect_error(classify_cases(fx, exclude_masked = v), "`exclude_masked` must be TRUE or FALSE")
    expect_error(relapse_fraction(fx, exclude_masked = v), "`exclude_masked` must be TRUE or FALSE")
    expect_error(force_of_infection(r, exclude_masked = v), "`exclude_masked` must be TRUE or FALSE")
  }
})

test_that("a run without the flag gives NA on its relapses, and cannot leave them out", {
  fx <- episodes_fixture()
  for (ev in list(data.table::copy(fx)[event_type == "hypnozoite_activation", into_patent_lineage := NA_integer_],
                  data.table::copy(fx)[, into_patent_lineage := NULL])) {
    cs <- classify_cases(ev)
    expect_true(all(is.na(cs[origin == "relapse"]$into_patent_lineage)))
    expect_false(anyNA(cs[origin != "relapse"]$into_patent_lineage))
    expect_error(classify_cases(ev, exclude_masked = TRUE), "4 of these events' activation rows lack it")
  }
})

test_that("on a run, exclude_masked removes exactly the flagged activations", {
  skip_if_no_slim()
  r <- suppressMessages(run_sim(vivax_params(N_HOSTS = 60, N_MOSQ = 150, SIM_DAYS = 400, CHROMOSOMES = "1",
                                             BLOODCURE_COVERAGE = 0.5, P_BLOODCURE_SURVIVAL = 0.3,
                                             P_RECRUD_GIVEN_FAIL = 0.5, MDA_TICK = 200, MDA_COVERAGE = 0.5),
                                seed = 3, outdir = tempfile("mfwd_masked_")))
  ev <- load_events(r)
  flagged <- ev[event_type == "hypnozoite_activation" & into_patent_lineage == 1L, .N]
  expect_gt(flagged, 0L)
  expect_identical(nrow(classify_cases(ev)) - nrow(classify_cases(ev, exclude_masked = TRUE)), flagged)
})

ep_classes <- function(d) vapply(d, function(x) class(x)[1], "")

test_that("classify_cases() is one row per clone onset, with its episode", {
  cs <- classify_cases(episodes_fixture())
  expect_identical(names(cs), c("tick", "human_id", "episode_id", "lineage_id", "oocyst_id", "pedigree_id",
                                "origin", "imported", "detected", "into_patent_lineage"))
  expect_identical(nrow(cs), 12L)
  expect_identical(cs[pedigree_id %in% c(101L, 102L)]$episode_id, c("1:10", "1:10"))
  expect_identical(cs[pedigree_id %in% c(403L, 404L)]$episode_id, c("4:45", "4:45"))
})

test_that("episodes() gives every blood-stage clone's span, founders included", {
  ep <- episodes(episodes_fixture())
  expect_identical(nrow(ep), 13L)
  expect_identical(anyDuplicated(ep$pedigree_id), 0L)
  f <- ep[pedigree_id == 201L]
  expect_identical(c(f$episode_id, f$onset_event, f$end_reason), c("2:1", "vh_inoculation_seed", "human_clearance"))
  expect_identical(f$end_tick, 15L)
  expect_identical(ep[pedigree_id == 301L, c(end_tick)], 12L)
  expect_identical(ep[pedigree_id == 301L, end_reason], "bloodcure_subpatent")
  og <- ep[end_reason == "ongoing"]
  expect_identical(sort(og$pedigree_id), c(105L, 404L, 503L))
  expect_true(all(is.na(og$end_tick)))
})

test_that("recurrences() pairs each recurrence with the earlier episodes of its host", {
  fx <- episodes_fixture()
  r <- recurrences(fx)
  expect_identical(nrow(r), 9L)
  expect_identical(nrow(recurrences(fx, prior = "previous")), 7L)
  expect_identical(nrow(recurrences(fx, from = 40)), 6L)
  expect_identical(nrow(recurrences(fx, exclude_masked = TRUE)), 6L)
  expect_identical(nrow(recurrences(fx, prior = "previous", exclude_masked = TRUE)), 6L)
  expect_false(any(r$prior_episode_id == r$episode_id))
  expect_true(all(r$days_since >= 1L))
  expect_identical(r[pedigree_id %in% c(403L, 404L)]$prior_pedigree_id, c(401L, 401L))
  expect_identical(recurrences(fx, prior = "previous")[pedigree_id == 105L]$prior_pedigree_id, 104L)
  expect_identical(recurrences(fx, prior = "previous", exclude_masked = TRUE)[pedigree_id == 105L]$prior_pedigree_id,
                   c(101L, 102L))
  expect_false(2L %in% r$human_id)
  rc <- r[pedigree_id == 302L]
  expect_identical(c(rc$prior_pedigree_id, rc$days_since), c(301L, 17L))
  expect_identical(rc$prior_origin, "local reinfection")
  im <- r[pedigree_id == 503L]
  expect_identical(c(im$origin, im$prior_origin), c("imported", "relapse"))
  expect_identical(c(im$prior_pedigree_id, im$days_since), c(502L, 30L))
  expect_true(im$imported)
})

test_that("recurrences() reads IBD from the run's object, NA for a genome it does not hold", {
  fx <- episodes_fixture()
  ids <- c(101L, 102L, 104L, 105L)
  pr <- t(utils::combn(ids, 2))
  stub <- malariafwd:::.new_vivax_ibd(
    pairs = data.table::data.table(a = pr[, 1], b = pr[, 2], ibd = c(0.1, 0.9, 0.95, 0.2, 0.3, 1)),
    samples = data.table::data.table(pedigree_id = ids, node_time = c(90, 90, 80, 50), tick = c(10L, 10L, 20L, 50L)),
    missing = integer(), segments = NULL, sequence_length = 1, min_span = 0, all_pairs = TRUE, end_tick = 100L,
    trees = "stub")
  r <- recurrences(fx, ibd = stub)
  h1 <- r[human_id == 1L]
  expect_identical(h1$ibd, ibd_lookup(stub, h1$pedigree_id, h1$prior_pedigree_id))
  expect_true(all(is.na(r[human_id != 1L]$ibd)))
})

test_that("the readers take a table, a run or its directory, and the same tables come back", {
  fx <- episodes_fixture(); run <- episodes_fixture_run()
  for (f in list(classify_cases, episodes, recurrences)) {
    expect_identical(f(fx), f(run))
    expect_identical(f(run), f(run$outdir))
  }
  expect_error(episodes(42), "`events` must be an events table")
})

test_that("a run with no parasite gives typed empty tables", {
  head_only <- data.table::fread(text = paste(ev_cols, collapse = "\t"), sep = "\t")
  cs <- classify_cases(episodes_fixture())
  expect_no_warning(z <- classify_cases(head_only))
  expect_identical(nrow(z), 0L); expect_identical(ep_classes(z), ep_classes(cs))
  ez <- episodes(head_only)
  expect_identical(nrow(ez), 0L); expect_identical(ep_classes(ez), ep_classes(episodes(episodes_fixture())))
  rz <- recurrences(head_only)
  expect_identical(nrow(rz), 0L); expect_identical(ep_classes(rz), ep_classes(recurrences(episodes_fixture())))
})

test_that("recurrences() refuses arguments it cannot read, naming them", {
  fx <- episodes_fixture()
  expect_error(recurrences(fx, prior = "last"), "should be one of")
  for (v in list(-1, 1.5, NA, c(1, 2))) expect_error(recurrences(fx, from = v), "`from` must be a single whole number")
  expect_error(recurrences(fx, ibd = 1), "`ibd` must be the run's IBD object")
})

test_that("plot_host_trajectory() draws each clone from episodes(), to the last day of state", {
  skip_if_not_installed("ggplot2")
  st <- data.table::data.table(tick = 1:100, entity_type = "human", entity_id = 1L, tag = 0L, moi_active = 0L,
                               moi_subpatent = 0L, moi_dormant = 0L, hyp_reservoir_size = 0L, lineage_ids = NA_character_,
                               n_transmissible = NA_integer_, mosq_age_days = NA_integer_, infectious_bites = 0L)
  p <- plot_host_trajectory(st, episodes_fixture(), 1)
  area <- p$layers[[1]]$data
  expect_identical(max(area$tick), 100L)
  d25 <- area[tick == 25L]
  expect_identical(d25[lineage_id == "1"]$n, 2L)
  expect_identical(d25[lineage_id == "2"]$n, 1L)
})

test_that("with co-transmission, an episode holds several clones and is never paired with itself", {
  skip_if_no_slim()
  r <- suppressMessages(run_sim(vivax_params(N_HOSTS = 100, N_MOSQ = 200, BITING_RATE = 0.11, SIM_DAYS = 548,
                                             LAMBDA_SPOROZOITES = 2, CHROMOSOMES = "1"),
                                seed = 1, outdir = tempfile("mfwd_cotx_")))
  ev <- load_events(r)
  cs <- classify_cases(ev)
  expect_gt(max(cs[, .N, by = episode_id]$N), 1L)
  rec <- recurrences(ev)
  expect_true(all(rec$days_since >= 1L))
  ## each clone against every earlier row of its host, in the order of the log, less the same-day pairs
  rows <- cs[, .(pedigree_id, tick, human_id, i = seq_len(.N))]
  pr <- rows[rows, on = .(human_id, i < i), nomatch = NULL, allow.cartesian = TRUE,
             .(pedigree_id = i.pedigree_id, prior_pedigree_id = x.pedigree_id, same_day = i.tick == x.tick)]
  expect_gt(sum(pr$same_day), 0L)
  expect_setequal(paste(rec$pedigree_id, rec$prior_pedigree_id),
                  paste(pr[!(same_day)]$pedigree_id, pr[!(same_day)]$prior_pedigree_id))
})

test_that("a clone's span covers the days it counts in moi_active", {
  skip_if_no_slim()
  r <- suppressMessages(run_sim(vivax_params(N_HOSTS = 60, N_MOSQ = 150, SIM_DAYS = 400, CHROMOSOMES = "1",
                                             BLOODCURE_COVERAGE = 0.5, P_BLOODCURE_SURVIVAL = 0.3,
                                             P_RECRUD_GIVEN_FAIL = 0.5, MDA_TICK = 200, MDA_COVERAGE = 0.5),
                                seed = 3, outdir = tempfile("mfwd_spans_")))
  ep <- episodes(r)
  st <- load_state(r)
  last <- max(st$tick)
  span <- ep[, .(tick = if (is.na(end_tick)) seq(onset_tick, last) else if (end_tick > onset_tick)
                   seq(onset_tick, end_tick - 1L) else integer()), by = .(human_id, pedigree_id, lineage_id)]
  have <- span[, .(n = data.table::uniqueN(lineage_id)), by = .(human_id, tick)]
  m <- merge(st[, .(human_id = entity_id, tick, moi_active)], have, by = c("human_id", "tick"), all.x = TRUE)
  m[is.na(n), n := 0L]
  expect_identical(m$n, m$moi_active)
})
