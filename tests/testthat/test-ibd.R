## ibd_pairs()'s vivax_ibd object, keyed by pedigree id, and its readers.

toy_ibd <- function() malariafwd:::.new_vivax_ibd(
  pairs = data.table::data.table(a = c(1L, 1L, 2L), b = c(2L, 3L, 3L), ibd = c(0.5, 0, 0.25)),
  samples = data.table::data.table(pedigree_id = 1:3, node_time = c(10, 5, 0), tick = c(90L, 95L, 100L)),
  missing = integer(), segments = NULL, sequence_length = 1e6, min_span = 1000, all_pairs = TRUE,
  end_tick = 100L, trees = "toy.trees")

## co-transmission and outcrossing, so the run holds every kind of pair; made once a session
.ibd_runs <- new.env()
ibd_fixture <- function() {
  if (is.null(.ibd_runs$run))
    .ibd_runs$run <- suppressMessages(run_sim(
      vivax_params(N_HOSTS = 50, N_MOSQ = 100, BITING_RATE = 0.12, SIM_DAYS = 250, N_FOUNDER_LINEAGES = 8,
                   LAMBDA_OOCYSTS = 2, LAMBDA_SPOROZOITES = 2), seed = 3, outdir = tempfile("mfwd_ibd_")))
  .ibd_runs$run
}

test_that("a vivax_ibd object reads back by pedigree id", {
  x <- toy_ibd()
  expect_identical(ibd_lookup(x, 1, 2), 0.5)
  expect_identical(ibd_lookup(x, 2, 1), 0.5)
  expect_identical(ibd_lookup(x, c(1, 2, 3), c(1, 3, 9)), c(1, 0.25, NA))
  expect_identical(ibd_lookup(x, 1, c(2, 3)), c(0.5, 0))
  expect_null(names(ibd_lookup(x, 1, 2)))
  expect_error(ibd_lookup(x, 1:2, 1:3), "same length")
  m <- as.matrix(x)
  expect_identical(dimnames(m), list(c("1", "2", "3"), c("1", "2", "3")))
  expect_identical(m, t(m))
  expect_identical(unname(diag(m)), c(1, 1, 1))
  expect_identical(m["1", "3"], 0)
  expect_identical(names(as.data.frame(x)), c("a", "b", "ibd"))
  out <- capture.output(print(x))
  expect_lte(length(out), 6L)
  expect_match(out[1], "^<vivax_ibd>")
  expect_false(any(grepl("function(", out, fixed = TRUE)))
})

test_that("ibd_pairs() returns a pedigree-keyed vivax_ibd", {
  skip_if_no_slim(); skip_if_no_tskit()
  r <- ibd_fixture(); ev <- load_events(r)
  ids <- classify_cases(ev)$pedigree_id
  x <- ibd_pairs(r, ids = ids)
  expect_s3_class(x, "vivax_ibd")
  expect_identical(vapply(x$pairs, function(v) class(v)[1], ""), c(a = "integer", b = "integer", ibd = "numeric"))
  expect_identical(names(x$samples), c("pedigree_id", "node_time", "tick"))
  n <- nrow(x$samples)
  expect_identical(nrow(x$pairs), as.integer(choose(n, 2)))
  expect_true(all(x$pairs$a < x$pairs$b))
  expect_true(all(x$pairs$ibd >= 0 & x$pairs$ibd <= 1))
  expect_identical(x$end_tick, max(load_state(r)$tick))
  expect_identical(x$samples$tick, as.integer(x$end_tick - x$samples$node_time))
  first <- ev[pedigree_id %in% x$samples$pedigree_id, .(t = min(tick)), by = pedigree_id][order(pedigree_id)]
  expect_identical(x$samples$tick, first$t)
})

test_that("ibd_lookup() gives the IBD a genome pair has", {
  skip_if_no_slim(); skip_if_no_tskit()
  r <- ibd_fixture(); ev <- load_events(r)
  rel <- ev[event_type == "hypnozoite_activation" & parent1_pedigree_id != -1][seq_len(min(.N, 30))]
  skip_if(nrow(rel) < 5, "too few relapses")
  x <- ibd_pairs(r, ids = c(rel$pedigree_id, rel$parent1_pedigree_id))
  expect_true(all(ibd_lookup(x, rel$pedigree_id, rel$parent1_pedigree_id) > 0.999))
  m <- as.matrix(x)
  expect_identical(m[cbind(as.character(rel$pedigree_id), as.character(rel$parent1_pedigree_id))],
                   ibd_lookup(x, rel$pedigree_id, rel$parent1_pedigree_id))
  expect_true(all(ibd_relapse_origin(ev, ibd_pairs(r, ids = classify_cases(ev)$pedigree_id))$ibd > 0.999))
})

test_that("segments carry their common ancestor and its age", {
  skip_if_no_slim(); skip_if_no_tskit()
  r <- ibd_fixture(); ev <- load_events(r)
  rel <- ev[event_type == "hypnozoite_activation" & parent1_pedigree_id != -1][seq_len(min(.N, 10))]
  skip_if(nrow(rel) < 2, "too few relapses")
  ids <- unique(c(rel$pedigree_id, rel$parent1_pedigree_id, utils::head(classify_cases(ev)$pedigree_id, 20)))
  x <- ibd_pairs(r, ids = ids, segments = TRUE)
  sg <- x$segments
  tt <- stats::setNames(x$samples$node_time, x$samples$pedigree_id)
  expect_equal(sg$depth, sg$mrca_time - pmax(tt[as.character(sg$a)], tt[as.character(sg$b)]), ignore_attr = TRUE)
  expect_true(all(sg$depth >= 0))
  expect_true(all(sg$mrca_id %in% ev$pedigree_id))
  pr <- data.table::data.table(a = pmin(rel$pedigree_id, rel$parent1_pedigree_id),
                               b = pmax(rel$pedigree_id, rel$parent1_pedigree_id), origin = rel$parent1_pedigree_id)
  own <- merge(sg, pr, by = c("a", "b"))
  expect_gt(nrow(own), 0)
  expect_true(all(own$mrca_id == own$origin & own$depth == 0))
})

test_that("ids not in the tree sequence are listed, look up NA, and are named in the errors", {
  skip_if_no_slim(); skip_if_no_tskit()
  r <- ibd_fixture()
  ids <- utils::head(classify_cases(load_events(r))$pedigree_id, 5)
  expect_warning(x <- ibd_pairs(r, ids = c(ids, 99999999L)), "1 of the requested ids are not in the tree sequence")
  expect_identical(x$missing, 99999999L)
  expect_identical(ibd_lookup(x, ids[1], 99999999L), NA_real_)
  expect_error(ibd_pairs(r, ids = 99999999L), "RETAIN_FULL_HISTORY")
  expect_error(ibd_pairs(r, ids = c(ids, NA)), "`ids` must be pedigree ids")
  expect_error(ibd_pairs(r, ids = ids, min_span = -1), "`min_span` must be")
  expect_error(ibd_pairs(r, ids = ids, min_span = NA), "`min_span` must be")
  expect_error(ibd_pairs(r, ids = ids, segments = NA), "`segments` must be TRUE or FALSE")
})

test_that("a run with no parasite gives an empty object", {
  skip_if_no_slim(); skip_if_no_tskit()
  r <- suppressMessages(run_sim(vivax_params(N_HOSTS = 20, N_MOSQ = 20, SIM_DAYS = 20, N_FOUNDER_HOSTS = 0,
                                             IMPORT_START = 1, IMPORT_END = 10, IMPORT_RATE_PER_TICK = 1e-6),
                                seed = 1, outdir = tempfile("mfwd_ibd_none_")))
  expect_true(r$never_started)
  x <- ibd_pairs(r)
  expect_identical(c(nrow(x$samples), nrow(x$pairs)), c(0L, 0L))
  expect_identical(names(x$pairs), c("a", "b", "ibd"))
  expect_output(print(x), "0 genomes, 0 pairs")
})
