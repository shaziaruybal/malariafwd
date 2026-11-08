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

## ---- segments (D74) ----

test_that("segments are maximal stretches from one ancestor, cut at the junctions, with min_span after the cut", {
  tr <- malariafwd:::.ibd_tracts
  lay <- data.table::data.table(chrom = c(1L, 2L), length = c(1000, 2000), start = c(0, 1000), end = c(1000, 3000))
  raw <- function(...) data.table::rbindlist(list(...))
  row <- function(left, right, mrca, a = 1L, b = 2L) data.table::data.table(a = a, b = b, left = left, right = right,
                                                                          mrca_id = as.integer(mrca), mrca_time = 10)
  ## abutting rows from one ancestor merge; two ancestors do not; a gap separates
  x <- tr(raw(row(0, 200, 7), row(200, 500, 7), row(500, 700, 8), row(800, 900, 8)), lay, 0)
  expect_identical(x$gstart, c(0, 500, 800)); expect_identical(x$gend, c(500, 700, 900))
  expect_identical(x$mrca_id, c(7L, 8L, 8L))
  ## a tract over the junction gives two pieces, in each chromosome's coordinates
  y <- tr(raw(row(900, 1100, 7)), lay, 0)
  expect_identical(y$chrom, 1:2); expect_identical(c(y$start, y$end), c(900, 0, 1000, 100))
  expect_identical(y$length, c(100, 100))
  ## min_span is applied after the cut, inclusively
  expect_identical(nrow(tr(raw(row(900, 1100, 7)), lay, 100)), 2L)
  expect_identical(nrow(tr(raw(row(900, 1100, 7)), lay, 150)), 0L)
  ## a tract that ends at the junction is not split into an empty piece
  expect_identical(tr(raw(row(500, 1000, 7)), lay, 0)$chrom, 1L)
  e <- tr(row(0, 1, 1)[0], lay, 0)
  expect_identical(names(e), c("a", "b", "chrom", "start", "end", "gstart", "gend", "length", "mrca_id", "mrca_time"))
  expect_identical(nrow(e), 0L)
  expect_identical(class(e$chrom), "integer")
})

test_that("a pair's segments are the same whatever other genomes are passed, from either function", {
  skip_if_no_slim(); skip_if_no_tskit()
  r <- ibd_fixture(); ids <- classify_cases(load_events(r))$pedigree_id
  x <- ibd_pairs(r, ids = ids, segments = TRUE)
  sub <- utils::head(ids, 40)
  y <- ibd_pairs(r, ids = sub, segments = TRUE)
  expect_equal(y$segments, x$segments[a %in% sub & b %in% sub], ignore_attr = TRUE)
  expect_equal(y$pairs$ibd, x$pairs[y$pairs, on = .(a, b)]$ibd)
  pp <- x$pairs[ibd > 0 & ibd < 1][seq_len(min(.N, 30))]
  s <- ibd_segments(r, pp$a, pp$b)
  ref <- x$segments[pp[, .(a, b, pair = .I)], on = .(a, b)][order(pair, gstart)]
  expect_equal(s[, !"pair"], ref[, names(s)[-1], with = FALSE], ignore_attr = TRUE)
  expect_identical(s$pair, ref$pair)
})

test_that("no segment crosses a chromosome junction or is shorter than min_span", {
  skip_if_no_slim(); skip_if_no_tskit()
  r <- ibd_fixture(); ids <- utils::head(classify_cases(load_events(r))$pedigree_id, 60)
  cm <- chrom_map(r)
  x <- ibd_pairs(r, ids = ids, segments = TRUE)
  pp <- x$pairs[ibd > 0][seq_len(min(.N, 40))]
  for (sg in list(x$segments, ibd_segments(r, pp$a, pp$b))) {
    c0 <- cm$start[match(sg$chrom, cm$chrom)]; c1 <- cm$end[match(sg$chrom, cm$chrom)]
    expect_true(all(sg$gstart >= c0 & sg$gend <= c1))
    expect_identical(sg$start, sg$gstart - c0)
    expect_identical(sg$length, sg$gend - sg$gstart)
    expect_true(all(sg$length >= 1000))
  }
})

test_that("at min_span 0 the segments add up to each pair's IBD", {
  skip_if_no_slim(); skip_if_no_tskit()
  r <- ibd_fixture(); ids <- utils::head(classify_cases(load_events(r))$pedigree_id, 60)
  x <- ibd_pairs(r, ids = ids, segments = TRUE, min_span = 0)
  s <- x$segments[, .(f = sum(length) / x$sequence_length), by = .(a, b)][x$pairs[ibd > 0], on = .(a, b)]
  expect_false(anyNA(s$f))
  expect_lt(max(abs(s$f - s$ibd)), 1e-12)
  y <- ibd_pairs(r, ids = ids, segments = TRUE)
  expect_identical(y$pairs$ibd, x$pairs$ibd)
  t <- y$segments[, .(f = sum(length) / y$sequence_length), by = .(a, b)][y$pairs, on = .(a, b), nomatch = NULL]
  expect_true(all(t$f <= t$ibd + 1e-12))
})

test_that("segment tables keep one schema, and a genome with itself is shared end to end", {
  skip_if_no_slim(); skip_if_no_tskit()
  r <- ibd_fixture(); ev <- load_events(r)
  ids <- classify_cases(ev)$pedigree_id
  rel <- ev[event_type == "hypnozoite_activation" & parent1_pedigree_id != -1][1]
  x <- ibd_pairs(r, ids = ids)
  none <- x$pairs[ibd == 0][1]
  s0 <- ibd_segments(r, none$a, none$b)
  s1 <- ibd_segments(r, rel$pedigree_id, rel$parent1_pedigree_id)
  expect_identical(nrow(s0), 0L)
  cls <- function(d) vapply(d, function(v) class(v)[1], "")
  expect_identical(cls(s0), cls(s1))
  both <- data.table::rbindlist(list(s1, s0, s1))
  expect_true(all(both$chrom %in% 1:14))
  ## a clonal relapse shares every chromosome end to end with its origin, at depth 0
  cm <- chrom_map(r)
  expect_identical(s1$chrom, as.integer(cm$chrom))
  expect_true(all(s1$depth == 0))
  me <- ibd_segments(r, ids[1], ids[1])
  expect_identical(c(nrow(me), unique(me$mrca_id), unique(me$depth)), c(nrow(cm), ids[1], 0))
  expect_identical(me$length, cm$length)
  ## pairs count from 1, a pair given twice has its rows twice
  d <- ibd_segments(r, c(rel$pedigree_id, rel$pedigree_id), c(rel$parent1_pedigree_id, rel$parent1_pedigree_id))
  expect_identical(unique(d$pair), 1:2)
  expect_identical(nrow(d), 2L * nrow(s1))
  expect_error(ibd_segments(r, 1:2, 1), "`ids_a` and `ids_b` must have the same length")
})

test_that("pairs mode computes only the listed pairs, with the values of the ids mode", {
  skip_if_no_slim(); skip_if_no_tskit()
  r <- ibd_fixture(); ids <- classify_cases(load_events(r))$pedigree_id
  x <- ibd_pairs(r, ids = ids, segments = TRUE)
  pp <- x$pairs[unique(round(seq(1, .N, length.out = 30)))]
  ## five listed twice, the other way round, and one genome with itself
  listed <- data.frame(a = c(pp$a, pp$b[1:5], ids[1]), b = c(pp$b, pp$a[1:5], ids[1]))
  y <- ibd_pairs(r, pairs = listed, segments = TRUE)
  expect_false(y$all_pairs)
  expect_identical(nrow(y$pairs), nrow(pp))
  expect_equal(y$pairs$ibd, x$pairs[y$pairs, on = .(a, b)]$ibd)
  expect_equal(y$segments, x$segments[y$pairs[, .(a, b)], on = .(a, b), nomatch = NULL][order(a, b, gstart)],
               ignore_attr = TRUE)
  expect_identical(ibd_lookup(y, ids[1], ids[1]), 1)
  sampled <- y$samples$pedigree_id
  other <- x$pairs[!y$pairs, on = .(a, b)][a %in% sampled & b %in% sampled][1]
  expect_identical(ibd_lookup(y, other$a, other$b), NA_real_)
  m <- as.matrix(y)
  expect_identical(sum(!is.na(m)), 2L * nrow(pp) + nrow(m))
  expect_output(print(y), paste(nrow(pp), "listed pairs"))
  expect_error(ibd_pairs(r, ids = ids, pairs = listed), "give `ids` or `pairs`, not both")
  expect_error(ibd_pairs(r, pairs = data.frame(x = 1)), "`pairs` must be a data frame")
  expect_error(ibd_pairs(r, pairs = data.frame(a = c(1, NA), b = 2:3)), "`pairs$a` must be pedigree ids", fixed = TRUE)
})

test_that("ibd_truth_pairs() computes only its pairs, with the values of every pair", {
  skip_if_no_slim(); skip_if_no_tskit()
  r <- ibd_fixture()
  tp <- ibd_truth_pairs(r)
  z <- ibd_pairs(r, ids = unique(c(tp$a, tp$b)))
  expect_equal(tp$ibd, ibd_lookup(z, tp$a, tp$b))
  expect_gt(nrow(tp), 0)
})
