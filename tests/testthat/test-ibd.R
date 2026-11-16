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
  expect_error(ibd_pairs(r, ids = 99999999L), "none of the 1 requested ids is in the tree sequence.", fixed = TRUE)
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

## ---- segments (D18) ----

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

test_that("new_infection_diversity()'s clonal_frac is the share of pairs at or above `clonal`", {
  skip_if_no_slim(); skip_if_no_tskit()
  r <- ibd_fixture()
  d <- new_infection_diversity(r, cap = 20)
  ## the same sample of the acquisitions after the burn-in (days 126 to 250), and its pairs
  pool <- classify_cases(load_events(r))[origin %in% c("imported", "local reinfection") & tick >= 126L,
                                         unique(pedigree_id)]
  x <- suppressWarnings(ibd_pairs(r, ids = malariafwd:::.sample_fixed(pool, 20, 1)))
  share <- mean(x$pairs$ibd >= 0.5)
  expect_gt(share, 0)
  expect_identical(d$n, length(x$samples$pedigree_id))
  expect_equal(d$clonal_frac, share)
  expect_identical(new_infection_diversity(r, cap = 20, clonal = 1.01)$clonal_frac, 0)
})

test_that("new_infection_diversity() gives one row for the window, or one for every window with bin", {
  skip_if_no_slim(); skip_if_no_tskit()
  r <- ibd_fixture()   # SIM_DAYS 250: the default window is days 126 to 250
  one <- new_infection_diversity(r, cap = 20)
  expect_identical(vapply(one, function(v) class(v)[1], ""),
                   c(win = "integer", days = "integer", n = "integer", n_pairs = "integer", mean_ibd = "numeric",
                     clonal_frac = "numeric", n_lineages = "integer"))
  expect_identical(c(one$win, one$days), c(126L, 125L))
  expect_identical(one$n_pairs, as.integer(choose(one$n, 2)))
  expect_identical(new_infection_diversity(r, cap = 20, bin = 251)[, !"win"], one[, !"win"])
  by60 <- new_infection_diversity(r, cap = 20, bin = 60)
  expect_identical(by60$win, c(120L, 180L, 240L))
  expect_identical(by60$days, c(54L, 60L, 11L))
  ## every genome is in the tree sequence, so each window compares its acquisitions, at most 20
  cs <- classify_cases(load_events(r))[origin %in% c("imported", "local reinfection") & tick >= 126L]
  expect_identical(length(suppressWarnings(ibd_pairs(r, ids = cs$pedigree_id))$missing), 0L)
  per <- cs[, .N, by = .(win = as.integer((tick %/% 60) * 60))]
  expect_identical(by60$n, pmin(data.table::fcoalesce(per$N[match(by60$win, per$win)], 0L), 20L))
  expect_identical(is.na(by60$mean_ibd), by60$n_pairs < 15L)
  nv <- new_infection_diversity(never_run(), bin = 10)
  expect_identical(nrow(nv), 4L)
  expect_true(all(nv$n == 0L) && all(is.na(nv$mean_ibd)))
  expect_error(new_infection_diversity(r, min_pairs = 0), "`min_pairs` must be a whole number, at least 1.",
               fixed = TRUE)
  expect_error(new_infection_diversity(r, bin = 0), "`bin` must be a whole number of days, at least 1")
  expect_error(new_infection_diversity(r, burn_frac = 0.2, from = 10), "give `burn_frac` or `from`, not both.",
               fixed = TRUE)
})

test_that("new_infection_diversity() compares only the pairs within a window, each as ibd_pairs() reads it", {
  skip_if_no_slim(); skip_if_no_tskit()
  r <- ibd_fixture()
  cs <- classify_cases(load_events(r))[origin %in% c("imported", "local reinfection") & tick >= 126L]
  g <- split(cs$pedigree_id, (cs$tick %/% 60) * 60)
  x <- malariafwd:::.ibd_compute(r, groups = g)
  ## every pair it holds is of one window
  grp <- rep(seq_along(g), lengths(g)); names(grp) <- unlist(g)
  expect_gt(nrow(x$got), 0L)
  expect_identical(unname(grp[as.character(x$got$a)]), unname(grp[as.character(x$got$b)]))
  by60 <- new_infection_diversity(r, cap = 1e6, bin = 60)
  for (k in seq_along(g)) {
    row <- by60[win == as.integer(names(g)[k])]
    if (row$n_pairs >= 15L) expect_equal(row$mean_ibd, mean(ibd_pairs(r, ids = g[[k]])$pairs$ibd), label = names(g)[k])
    else expect_true(is.na(row$mean_ibd))
  }
})

## a relapse pair twice under one label, a sibling pair, and two founders: four pairs, three groups
grouped_pairs <- function(run) {
  tp <- ibd_truth_pairs(run)
  rbind(tp[cls == "relapse vs its origin"][1:2], tp[cls == "same-oocyst sibs"][1],
        tp[cls == "distinct founders"][1])
}

test_that("plot_relatedness_ladder(style = \"bars\") keys every band's colour, used or not", {
  skip_if_not_installed("ggplot2", "3.5.0")
  ## two pairs, so four of the six bands hold no pair
  tp <- data.table::data.table(cls = factor(c("relapse vs its origin", "distinct founders"),
                                            levels = malariafwd:::.LADDER), ibd = c(1, 0))
  p <- plot_relatedness_ladder(tp, style = "bars")
  expect_true(isTRUE(p$layers[[1]]$show.legend))
  expect_identical(nrow(ggplot2::get_guide_data(p, "fill")), 6L)
})

test_that("plot_ibd_pair() refuses labels that are not one per pair", {
  skip_if_not_installed("ggplot2")
  expect_error(plot_ibd_pair(NULL, 1:3, 4:6, labels = c("x", "y")),
               "`labels` must give one label per pair.*there are 3 pairs, and it has 2 values")
  expect_error(plot_ibd_pair(NULL, 1:2, 3:4, labels = c("x", NA)), "some NA", fixed = TRUE)
})

test_that("plot_ibd_pair() shows NA, not 0, for a pair with a genome the tree sequence lacks", {
  skip_if_not_installed("ggplot2"); skip_if_no_tskit()
  run <- example_run()
  tp <- ibd_truth_pairs(run)[cls == "relapse vs its origin"][1]
  ids_a <- c(tp$a, tp$a); ids_b <- c(tp$b, 999999999)
  expect_warning(p <- plot_ibd_pair(run, ids_a, ids_b, labels = c("real", "absent")),
                 "not in the tree sequence")
  txt <- ggplot2::layer_data(p, 3)$label
  expect_identical(sort(txt), sort(c("IBD = 1.00", "IBD = NA")))
  expect_warning(pg <- plot_ibd_pair(run, ids_a, ids_b, labels = c("real", "absent"), group = c("g", "g")),
                 "not in the tree sequence")
  expect_identical(pg$scales$get_scales("y")$labels(c("1", "2")), c("real (IBD 1.00)", "absent (IBD NA)"))
})

test_that("plot_ibd_pair(group =) draws a panel per group, a row per pair and a colour per group", {
  skip_if_no_tskit(); skip_if_not_installed("ggplot2")
  run <- example_run(); pk <- grouped_pairs(run)
  grp <- c("relapse", "relapse", "siblings", "founders")
  lab <- c("r", "r", "s", "f")   # a label given twice in a group still gets two rows
  seg <- ibd_segments(run, pk$a, pk$b)
  for (layout in c("chromosome", "concatenated")) {
    p <- plot_ibd_pair(run, pk$a, pk$b, labels = lab, group = grp, layout = layout)
    b <- ggplot2::ggplot_build(p)
    lay <- data.table::as.data.table(b$layout$layout)
    expect_identical(levels(lay$grp), c("relapse", "siblings", "founders"), label = layout)
    expect_identical(nrow(lay), if (layout == "chromosome") 3L * 14L else 3L, label = layout)
    ## rows top to bottom in the order given, two rows for the relapse pairs
    first <- lay[, .SD[1], by = grp]
    rows <- lapply(first$PANEL, function(k) b$layout$panel_params[[k]]$y$get_labels())
    if (layout == "chromosome") rows <- lapply(rows, sub, pattern = " \\(IBD .*", replacement = "")
    expect_identical(lapply(rows, rev), list(c("r", "r"), "s", "f"), label = layout)
    ## every segment drawn, each in its group's colour
    sg <- data.table::as.data.table(ggplot2::layer_data(p, 2))
    expect_identical(nrow(sg), nrow(seg), label = layout)
    sg <- merge(sg, lay[, .(PANEL, grp)], by = "PANEL")
    expect_identical(sort(unique(sg[, paste(grp, fill)])),
                     sort(paste(unique(grp[seg$pair]), c(relapse = "#2166ac", siblings = "#c8971f")[unique(grp[seg$pair])])),
                     label = layout)
  }
  ## a factor's levels set the order, unused ones dropped; colours one per group, or named
  f <- factor(grp, c("unused", "founders", "siblings", "relapse"))
  p <- plot_ibd_pair(run, pk$a, pk$b, group = f, fill = c(relapse = "black", siblings = "red", founders = "grey"),
                     layout = "concatenated")
  lay <- data.table::as.data.table(ggplot2::ggplot_build(p)$layout$layout)
  expect_identical(levels(lay$grp), c("founders", "siblings", "relapse"))
  sg <- merge(data.table::as.data.table(ggplot2::layer_data(p, 2)), lay[, .(PANEL, grp)], by = "PANEL")
  expect_identical(sort(unique(sg[, paste(grp, fill)])), c("relapse black", "siblings red"))
  p <- plot_ibd_pair(run, pk$a, pk$b, group = grp, fill = "black", layout = "concatenated")
  expect_identical(unique(ggplot2::layer_data(p, 2)$fill), "black")
})

test_that("plot_ibd_pair() refuses a group or a fill that does not fit the pairs", {
  skip_if_no_tskit(); skip_if_not_installed("ggplot2")
  run <- example_run(); pk <- grouped_pairs(run)
  grp <- c("relapse", "relapse", "siblings", "founders")
  expect_error(plot_ibd_pair(run, pk$a, pk$b, group = grp[1:3]), "one group per pair")
  expect_error(plot_ibd_pair(run, pk$a, pk$b, group = replace(grp, 2, NA)), "one group per pair")
  expect_error(plot_ibd_pair(run, pk$a, pk$b, group = grp, fill = c("red", "blue")), "one for each of the 3 groups")
  expect_error(plot_ibd_pair(run, pk$a, pk$b, group = grp, fill = c(relapse = "red", siblings = "blue")),
               "no colour for group \"founders\"")
})

test_that("plot_ibd_pair() without group keeps a panel per label", {
  skip_if_no_tskit(); skip_if_not_installed("ggplot2")
  run <- example_run(); pk <- grouped_pairs(run)
  for (layout in c("chromosome", "concatenated")) {
    lay <- ggplot2::ggplot_build(plot_ibd_pair(run, pk$a, pk$b, labels = c("r", "r", "s", "f"), layout = layout))$layout$layout
    expect_identical(as.character(lay$lab), c("r", "s", "f"), label = layout)
  }
})
