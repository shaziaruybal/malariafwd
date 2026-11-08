## transmission_links(): who infected whom, through a bite's blood stage or a relapse of the batch it
## left.

tl_classes <- function(x) vapply(x, function(v) class(v)[1], "")

test_that("transmission_links() gives every link of the fixture, direct and liver", {
  l <- transmission_links(load_events(transmission_fixture_run()))
  want <- data.table::fread(text = "
source_pedigree_id,target_pedigree_id,kind,source_human_id,target_human_id,mosq_id,oocyst_id,source_tick,uptake_tick,inoculation_tick,target_tick,weight
0,2,direct,0,1,100,1,1,5,16,16,1
30,31,direct,7,8,106,8,1,10,22,22,1
2,14,direct,1,1,104,5,16,25,38,38,1
2,15,direct,1,6,105,6,16,26,48,48,1
7,16,direct,3,6,105,7,30,36,48,48,1
0,4,liver,0,1,100,1,1,5,16,60,1
15,6,direct,6,2,101,2,48,70,85,85,0.5
16,6,direct,6,2,101,2,48,70,85,85,0.5
0,5,liver,0,1,100,1,1,5,16,120,1
12,13,direct,2,5,103,4,105,110,122,122,1
15,18,liver,6,2,101,2,48,70,85,130,0.5
16,18,liver,6,2,101,2,48,70,85,130,0.5
7,11,liver,3,4,102,3,30,35,50,140,1", colClasses = c(rep("integer", 2), "character", rep("integer", 8), "numeric"))
  expect_equal(l[, 1:12], want, ignore_attr = TRUE)
  expect_identical(tl_classes(l)[1:12], tl_classes(want))
})

test_that("only local infections are targets, and no host's relapse is linked to its own infection", {
  l <- transmission_links(load_events(transmission_fixture_run()))
  expect_setequal(unique(l$target_pedigree_id), c(2L, 4L, 5L, 6L, 11L, 13L, 14L, 15L, 16L, 18L, 31L))
  ## founders, the import, its batch's relapse, the recrudescence and the batches
  expect_false(any(c(0L, 30L, 7L, 9L, 12L, 3L, 8L, 10L, 17L) %in% l$target_pedigree_id))
  expect_false(any(l$source_pedigree_id %in% c(2L, 14L) & l$target_pedigree_id %in% c(4L, 5L)))
})

test_that("each parent of an outcrossed oocyst carries half a link, and a mosquito can bite where it fed", {
  l <- transmission_links(load_events(transmission_fixture_run()))
  for (t in c(6L, 18L)) {
    x <- l[target_pedigree_id == t]
    expect_identical(x$source_pedigree_id, c(15L, 16L))
    expect_identical(x$weight, c(0.5, 0.5))
    expect_identical(x$source_human_id, c(6L, 6L))
  }
  expect_true(all(l[, sum(weight), by = target_pedigree_id]$V1 == 1))
  self <- l[source_pedigree_id == 2L & target_pedigree_id == 14L]
  expect_identical(c(self$source_human_id, self$target_human_id), c(1L, 1L))
})

test_that("kind keeps the direct or the liver links", {
  ev <- load_events(transmission_fixture_run())
  a <- transmission_links(ev)
  d <- transmission_links(ev, kind = "direct"); v <- transmission_links(ev, kind = "liver")
  ## the same links; their chains follow kind
  expect_equal(d[, 1:12], a[kind == "direct", 1:12], ignore_attr = TRUE)
  expect_equal(v[, 1:12], a[kind == "liver", 1:12], ignore_attr = TRUE)
  expect_identical(nrow(d) + nrow(v), nrow(a))
  expect_true(all(v$inoculation_tick < v$target_tick))
  expect_true(all(d$inoculation_tick == d$target_tick))
  ## the direct links are the oocyst's parents joined to the blood stage its products started
  up <- ev[event_type %in% c("hv_uptake_clonal", "hv_uptake_recombinant"),
           .(oocyst_id, p1 = parent1_pedigree_id, p2 = parent2_pedigree_id)]
  m <- merge(ev[event_type == "vh_inoculation" & is_dormant == 0 & oocyst_id != -1, .(oocyst_id, target = pedigree_id)],
             up, by = "oocyst_id")
  old <- rbind(m[p1 != -1, .(source = p1, target)], m[p2 != -1, .(source = p2, target)])
  expect_setequal(paste(old$source, old$target), paste(d$source_pedigree_id, d$target_pedigree_id))
  expect_error(transmission_links(ev, kind = "both"), "'arg' should be one of")
})

test_that("a log with no transmission gives no link, in the same columns", {
  full <- transmission_links(load_events(transmission_fixture_run()))
  no_take <- make_events(list(list(tick = 3L, event_type = "vh_inoculation_import_no_take", human_id = 1L)))
  for (ev in list(make_events(), no_take, load_events(fake_run()))) {
    l <- transmission_links(ev)
    expect_identical(nrow(l), 0L)
    expect_identical(tl_classes(l), tl_classes(full))
  }
})

test_that("a cut event log warns, and a run or a table without the columns is refused", {
  r <- transmission_fixture_run(); ev <- load_events(r)
  expect_warning(l <- transmission_links(ev[tick >= 20]),
                 "3 infected genomes have no feed, 2 links lack a bite or onset day): pass the whole event log",
                 fixed = TRUE)
  expect_false(31L %in% l$target_pedigree_id)
  expect_no_warning(transmission_links(ev))
  expect_error(transmission_links(r), "for a run, pass load_events(run)", fixed = TRUE)
  expect_error(transmission_links(data.table::copy(ev)[, oocyst_id := NULL]), "lacks the column(s) oocyst_id",
               fixed = TRUE)
})

test_that("on a run, every local case is the target of one transmission, and each link holds together", {
  skip_if_no_slim()
  r <- links_run()
  ev <- load_events(r); l <- transmission_links(ev); cs <- classify_cases(ev)
  batch <- ev[event_type %in% c("vh_inoculation", "vh_inoculation_import") & is_dormant == 1]
  act <- ev[event_type == "hypnozoite_activation"]
  primary <- ev[event_type == "vh_inoculation" & is_dormant == 0]
  batch_only <- batch[event_type == "vh_inoculation" & !paste(human_id, oocyst_id) %in%
                        paste(primary$human_id, primary$oocyst_id), pedigree_id]
  ## the run holds every kind of link
  expect_gt(l[kind == "liver", .N], 50)
  expect_gt(l[weight == 0.5, uniqueN(target_pedigree_id)], 10)
  expect_gt(act[parent1_pedigree_id %in% batch[event_type == "vh_inoculation_import", pedigree_id], .N], 10)
  expect_gt(l[source_pedigree_id %in% ev[event_type == "recrudescence", pedigree_id], .N], 10)
  expect_gt(l[target_pedigree_id %in% act[parent1_pedigree_id %in% batch_only, pedigree_id], .N], 10)
  expect_gt(l[source_human_id == target_human_id, .N], 0)
  local <- cs[origin == "local reinfection" | (origin == "relapse" & !imported), pedigree_id]
  expect_setequal(unique(l$target_pedigree_id), local)
  expect_true(all(l[, sum(weight), by = target_pedigree_id]$V1 == 1))
  expect_true(all(l$weight %in% c(0.5, 1)))
  expect_setequal(l[kind == "liver", unique(target_pedigree_id)], act[oocyst_id != -1, pedigree_id])
  expect_true(all(l$source_tick <= l$uptake_tick))
  expect_true(all(l$inoculation_tick - l$uptake_tick >= load_params(r)$EIP_DAYS))
  expect_true(all(l[kind == "direct", inoculation_tick == target_tick]))
  expect_true(all(l[kind == "liver", target_tick > inoculation_tick]))
  ## the mosquito that fed is the one that bit, and the source's onset is in the host it fed on
  bit <- rbind(primary[, .(target_pedigree_id = pedigree_id, bit = mosq_id)],
               batch[act, on = .(pedigree_id = parent1_pedigree_id), nomatch = NULL,
                     .(target_pedigree_id = i.pedigree_id, bit = x.mosq_id)])
  expect_true(all(bit[l, on = "target_pedigree_id"][, bit == mosq_id]))
  onset <- ev[event_type %in% malariafwd:::.onset_events() & is_dormant == 0, .(source_pedigree_id = pedigree_id, h = human_id)]
  expect_true(all(onset[l, on = "source_pedigree_id"][, h == source_human_id]))
})

## rc_by_origin(): the case reproduction number by true origin.
rc_row <- function(rc, o) unlist(rc[origin == o, .(n_cases, n_secondary, Rc, onward_frac)])
rc_want <- function(...) stats::setNames(c(...), c("n_cases", "n_secondary", "Rc", "onward_frac"))

test_that("rc_by_origin() counts each secondary case once, through bites and the relapses they leave", {
  r <- transmission_fixture_run()
  rc <- rc_by_origin(r)
  expect_identical(names(rc), c("origin", "n_cases", "n_secondary", "Rc", "onward_frac"))
  expect_identical(levels(rc$origin), c("imported", "local reinfection", "relapse", "recrudescence"))
  expect_identical(rc$n_cases, c(1L, 7L, 5L, 1L))
  expect_equal(rc_row(rc, "imported"), rc_want(1, 2, 2, 1))
  expect_equal(rc_row(rc, "local reinfection"), rc_want(7, 4, 4 / 7, 3 / 7))
  expect_equal(rc_row(rc, "relapse"), rc_want(5, 0, 0, 0))
  expect_equal(rc_row(rc, "recrudescence"), rc_want(1, 1, 1, 1))
  ## a recrudescence's transmission belongs to the clone it came from
  ac <- rc_by_origin(r, unit = "acquisition")
  expect_identical(levels(ac$origin), c("imported", "local reinfection", "relapse"))
  expect_equal(rc_row(ac, "imported"), rc_want(1, 2, 2, 1))
  expect_equal(rc_row(ac, "local reinfection"), rc_want(7, 5, 5 / 7, 4 / 7))
  expect_equal(rc_row(ac, "relapse"), rc_want(5, 0, 0, 0))
})

test_that("the direct and liver links add up to all of them", {
  r <- transmission_fixture_run()
  for (u in c("episode", "acquisition"))
    expect_equal(rc_by_origin(r, u, kind = "direct")$n_secondary + rc_by_origin(r, u, kind = "liver")$n_secondary,
                 rc_by_origin(r, u)$n_secondary)
  lv <- rc_by_origin(r, kind = "liver")
  expect_equal(lv$n_secondary, c(1, 1, 0, 0))
  expect_equal(rc_row(lv, "local reinfection"), rc_want(7, 1, 1 / 7, 2 / 7))
})

test_that("a run with no case gives every origin, with no case and no Rc, without a warning", {
  no_take <- make_events(list(list(tick = 3L, event_type = "vh_inoculation_import_no_take", human_id = 1L)))
  for (r in list(fake_run(), transmission_fixture_run(no_take))) {
    expect_no_warning(rc <- rc_by_origin(r))
    expect_identical(rc$n_cases, rep(0L, 4))
    expect_identical(rc$n_secondary, rep(0, 4))
    expect_identical(c(rc$Rc, rc$onward_frac), rep(NA_real_, 8))
    expect_identical(levels(rc$origin), c("imported", "local reinfection", "relapse", "recrudescence"))
  }
})

test_that("from and followup choose the onsets counted", {
  r <- transmission_fixture_run()
  f <- rc_by_origin(r, followup = 60)
  expect_equal(rc_row(f, "imported"), rc_want(1, 2, 2, 1))
  expect_equal(rc_row(f, "local reinfection"), rc_want(6, 4, 2 / 3, 1 / 2))
  expect_equal(rc_row(f, "relapse"), rc_want(2, 0, 0, 0))
  expect_equal(rc_row(f, "recrudescence"), rc_want(0, 0, NA, NA))
  f <- rc_by_origin(r, from = 40)
  expect_equal(rc_row(f, "imported"), rc_want(0, 0, NA, NA))
  expect_equal(rc_row(f, "local reinfection"), rc_want(4, 2, 0.5, 0.5))
  expect_equal(rc_row(f, "relapse"), rc_want(5, 0, 0, 0))
  expect_equal(rc_row(f, "recrudescence"), rc_want(1, 1, 1, 1))
  expect_error(rc_by_origin(r, followup = 150), "no case can be counted: onsets from day 1 to day 0")
  expect_error(rc_by_origin(r, from = 151), "no case can be counted")
  for (v in list(-1, 2.5, c(1, 2))) expect_error(rc_by_origin(r, followup = v), "`followup` must be")
  expect_error(rc_by_origin(r, from = 2.5), "`from` must be")
  expect_error(rc_by_origin(load_events(r)), "`run` must be a vivax_run")
})

test_that("on a run, Rc's counts are the links' weights, and a window only removes cases", {
  skip_if_no_slim()
  r <- links_run(); ev <- load_events(r)
  l <- transmission_links(ev); cs <- classify_cases(ev)
  for (u in c("episode", "acquisition"))
    expect_equal(rc_by_origin(r, u, kind = "direct")$n_secondary + rc_by_origin(r, u, kind = "liver")$n_secondary,
                 rc_by_origin(r, u)$n_secondary)
  rc <- rc_by_origin(r)
  expect_equal(sum(rc$n_secondary), sum(l[source_pedigree_id %in% cs$pedigree_id]$weight))
  expect_identical(sum(rc$n_cases), nrow(cs))
  n <- vapply(c(0, 30, 90), function(f) sum(rc_by_origin(r, followup = f)$n_cases), 0L)
  expect_true(all(diff(n) <= 0))
  n <- vapply(c(0, 50, 100), function(f) sum(rc_by_origin(r, from = f)$n_cases), 0L)
  expect_true(all(diff(n) <= 0))
})

## transmission_links()'s chains.

## a plain union-find with path halving, as a reference; each node's root, in the order of
## unique(c(a, b))
ref_components <- function(a, b) {
  node <- unique(c(a, b)); p <- seq_along(node)
  find <- function(x) { while (p[x] != x) { p[x] <<- p[p[x]]; x <- p[x] }; x }
  ia <- match(a, node); ib <- match(b, node)
  for (k in seq_along(ia)) { ra <- find(ia[k]); rb <- find(ib[k]); if (ra != rb) p[max(ra, rb)] <- min(ra, rb) }
  vapply(seq_along(node), find, integer(1))
}

test_that("the links carry their chains, which follow kind and keep a recrudescence with its clone", {
  ev <- load_events(transmission_fixture_run())
  l <- transmission_links(ev)
  expect_identical(names(l)[13:14], c("chain_id", "chain_size"))
  ch <- function(k) unique(transmission_links(ev, kind = k)[, .(chain_id, chain_size)])
  ## the outcrossed oocyst joins 15's chain and 16's, and the recrudescence 12 is in 6's
  expect_identical(ch("all")$chain_id, 1:2); expect_identical(ch("all")$chain_size, c(12L, 2L))
  expect_identical(ch("direct")$chain_size, c(8L, 2L))
  expect_identical(ch("liver")$chain_id, 1:3); expect_identical(ch("liver")$chain_size, c(3L, 3L, 2L))
  expect_identical(l[target_pedigree_id == 31L]$chain_id, 2L)
  expect_true(all(l[target_pedigree_id != 31L]$chain_id == 1L))
  e <- transmission_links(make_events())
  expect_identical(c(nrow(e), class(e$chain_id), class(e$chain_size)), c("0", "integer", "integer"))
})

test_that(".link_components() finds the components of long and of random graphs", {
  lc <- malariafwd:::.link_components
  set.seed(1)
  n <- 100000; ord <- sample(n)
  t0 <- proc.time()[["elapsed"]]
  expect_identical(data.table::uniqueN(lc(ord[-n], ord[-1])$chain), 1L)
  expect_lt(proc.time()[["elapsed"]] - t0, 10)
  expect_identical(data.table::uniqueN(lc(c(1, 2, 10, 11), c(2, 3, 11, 12))$chain), 2L)
  expect_identical(data.table::uniqueN(lc(c(1, 2, 10, 11, 3), c(2, 3, 11, 12, 10))$chain), 1L)
  expect_identical(data.table::uniqueN(lc(rep(0, 5), 1:5)$chain), 1L)
  for (g in 1:200) {
    k <- sample(2:300, 1); a <- sample(0:(k - 1), k, TRUE); b <- sample(0:(k - 1), k, TRUE)
    x <- lc(a, b)$chain; y <- ref_components(a, b)
    expect_identical(match(x, unique(x)), match(y, unique(y)))
  }
})

test_that("on a run, the chains partition the linked genomes, a recrudescence with its clone", {
  skip_if_no_slim()
  ev <- load_events(links_run()); l <- transmission_links(ev)
  src <- malariafwd:::.acquisition_root(l$source_pedigree_id, ev)
  node <- unique(c(src, l$target_pedigree_id))
  ch <- unique(l[, .(chain_id, chain_size)])
  expect_identical(sum(ch$chain_size), length(node))
  expect_true(all(ch$chain_size >= 2L))
  expect_false(any(node %in% ev[event_type == "recrudescence", pedigree_id]))
  y <- ref_components(src, l$target_pedigree_id)[match(src, node)]
  expect_identical(l$chain_id, match(y, unique(y)))
})
