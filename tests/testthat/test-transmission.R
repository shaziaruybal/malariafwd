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
  expect_equal(l, want, ignore_attr = TRUE)
  expect_identical(tl_classes(l), tl_classes(want))
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
  expect_equal(d, a[kind == "direct"], ignore_attr = TRUE)
  expect_equal(v, a[kind == "liver"], ignore_attr = TRUE)
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
  r <- suppressMessages(run_sim(vivax_params(
    N_HOSTS = 60, N_MOSQ = 120, BITING_RATE = 0.16, SIM_DAYS = 300, LAMBDA_SPOROZOITES = 2, IMPORT_START = 1L,
    IMPORT_END = 300L, IMPORT_RATE_PER_TICK = 0.1, BLOODCURE_COVERAGE = 0.5, P_BLOODCURE_SURVIVAL = 0.3,
    P_RECRUD_GIVEN_FAIL = 0.7, P_PRIMARY = 0.8, CHROMOSOMES = "1"), seed = 2, outdir = tempfile("mfwd_links_")))
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
