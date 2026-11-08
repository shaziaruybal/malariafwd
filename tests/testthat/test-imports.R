## classify_chain(): each case's WHO class from its true chain of transmission from an import;
## classify_travel(): what a travel rule says instead.

chain_fixture <- function() {
  ## the transmission fixture, with the import 7 suppressed on day 40 and recrudescing as 19
  extra <- make_events(list(
    list(tick = 40L, event_type = "bloodcure_subpatent", human_id = 3L, lineage_id = 7L, pedigree_id = 7L,
         is_imported = 1L),
    list(tick = 52L, event_type = "recrudescence", human_id = 3L, lineage_id = 7L, pedigree_id = 19L,
         parent1_pedigree_id = 7L, is_imported = 1L, detected = 1L)))
  fx <- rbind(transmission_fixture(), extra)[order(tick)]
  fx[, event_id := paste0(tick, "_", seq_len(.N)), by = tick]
  fx[]
}

test_that("classify_chain() classes each case by its shortest chain from an import", {
  cc <- classify_chain(chain_fixture())
  expect_identical(cc$pedigree_id, c(2L, 31L, 7L, 14L, 15L, 16L, 19L, 4L, 6L, 9L, 12L, 5L, 13L, 18L, 11L))
  lc <- "indigenous (local lineage)"; fi <- "indigenous (from an import)"
  expect_identical(cc$who_class, c(lc, lc, "imported", lc, lc, "introduced", "imported", lc, fi, "imported", fi, lc,
                                   fi, fi, "introduced"))
  expect_identical(cc$generation, c(NA, NA, 0L, NA, NA, 1L, 0L, NA, 2L, 0L, 2L, NA, 3L, 2L, 1L))
  ## an imported relapse begins from its batch's import row
  expect_identical(cc$import_pedigree_id, c(NA, NA, 7L, NA, NA, 7L, 7L, NA, 7L, 8L, 7L, NA, 7L, 7L, 7L))
  expect_identical(cc$from_import, !is.na(cc$generation))
  expect_identical(cc$who_class == "imported", cc$imported)
  expect_identical(names(cc), c(names(classify_cases(chain_fixture())), "from_import", "generation",
                                "import_pedigree_id", "who_class"))
})

test_that("classify_chain() of a log with no case gives no row, in the same columns", {
  e <- classify_chain(make_events())
  full <- classify_chain(chain_fixture())
  expect_identical(nrow(e), 0L)
  expect_identical(vapply(e, function(v) class(v)[1], ""), vapply(full, function(v) class(v)[1], ""))
  expect_identical(classify_chain(transmission_fixture_run()), classify_chain(transmission_fixture()))
})

test_that("on a run, each generation is one link from the one before", {
  skip_if_no_slim()
  ev <- load_events(links_run()); cc <- classify_chain(ev)
  expect_identical(cc$who_class == "imported", cc$imported)
  expect_gt(cc[generation >= 1L, .N], 10)
  l <- transmission_links(ev)
  src <- malariafwd:::.acquisition_root(l$source_pedigree_id, ev)
  ## each source's generation (a founder, no case, has none); a case's is one more than its best
  g <- data.table::data.table(target = l$target_pedigree_id, sg = as.numeric(cc$generation[match(src, cc$pedigree_id)]))
  best <- g[!is.na(sg), .(sg = min(sg)), by = target]
  chk <- best[cc[generation >= 1L & origin != "recrudescence", .(target = pedigree_id, generation)], on = "target"]
  expect_true(all(chk$generation == chk$sg + 1))
})

## Host 1: an import (1) on day 10, local cases 2 and 3 on days 50 and 100; host 2: an import that
## left hypnozoites only (4) on day 20 and their relapse (5) on day 30; host 3: a local case (6) on
## day 40.
travel_fixture <- function() make_events(list(
  list(tick = 10L, event_type = "vh_inoculation_import", human_id = 1L, lineage_id = 1L, pedigree_id = 1L,
       is_imported = 1L, detected = 1L),
  list(tick = 50L, event_type = "vh_inoculation", human_id = 1L, lineage_id = 2L, pedigree_id = 2L, detected = 1L),
  list(tick = 100L, event_type = "vh_inoculation", human_id = 1L, lineage_id = 3L, pedigree_id = 3L, detected = 1L),
  list(tick = 20L, event_type = "vh_inoculation_import", human_id = 2L, lineage_id = 4L, pedigree_id = 4L,
       is_imported = 1L, is_dormant = 1L),
  list(tick = 30L, event_type = "hypnozoite_activation", human_id = 2L, lineage_id = 4L, pedigree_id = 5L,
       parent1_pedigree_id = 4L, is_imported = 1L, detected = 1L),
  list(tick = 40L, event_type = "vh_inoculation", human_id = 3L, lineage_id = 6L, pedigree_id = 6L, detected = 1L)))

test_that("classify_travel() flags a case whose host had an import in the window up to it", {
  ev <- travel_fixture()
  imp <- "imported (recent travel)"; loc <- "local (no recent travel)"
  by_ped <- function(w) { x <- classify_travel(ev, window = w); x$surveillance[order(x$pedigree_id)] }
  ## cases 1, 2, 3, 5 and 6
  expect_identical(by_ped(60), c(imp, imp, loc, imp, loc))
  expect_identical(by_ped(120), c(imp, imp, imp, imp, loc))
  expect_identical(by_ped(0), rep(loc, 5))
  expect_identical(by_ped(1), c(imp, loc, loc, loc, loc))
  x <- classify_travel(ev)
  cs <- classify_cases(ev)
  expect_identical(x[, names(cs), with = FALSE], cs)
  expect_identical(names(x), c(names(cs), "recent_travel", "surveillance"))
  expect_identical(x$recent_travel, x$surveillance == imp)
  expect_identical(classify_travel(transmission_fixture_run()), classify_travel(transmission_fixture()))
})

test_that("classify_travel() gives the same columns without imports or cases, and checks its window", {
  ev <- travel_fixture()
  cls <- function(d) vapply(d, function(v) class(v)[1], "")
  full <- classify_travel(ev)
  none <- classify_travel(ev[event_type != "vh_inoculation_import"])
  expect_identical(none$surveillance, rep("local (no recent travel)", 4))
  expect_identical(cls(none), cls(full))
  e <- classify_travel(make_events())
  expect_identical(nrow(e), 0L)
  expect_identical(cls(e), cls(full))
  for (w in list(-1, 1.5, NA, c(30, 60)))
    expect_error(classify_travel(ev, window = w), "`window` must be a single whole number of days, 0 or more.",
                 fixed = TRUE)
})

test_that("classify_travel() counts an import once, on the day it arrived", {
  ## host 4: an import (inoculation 7) on day 100, its hypnozoite batch logged on day 101, and the
  ## batch's relapse on day 160
  ev <- make_events(list(
    list(tick = 100L, event_type = "vh_inoculation_import", human_id = 4L, lineage_id = 9L, pedigree_id = 41L,
         inoculation_id = 7L, is_imported = 1L, detected = 1L),
    list(tick = 101L, event_type = "vh_inoculation_import", human_id = 4L, lineage_id = 9L, pedigree_id = 42L,
         inoculation_id = 7L, is_imported = 1L, is_dormant = 1L),
    list(tick = 160L, event_type = "hypnozoite_activation", human_id = 4L, lineage_id = 9L, pedigree_id = 43L,
         parent1_pedigree_id = 42L, is_imported = 1L, detected = 1L)))
  relapse <- function(w) classify_travel(ev, window = w)[pedigree_id == 43L, surveillance]
  expect_identical(relapse(60), "local (no recent travel)")   # days 101 to 160 leave out day 100
  expect_identical(relapse(61), "imported (recent travel)")
  expect_identical(classify_travel(ev, window = 1)[pedigree_id == 41L, surveillance], "imported (recent travel)")
})
