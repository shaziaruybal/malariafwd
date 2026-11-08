## classify_chain(): each case's WHO class from its true chain of transmission from an import.

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
