## host_state_prevalence(): every state on every day of the state table, subpatent carriers a
## state of their own, and each day's mean reservoir.

hs_state <- function() data.table::data.table(
  tick = c(1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L),
  entity_type = c(rep("human", 5), "mosquito", rep("human", 5)),
  entity_id = c(1:5, 1L, 1:5),
  moi_active = c(1L, 0L, 0L, 0L, 0L, 1L, 1L, 0L, 0L, 0L, 0L),
  moi_subpatent = c(0L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
  moi_dormant = c(0L, 0L, 0L, 2L, 1L, 0L, 1L, 0L, 0L, 0L, 0L),
  hyp_reservoir_size = c(0L, 0L, 0L, 5L, 2L, 9L, 3L, 0L, 0L, 0L, 0L))

test_that("host_state_prevalence() has every state on every day, subpatent carriers their own", {
  hs <- host_state_prevalence(hs_state())
  lv <- c("uninfected", "active only", "active + reservoir", "subpatent", "reservoir only (invisible)")
  expect_identical(names(hs), c("tick", "cls", "n", "frac", "hypnozoites_per_host"))
  expect_identical(hs$tick, rep(1:2, each = 5L))
  expect_identical(hs$cls, factor(rep(lv, 2), levels = lv))
  ## day 1: host 3 subpatent alone, host 4 subpatent with hypnozoites; day 2 has three empty states
  expect_identical(hs$n, c(1L, 1L, 0L, 2L, 1L, 4L, 0L, 1L, 0L, 0L))
  expect_equal(hs$frac, hs$n / 5)
  ns <- host_state_prevalence(hs_state()[, moi_subpatent := NULL])
  expect_identical(ns$n[1:5], c(2L, 1L, 0L, 0L, 2L))
  e <- host_state_prevalence(hs_state()[0])
  expect_identical(nrow(e), 0L)
  expect_identical(levels(e$cls), lv)
})

test_that("host_state_prevalence() gives each day's mean reservoir per host on each of its rows", {
  hs <- host_state_prevalence(hs_state())
  ## day 1: 5 and 2 hypnozoites in hosts 4 and 5; day 2: 3 in host 1; the mosquito row's 9 is left out
  expect_equal(hs$hypnozoites_per_host, rep(c(7, 3) / 5, each = 5L))
  expect_equal(unique(hs[, .(tick, hypnozoites_per_host)]),
               data.table::data.table(tick = 1:2, hypnozoites_per_host = c(7, 3) / 5))
  expect_identical(class(host_state_prevalence(hs_state()[0])$hypnozoites_per_host), "numeric")
})

test_that("plot_host_state_prevalence() draws every state", {
  skip_if_not_installed("ggplot2")
  expect_no_warning(b <- ggplot2::ggplot_build(plot_host_state_prevalence(host_state_prevalence(hs_state()))))
  expect_identical(length(unique(b$data[[1]]$fill)), 5L)
})

test_that("coi_distribution() counts the infected host-days by COI, in a window", {
  st <- fixture_state()
  cd <- coi_distribution(st)
  expect_identical(names(cd), c("moi_active", "n", "frac"))
  expect_identical(vapply(cd, function(v) class(v)[1], ""), c(moi_active = "integer", n = "integer", frac = "numeric"))
  h <- st[entity_type == "human" & moi_active > 0]
  expect_identical(cd$n, h[, .N, by = moi_active][order(moi_active)]$N)
  expect_equal(sum(cd$frac), 1)
  w <- coi_distribution(st, from = 50, to = 70)
  expect_identical(w$n, h[tick >= 50 & tick <= 70, .N, by = moi_active][order(moi_active)]$N)
  expect_identical(coi_distribution(st, from = 50)$n, h[tick >= 50, .N, by = moi_active][order(moi_active)]$N)
  expect_error(coi_distribution(st, from = 70, to = 50), "`from` = 70 is after `to` = 50.", fixed = TRUE)
  expect_identical(nrow(coi_distribution(st[0])), 0L)
})

test_that("event_types() is a data.table", {
  expect_true(data.table::is.data.table(event_types()))
  expect_identical(names(event_types()), c("event_type", "stage", "ends_episode", "is_onset", "description"))
})
