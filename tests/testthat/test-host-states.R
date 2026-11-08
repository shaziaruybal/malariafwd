## host_state_prevalence(): every state on every day of the state table, subpatent carriers a
## state of their own.

hs_state <- function() data.table::data.table(
  tick = c(1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L),
  entity_type = c(rep("human", 5), "mosquito", rep("human", 5)),
  entity_id = c(1:5, 1L, 1:5),
  moi_active = c(1L, 0L, 0L, 0L, 0L, 1L, 1L, 0L, 0L, 0L, 0L),
  moi_subpatent = c(0L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
  moi_dormant = c(0L, 0L, 0L, 2L, 1L, 0L, 1L, 0L, 0L, 0L, 0L))

test_that("host_state_prevalence() has every state on every day, subpatent carriers their own", {
  hs <- host_state_prevalence(hs_state())
  lv <- c("uninfected", "active only", "active + reservoir", "subpatent", "reservoir only (invisible)")
  expect_identical(names(hs), c("tick", "cls", "n", "frac"))
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

test_that("plot_host_state_prevalence() draws every state", {
  skip_if_not_installed("ggplot2")
  expect_no_warning(b <- ggplot2::ggplot_build(plot_host_state_prevalence(host_state_prevalence(hs_state()))))
  expect_identical(length(unique(b$data[[1]]$fill)), 5L)
})
