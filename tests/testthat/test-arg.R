## The ARG readers, also on a run that recorded no genome.

test_that("the ARG readers give their typed empty results on a run that recorded no genome", {
  skip_if_no_slim(); skip_if_no_tskit()
  nv <- never_run(); full <- links_run()
  cls <- function(d) vapply(d, function(v) class(v)[1], "")
  s <- arg_summary(nv)
  expect_identical(c(s$num_nodes, s$num_edges, s$num_samples, s$roots_min, s$roots_max), c(0L, 0L, 0L, 0L, 0L))
  expect_identical(s$roots_median, 0)
  expect_identical(c(s$pct_trees_coalesced, s$max_root_time), c(NA_real_, NA_real_))
  expect_identical(cls(s), cls(arg_summary(full)))
  t1 <- arg_tmrca(nv, 1L, 2L)
  expect_identical(nrow(t1), 1L)
  expect_true(all(is.na(unlist(t1[, 3:5]))))
  t0 <- arg_tmrca(nv, integer(), integer())
  expect_identical(names(t0), c("a", "b", "coalesced", "tmrca", "tmrca_node_time"))
  expect_identical(nrow(t0), 0L)
  rel <- load_events(full)[event_type == "hypnozoite_activation"][1:3]
  t2 <- arg_tmrca(full, rel$pedigree_id, rel$parent1_pedigree_id)
  expect_identical(cls(t0), cls(t2))
  expect_true(all(t2$coalesced > 0.99) && all(abs(t2$tmrca) < 1e-9))
})

test_that("ids none of which is a genome are refused, as are pairs of unequal length", {
  skip_if_no_slim(); skip_if_no_tskit()
  full <- links_run()
  expect_error(arg_summary(full, ids = -5L), "arg_summary(): none of the 1 ids given is a genome", fixed = TRUE)
  expect_error(arg_tmrca(full, -5L, -6L), "arg_tmrca(): none of the 2 ids given is a genome", fixed = TRUE)
  expect_error(arg_tmrca(full, 1:2, 1), "`ids_a` and `ids_b` must have the same length")
  expect_error(arg_tmrca(full, c(1, NA), 1:2), "`ids_a` must be pedigree ids")
})
