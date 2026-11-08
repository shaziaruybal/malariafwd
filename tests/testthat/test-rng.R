## Seeded draws: the same under any RNG kind, and the caller's RNG left as it was.

test_that("a seeded draw is the same under any RNG kind, and leaves the caller's kinds and stream alone", {
  kinds <- RNGkind()
  on.exit(suppressWarnings(RNGkind(kinds[1], kinds[2], kinds[3])), add = TRUE)
  f <- function() malariafwd:::.sample_fixed(1:1000, 5, seed = 42)
  set.seed(1); before <- .Random.seed
  a <- f()
  expect_identical(.Random.seed, before)
  RNGkind("L'Ecuyer-CMRG"); set.seed(1); before <- .Random.seed
  expect_identical(f(), a)
  expect_identical(RNGkind()[1], "L'Ecuyer-CMRG")
  expect_identical(.Random.seed, before)
  suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rounding")); set.seed(1)
  expect_identical(f(), a)
  expect_identical(RNGkind()[3], "Rounding")
})

test_that("a seeded draw leaves no RNG state where the caller had none", {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv())
  on.exit(if (had) assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  if (had) rm(".Random.seed", envir = globalenv())
  malariafwd:::.sample_fixed(1:1000, 5, seed = 42)
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
})

test_that("new_infection_diversity() gives the same row inside a future worker's RNG kind", {
  skip_if_no_slim(); skip_if_no_tskit()
  kinds <- RNGkind()
  on.exit(suppressWarnings(RNGkind(kinds[1], kinds[2], kinds[3])), add = TRUE)
  r <- links_run()
  a <- new_infection_diversity(r, cap = 20)
  RNGkind("L'Ecuyer-CMRG")
  expect_identical(new_infection_diversity(r, cap = 20), a)
})
