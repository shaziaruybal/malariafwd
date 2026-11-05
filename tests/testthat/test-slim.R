test_that("a SLiM older than 5.2 is refused, and an unreadable version is let through", {
  expect_false(malariafwd:::.slim_version_ok("5.1"))
  expect_false(malariafwd:::.slim_version_ok("4.3"))
  expect_true(malariafwd:::.slim_version_ok("5.2"))
  expect_true(malariafwd:::.slim_version_ok("6.0"))
  expect_true(malariafwd:::.slim_version_ok(NA_character_))
})

test_that("params.tsv records the seed and the SLiM and package versions of the run", {
  skip_if_not(slim_available(), "SLiM not on PATH")
  run <- run_sim(vivax_params(N_HOSTS = 20, N_MOSQ = 40, SIM_DAYS = 5), seed = 5,
                 outdir = file.path(tempdir(), "mfwd_test_versions"))
  p <- load_params(run)
  expect_equal(p$ACTUAL_SEED, "5")
  expect_equal(p$SLIM_VERSION, slim_version())
  expect_equal(p$MALARIAFWD_VERSION, as.character(utils::packageVersion("malariafwd")))
})
