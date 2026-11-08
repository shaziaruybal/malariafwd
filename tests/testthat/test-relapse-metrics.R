## Relapse metrics: relapse_fraction()'s two readings.

test_that("relapse_fraction() gives relapses over every case and over the acquisitions", {
  rf <- relapse_fraction(episodes_fixture())
  expect_identical(names(rf), c("n_relapse", "n_new_infection", "n_recrudescence",
                                "relapse_frac_episodes", "relapse_frac_acquisitions"))
  expect_identical(c(rf$n_relapse, rf$n_new_infection, rf$n_recrudescence), c(4L, 7L, 1L))
  expect_identical(c(rf$relapse_frac_episodes, rf$relapse_frac_acquisitions), c(4 / 12, 4 / 11))
})

test_that("relapse_fraction() leaves the masked relapse out of both readings when asked", {
  rf <- relapse_fraction(episodes_fixture(), exclude_masked = TRUE)
  expect_identical(rf$n_relapse, 3L)
  expect_identical(c(rf$relapse_frac_episodes, rf$relapse_frac_acquisitions), c(3 / 11, 3 / 10))
})

test_that("a run with no case gives zero counts and no fraction, without a warning", {
  expect_no_warning(rf <- relapse_fraction(fake_run()))
  expect_identical(c(rf$n_relapse, rf$n_new_infection, rf$n_recrudescence), c(0L, 0L, 0L))
  expect_identical(c(rf$relapse_frac_episodes, rf$relapse_frac_acquisitions), c(NA_real_, NA_real_))
})

test_that("relapse_fraction() reads a run as it reads its events", {
  expect_identical(relapse_fraction(episodes_fixture_run()), relapse_fraction(episodes_fixture()))
  expect_identical(relapse_fraction(episodes_fixture_run(), exclude_masked = TRUE),
                   relapse_fraction(episodes_fixture(), exclude_masked = TRUE))
})
