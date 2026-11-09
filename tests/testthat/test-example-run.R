## The bundled run: whole, made by this package's model and defaults, and holding every kind of
## case, recurrence and link the examples show. When the model or a default changes, make it again
## with data-raw/example_run.R.

test_that("example_run() gives a whole run, the same one on every call", {
  run <- example_run()
  expect_s3_class(run, "vivax_run")
  expect_true(run$complete)
  expect_false(run$never_started)
  expect_true(is.na(run$extinct_tick))
  expect_identical(example_run()$outdir, run$outdir)
})

test_that("the example run is this model's run of its settings (or run data-raw/example_run.R)", {
  run <- example_run()
  spec <- malariafwd:::.example_run_spec
  expect_identical(as.integer(run$seed), spec$seed)
  want <- unclass(do.call(vivax_params, spec$params))
  got <- unclass(load_params(run))
  keep <- setdiff(names(want), c("SEED", "OUTDIR"))
  expect_equal(got[keep], want[keep])
  done <- malariafwd:::.read_done(run$outdir)
  expect_identical(unname(done[["model_md5"]]), unname(tools::md5sum(slim_model())))
})

test_that("the example run holds every kind of case, recurrence and transmission link", {
  ev <- load_events(example_run())
  cs <- classify_cases(ev)
  expect_setequal(unique(as.character(cs$origin)), c("imported", "local reinfection", "relapse", "recrudescence"))
  expect_true(any(cs$origin == "relapse" & cs$imported))
  expect_setequal(unique(transmission_links(ev)$kind), c("direct", "liver"))
  expect_length(unique(classify_chain(ev)$who_class), 4L)
  expect_true(any(ev$is_recombinant == 1L))
})

test_that("the example run's genealogy gives every rung of the relatedness ladder", {
  skip_if_no_tskit()
  tp <- ibd_truth_pairs(example_run())
  expect_setequal(as.character(unique(tp$cls)),
                  c("self (same genome)", "clonal copies", "relapse vs its origin", "same lineage, 2 hosts",
                    "same-oocyst sibs", "distinct founders"))
})
