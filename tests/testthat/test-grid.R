## run_grid() checks every set and cell before any cell starts, resumes from what is on disk and
## records the grid in grid.tsv, which load_grid() reads.

tiny <- function(...)
  do.call(vivax_params, utils::modifyList(list(N_HOSTS = 20, N_MOSQ = 40, SIM_DAYS = 3, CHROMOSOMES = "1"), list(...)))

## run_sim() stood in for, counting its calls in `calls` of the calling test
mock_cells <- function(env = parent.frame()) {
  assign("calls", 0L, envir = env)
  local_mocked_bindings(run_sim = function(...) {
    assign("calls", get("calls", envir = env) + 1L, envir = env)
    list(extinct_tick = NA_integer_, never_started = FALSE)
  }, .require_slim = function() "5.2", .env = env)
}

test_that("set names must differ, also ignoring case, before any cell runs", {
  mock_cells()
  p <- vivax_params()
  expect_error(run_grid(list(a = p, a = p), outdir_base = tempfile("mfwd_grid_")), "set names must differ.*: a\\.")
  expect_error(run_grid(list(a = p, A = p), outdir_base = tempfile("mfwd_grid_")), "ignoring case.*: a, A\\.")
  expect_identical(calls, 0L)
})

test_that("blank names are filled by position, and a name that is no directory's stops", {
  mock_cells()
  p <- vivax_params()
  g <- suppressMessages(run_grid(list(a = p, p), outdir_base = tempfile("mfwd_grid_")))
  expect_identical(unique(g$set), c("a", "set02"))
  for (nm in c("a/b", "..", " a", "a.", "a:b"))
    expect_error(run_grid(stats::setNames(list(p), nm), outdir_base = tempfile("mfwd_grid_")),
                 "must not hold", label = nm)
})

test_that("seeds must be distinct whole numbers", {
  mock_cells()
  for (s in list(2.5, c(1, 1), NA, numeric(0), "1", Inf))
    expect_error(run_grid(vivax_params(), seeds = s, outdir_base = tempfile("mfwd_grid_")), "`seeds` must",
                 label = deparse(s))
  expect_identical(calls, 0L)
})

test_that("a grid leaves the caller's random-number stream as it was", {
  mock_cells()
  set.seed(42); before <- .Random.seed
  suppressMessages(run_grid(list(a = vivax_params(), b = vivax_params(N_HOSTS = 50)), seeds = 1:3,
                            outdir_base = tempfile("mfwd_grid_")))
  expect_identical(.Random.seed, before)
  expect_identical(calls, 6L)
})

test_that("a grid resumes: cells that hold their run are not run again", {
  skip_if_no_slim()
  base <- tempfile("mfwd_grid_")
  sets <- list(a = tiny(), b = tiny(N_HOSTS = 25))
  g1 <- suppressMessages(run_grid(sets, seeds = 1:2, outdir_base = base))
  expect_true(all(g1$ok) && !any(g1$skipped))
  t1 <- file.mtime(file.path(g1$outdir, "DONE"))
  expect_message(g2 <- run_grid(sets, seeds = 1:2, outdir_base = base), "4 of 4 grid cell\\(s\\) already held their run")
  expect_true(all(g2$skipped))
  expect_identical(file.mtime(file.path(g2$outdir, "DONE")), t1)
  g3 <- suppressMessages(run_grid(sets, seeds = 1:3, outdir_base = base))
  expect_identical(g3$skipped, rep(c(TRUE, TRUE, FALSE), 2))
  ## an edited set holds another run in its cells: stopped up front, and run again only when asked
  sets$b <- tiny(N_HOSTS = 30)
  expect_error(run_grid(sets, seeds = 1:3, outdir_base = base),
               "did not start: 3 of 6 cells already hold a different run.*b_seed1, b_seed2, b_seed3")
  g4 <- suppressMessages(run_grid(sets, seeds = 1:3, outdir_base = base, overwrite = TRUE))
  expect_identical(g4$skipped, rep(c(TRUE, FALSE), each = 3))
  expect_identical(unique(g4$N_HOSTS), c(20L, 30L))
})

test_that("grid.tsv records every cell and parameter, and load_grid() reads it back", {
  skip_if_no_slim()
  base <- tempfile("mfwd_grid_")
  sets <- list(low = tiny(BITING_RATE = 0.05), high = tiny(BITING_RATE = 0.12, SIM_ID = "NA"))
  g <- suppressMessages(run_grid(sets, seeds = 1:2, outdir_base = base))
  expect_s3_class(g, "data.table")
  expect_identical(names(g), c("set", "seed", "outdir", "ok", "skipped", "error", "extinct_tick", "never_started",
                               "elapsed_s", "finished", "slim_version", "malariafwd_version", "SIM_ID", "BITING_RATE"))
  expect_identical(vapply(g, function(x) class(x)[1], ""),
                   c(set = "character", seed = "integer", outdir = "character", ok = "logical", skipped = "logical",
                     error = "character", extinct_tick = "integer", never_started = "logical", elapsed_s = "numeric",
                     finished = "character", slim_version = "character", malariafwd_version = "character",
                     SIM_ID = "character", BITING_RATE = "numeric"))
  m <- data.table::fread(file.path(base, "grid.tsv"), colClasses = "character")
  expect_identical(nrow(m), 4L)
  expect_identical(names(m)[1:14], c("set", "seed", "dir", "ok", "skipped", "error", "extinct_tick", "never_started",
                                     "elapsed_s", "finished", "slim_version", "malariafwd_version", "run_key", "model_md5"))
  expect_identical(names(m)[-(1:14)], names(vivax_params()))
  expect_equal(load_grid(base), g)
  expect_identical(names(load_grid(base, params = "all"))[-(1:12)], names(vivax_params()))
  ## a grid moved elsewhere still loads, and a call with other seeds keeps the cells of the first
  moved <- tempfile("mfwd_moved_")
  expect_true(file.rename(base, moved))
  lg <- load_grid(moved)
  expect_true(all(dir.exists(lg$outdir)))
  suppressMessages(run_grid(sets, seeds = 3, outdir_base = moved))
  expect_identical(sort(paste(load_grid(moved)$set, load_grid(moved)$seed)),
                   sort(paste(rep(c("high", "low"), each = 3), 1:3)))
})

test_that("load_grid() says what to do without grid.tsv, and warns for a cell that lost its run", {
  d <- tempfile("mfwd_nogrid_"); dir.create(d)
  expect_error(load_grid(d), "no grid.tsv in .*Call run_grid\\(\\) again with the same sets and seeds")
  skip_if_no_slim()
  base <- tempfile("mfwd_grid_")
  g <- suppressMessages(run_grid(tiny(), seeds = 1:2, outdir_base = base))
  unlink(file.path(g$outdir[1], "DONE"))
  expect_warning(load_grid(base), "1 cell\\(s\\) recorded as ok no longer hold their whole run: set01_seed1")
})

test_that("a cell that fails leaves its directory as it was, and its log in <cell>.failed", {
  skip_on_os("windows")
  fake <- tempfile("fakeslim")
  writeLines(c("#!/bin/sh", 'if [ "$1" = "-v" ]; then echo "SLiM version 5.2"; exit 0; fi', "echo boom >&2", "exit 1"), fake)
  Sys.chmod(fake, "755")
  old <- Sys.getenv("MALARIAFWD_SLIM", NA)
  Sys.setenv(MALARIAFWD_SLIM = fake)
  on.exit(if (is.na(old)) Sys.unsetenv("MALARIAFWD_SLIM") else Sys.setenv(MALARIAFWD_SLIM = old), add = TRUE)
  base <- tempfile("mfwd_grid_")
  expect_warning(g <- run_grid(vivax_params(), seeds = 1, outdir_base = base),
                 "1 of 1 grid cell\\(s\\) FAILED: their directories were left as they were")
  expect_false(g$ok)
  expect_match(g$error, "set01_seed1\\.failed/run\\.log")
  expect_true(file.exists(file.path(base, "set01_seed1.failed", "run.log")))
  expect_false(dir.exists(file.path(base, "set01_seed1")))
})

test_that("grid.tsv keeps a quote in a text field, however often it is written", {
  base <- tempfile("mfwd_gridq_"); dir.create(base)
  row <- function(seed) data.table::data.table(set = "lo", seed = seed, dir = paste0("lo_seed", seed), ok = FALSE,
    skipped = FALSE, error = 'a "quoted" failure: ""', extinct_tick = NA_integer_, never_started = NA, elapsed_s = 1.5,
    finished = "2026-11-08 03:00:00", slim_version = "5.2", malariafwd_version = "0.1.0", run_key = "k",
    model_md5 = "m", SIM_ID = 'say "hi"')
  malariafwd:::.write_grid_manifest(base, row(1L))
  r1 <- malariafwd:::.read_grid_manifest(base)
  expect_identical(c(r1$error, r1$SIM_ID), c('a "quoted" failure: ""', 'say "hi"'))
  malariafwd:::.write_grid_manifest(base, row(2L))
  malariafwd:::.write_grid_manifest(base, row(3L))
  r3 <- malariafwd:::.read_grid_manifest(base)
  expect_identical(r3$SIM_ID, rep('say "hi"', 3))
  expect_identical(r3$error, rep('a "quoted" failure: ""', 3))
})
