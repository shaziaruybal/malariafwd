## The run directory: a run is written beside its outdir and moved into place when it has finished,
## DONE marks it, and an outdir that already holds a run is reused or replaced only when asked.

md5s <- function(d) tools::md5sum(file.path(d, sort(list.files(d))))

## A stand-in for SLiM that reports version 5.2 and then runs `body` (sh), for the paths that need
## SLiM to fail in a given way.
fake_slim <- function(body) {
  fake <- tempfile("fakeslim")
  writeLines(c("#!/bin/sh", 'if [ "$1" = "-v" ]; then echo "SLiM version 5.2"; exit 0; fi', body), fake)
  Sys.chmod(fake, "755")
  fake
}

test_that("a run stopped part-way leaves its outdir as it was, and its files in <outdir>.failed", {
  skip_if_no_slim()
  od <- file.path(tempfile("mfwd_dir_"), "run")
  suppressMessages(run_sim(vivax_params(N_HOSTS = 40, N_MOSQ = 60, SIM_DAYS = 60, CHROMOSOMES = "1"),
                           seed = 1, outdir = od))
  before <- md5s(od)
  big <- vivax_params(N_HOSTS = 400, N_MOSQ = 800, SIM_DAYS = 3000, CHROMOSOMES = "1")
  expect_error(run_sim(big, seed = 2, outdir = od, overwrite = TRUE, timeout = 2), "timed out.*\\.failed")
  expect_identical(md5s(od), before)
  expect_true(file.exists(file.path(paste0(od, ".failed"), "run.log")))
})

test_that("outdir appears only when the run has finished", {
  skip_on_os("windows")
  old <- Sys.getenv("MALARIAFWD_SLIM", NA)
  Sys.setenv(MALARIAFWD_SLIM = fake_slim(c('echo "half a run" > events.tsv', "exit 3")))
  on.exit(if (is.na(old)) Sys.unsetenv("MALARIAFWD_SLIM") else Sys.setenv(MALARIAFWD_SLIM = old), add = TRUE)
  od <- file.path(tempfile("mfwd_dir_"), "run")
  expect_error(run_sim(vivax_params(), outdir = od), "SLiM run failed \\(exit 3\\)\\. Log: .*run\\.failed/run\\.log")
  expect_false(dir.exists(od))
  expect_true(all(file.exists(file.path(paste0(od, ".failed"), c("events.tsv", "run.log", "model.slim", "cmdline.txt")))))
  expect_identical(list.files(dirname(od)), "run.failed")
})

test_that("SLiM's error output reaches run.log", {
  skip_on_os("windows")
  old <- Sys.getenv("MALARIAFWD_SLIM", NA)
  Sys.setenv(MALARIAFWD_SLIM = fake_slim(c('echo "boom on stderr" >&2', "exit 1")))
  on.exit(if (is.na(old)) Sys.unsetenv("MALARIAFWD_SLIM") else Sys.setenv(MALARIAFWD_SLIM = old), add = TRUE)
  od <- tempfile("mfwd_dir_")
  expect_error(run_sim(vivax_params(), outdir = od), "boom on stderr")
  expect_match(paste(readLines(file.path(paste0(od, ".failed"), "run.log")), collapse = "\n"), "boom on stderr")
})

test_that("an existing run is reused or replaced only when asked", {
  skip_if_no_slim()
  p <- vivax_params(N_HOSTS = 20, N_MOSQ = 40, SIM_DAYS = 5, CHROMOSOMES = "1")
  od <- file.path(tempfile("mfwd_dir_"), "run")
  r1 <- suppressMessages(run_sim(p, seed = 1, outdir = od))
  done <- tools::md5sum(file.path(od, "DONE"))
  e <- tryCatch(run_sim(p, seed = 1, outdir = od), error = identity)
  expect_s3_class(e, "malariafwd_outdir_in_use")
  expect_match(conditionMessage(e), "already holds this run: skip_existing = TRUE returns it")
  expect_message(r2 <- run_sim(p, seed = 1, outdir = od, skip_existing = TRUE), "already holds this run .*not run again")
  expect_identical(r2$seed, 1L)
  expect_identical(unname(tools::md5sum(file.path(od, "DONE"))), unname(done))
  expect_error(run_sim(p, seed = 2, outdir = od, skip_existing = TRUE), "already holds a run of other parameters, seed or model",
               class = "malariafwd_outdir_in_use")
  r3 <- suppressMessages(run_sim(p, seed = 2, outdir = od, overwrite = TRUE))
  expect_identical(r3$seed, 2L)
  expect_identical(list.files(dirname(od), all.files = TRUE, no.. = TRUE), "run")
})

test_that("a directory with files that are not a run's is never used, and ignored files do not count", {
  r <- fake_run()
  writeLines("my notes", file.path(r, "notes.txt"))
  expect_error(run_sim(vivax_params(), outdir = r, overwrite = TRUE),
               "holds files that are not part of a run \\(notes\\.txt\\)", class = "malariafwd_outdir_in_use")
  expect_true(file.exists(file.path(r, "notes.txt")))
  skip_if_no_slim()
  od <- tempfile("mfwd_dir_"); dir.create(od); writeLines("", file.path(od, ".DS_Store"))
  r <- suppressMessages(run_sim(vivax_params(N_HOSTS = 20, N_MOSQ = 40, SIM_DAYS = 5, CHROMOSOMES = "1"),
                                seed = 1, outdir = od))
  expect_true(r$complete)
})

test_that("DONE is written last, keys the run, and lists every file with its size", {
  skip_if_no_slim()
  p <- vivax_params(N_HOSTS = 20, N_MOSQ = 40, SIM_DAYS = 5, CHROMOSOMES = "1")
  r <- suppressMessages(run_sim(p, seed = 3, outdir = tempfile("mfwd_dir_")))
  expect_true(r$complete)
  d <- malariafwd:::.read_done(r$outdir)
  expect_identical(names(d)[1:7], c("finished", "elapsed_s", "seed", "run_key", "model_md5", "malariafwd_version",
                                    "slim_version"))
  files <- sub("^file:", "", names(d)[-(1:7)])
  expect_setequal(files, setdiff(list.files(r$outdir), "DONE"))
  expect_identical(as.numeric(d[-(1:7)]), as.numeric(file.size(file.path(r$outdir, files))))
  expect_identical(d[["seed"]], "3")
  expect_identical(d[["run_key"]], malariafwd:::.run_key(p, 3))
  expect_identical(d[["model_md5"]], unname(tools::md5sum(slim_model())))
  ## the key ignores narration, and nothing else
  expect_identical(malariafwd:::.run_key(vivax_params(N_HOSTS = 20, N_MOSQ = 40, SIM_DAYS = 5, CHROMOSOMES = "1",
                                                      VERBOSE = TRUE), 3), d[["run_key"]])
  expect_false(identical(malariafwd:::.run_key(p, 4), d[["run_key"]]))
  expect_false(identical(malariafwd:::.run_key(vivax_params(N_HOSTS = 21, N_MOSQ = 40, SIM_DAYS = 5,
                                                            CHROMOSOMES = "1"), 3), d[["run_key"]]))
})

test_that("a directory without DONE, or changed since its run finished, is reported once", {
  r <- fake_run(done = FALSE)
  w <- list()
  withCallingHandlers({ load_events(r); load_state(r) },
                      warning = function(x) { w[[length(w) + 1]] <<- x; invokeRestart("muffleWarning") })
  expect_length(w, 1L)
  expect_s3_class(w[[1]], "malariafwd_incomplete_run")
  expect_match(conditionMessage(w[[1]]), "has no DONE file")
  expect_false(suppressWarnings(load_run(r))$complete)
  r <- fake_run()
  writeLines("tick", file.path(r, "state.tsv"))
  expect_warning(run <- load_run(r), "has changed since its run finished: state.tsv has another size",
                 class = "malariafwd_incomplete_run")
  expect_false(run$complete)
  expect_output(print(run), "not marked complete")
  expect_true(load_run(fake_run())$complete)
})

test_that("the run directory repeats its own run", {
  skip_if_no_slim()
  r <- suppressMessages(run_sim(vivax_params(N_HOSTS = 20, N_MOSQ = 40, SIM_DAYS = 20, CHROMOSOMES = "1"),
                                seed = 5, outdir = tempfile("mfwd_dir_")))
  cmd <- readLines(file.path(r$outdir, "cmdline.txt"), encoding = "UTF-8")
  args <- cmd[!startsWith(cmd, "#")]
  expect_true("OUTDIR='.'" %in% args)
  expect_identical(args[length(args)], "model.slim")
  again <- tempfile("mfwd_again_"); dir.create(again)
  file.copy(file.path(r$outdir, "model.slim"), again)
  res <- processx::run(slim_path(), args, wd = again, error_on_status = FALSE)
  expect_identical(res$status, 0L)
  f <- c("events.tsv", "state.tsv", "params.tsv")
  expect_identical(unname(tools::md5sum(file.path(again, f))), unname(tools::md5sum(file.path(r$outdir, f))))
})

test_that("a run's key takes any SIM_ID, in any locale", {
  p <- vivax_params(SIM_ID = "Côte d’Ivoire 1-3")
  k <- malariafwd:::.run_key(p, 1)
  expect_match(k, "^[0-9a-f]{32}$")
  old <- Sys.getlocale("LC_COLLATE")
  on.exit(Sys.setlocale("LC_COLLATE", old), add = TRUE)
  Sys.setlocale("LC_COLLATE", "C")
  expect_identical(malariafwd:::.run_key(p, 1), k)
})
