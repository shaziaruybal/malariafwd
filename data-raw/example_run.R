## Makes inst/extdata/example_run.tar.xz, the run example_run() returns, from the settings and seed
## in R/example.R. Run it from the package's root, with SLiM 5.2 on the PATH, whenever the model or
## a default changes (the tests in tests/testthat/test-example-run.R say when):
##
##   Rscript data-raw/example_run.R

pkgload::load_all(".", quiet = TRUE)
spec <- malariafwd:::.example_run_spec
out <- file.path(tempdir(), "example-run")
unlink(out, recursive = TRUE)
run <- run_sim(do.call(vivax_params, spec$params), seed = spec$seed, outdir = out)

## every kind of case, recurrence and transmission link the package names occurs in it
ev <- load_events(run)
cs <- classify_cases(ev)
stopifnot(is.na(run$extinct_tick), !run$never_started,
          setequal(cs$origin, c("imported", "local reinfection", "relapse", "recrudescence")),
          any(cs$origin == "relapse" & cs$imported),
          setequal(transmission_links(ev)$kind, c("direct", "liver")),
          length(unique(classify_chain(ev)$who_class)) == 4L,
          any(ev$is_recombinant == 1L))

## the run's files, in name order, as a plain ustar archive (which R's own untar() reads, as
## example_run() does): no owner's name or number, no extended attributes, no macOS metadata
files <- sort(list.files(out, all.files = TRUE, no.. = TRUE), method = "radix")
tarball <- normalizePath(file.path("inst", "extdata", "example_run.tar.xz"), mustWork = FALSE)
dir.create(dirname(tarball), showWarnings = FALSE, recursive = TRUE)
bsd <- any(grepl("bsdtar", system2("tar", "--version", stdout = TRUE)))
owner <- if (bsd) {
  c("--format", "ustar", "--uid", "0", "--gid", "0", "--numeric-owner", "--no-xattrs", "--no-acls",
    "--no-fflags")
} else {
  c("--format=ustar", "--owner=0", "--group=0", "--numeric-owner")
}
if (bsd && Sys.info()[["sysname"]] == "Darwin") owner <- c(owner, "--no-mac-metadata")
Sys.setenv(COPYFILE_DISABLE = "1", XZ_OPT = "-9e")
unlink(tarball)
status <- system2("tar", c(owner, "-cJf", shQuote(tarball), "-C", shQuote(out), shQuote(files)))
stopifnot(status == 0L)
cat(sprintf("%s: %d files, %.0f KB\n", tarball, length(files), file.size(tarball) / 1024))
