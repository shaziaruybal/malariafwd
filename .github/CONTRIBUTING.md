# Contributing to malariafwd

Thank you for your interest in malariafwd. It is a pre-release (0.1.0): its interface and its model
may still change, so please get in touch before building on it.

## Reporting a bug or asking a question

Open an [issue](https://github.com/shaziaruybal/malariafwd/issues) and choose the template that
fits. For a bug, please include:

- the smallest example that shows it, ideally a run with few hosts and few days;
- the output of `malariafwd::check_setup()`, which gives the versions of malariafwd, R, SLiM,
  reticulate, Python, tskit and numpy, and where it found SLiM and the model;
- the output of `sessionInfo()`;
- for a run that did something unexpected, its `params.tsv` and `cmdline.txt`, which repeat it
  (with both in a directory `dir`, `run_sim(load_params(dir), seed = s, outdir = tempfile())`, where
  `s` is the run's seed, `ACTUAL_SEED` in `params.tsv`), and its `run.log` if it failed.

## Proposing a change

Please open an issue first, so that we can agree on the change before you spend time on it. Then,
on a branch of your fork:

- add tests for what you change in `tests/testthat/`, and run `devtools::test()` and
  `devtools::check(cran = FALSE)`: with SLiM installed, the `--as-cran` check stops at the first
  example that runs it, on two FIFOs that processx keeps open so that SLiM stops if R dies;
- document exported functions with roxygen, and run `devtools::document()`;
- add a line to `NEWS.md`;
- for a change to what the model does, update `MODEL_DESCRIPTION.md` and give the reason in
  `DECISIONS.md`: change the entry the choice belongs to, or add one to its topic's section with
  the next free number (entries keep their numbers, which the code and the docs cite);
- after a change to the model file or to a default, make the example run again with
  `Rscript data-raw/example_run.R` (it needs SLiM 5.2 on the PATH), run the validation scripts
  with `Rscript dev/validation-results.R`, which rewrites `inst/validation/RESULTS.md`, and remake
  the documents' figures and the numbers they quote with `Rscript dev/make-doc-figures.R`;
- after a change to the parameter dictionary (`R/param-dictionary.R`) or the dictionary of the
  output columns (`R/outputs.R`), update the documentation's tables with
  `Rscript dev/update-doc-tables.R`.

## Code of conduct

Please note that the malariafwd project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this project, you agree to
abide by its terms.
