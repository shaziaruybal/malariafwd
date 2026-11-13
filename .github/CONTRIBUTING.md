# Contributing to malariafwd

Thank you for your interest in malariafwd. It is a pre-release (0.1.0): its interface and its model
may still change, so please get in touch before building on it.

## Reporting a bug or asking a question

Open an [issue](https://github.com/shaziaruybal/malariafwd/issues) and choose the template that
fits. For a bug, please include:

- the smallest example that shows it, ideally a run with few hosts and few days;
- the output of `malariafwd::check_setup()`, which gives the versions of R, SLiM, the model,
  Python, tskit and numpy that the package finds;
- the output of `sessionInfo()`;
- for a run that did something unexpected, its `params.tsv` and `cmdline.txt`, from which it can
  be repeated (`run_sim(load_params(run), seed = run$seed)`), and its `run.log` if it failed.

## Proposing a change

Please open an issue first, so that we can agree on the change before you spend time on it. Then,
on a branch of your fork:

- add tests for what you change in `tests/testthat/`, and run `devtools::test()` and
  `devtools::check()`;
- document exported functions with roxygen, and run `devtools::document()`;
- add a line to `NEWS.md`;
- for a change to what the model does, add an entry to `DECISIONS.md`, which records the reasons
  for the model's choices, and update `MODEL_DESCRIPTION.md`;
- after a change to the model file or to a default, make the example run again with
  `Rscript data-raw/example_run.R` (it needs SLiM 5.2 on the PATH);
- after a change to the parameter dictionary, update the documentation's tables with
  `Rscript dev/update-doc-tables.R`.

## Code of conduct

Please note that the malariafwd project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this project, you agree to
abide by its terms.
