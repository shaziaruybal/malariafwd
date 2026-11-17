# Run one simulation

Launches SLiM with the given parameters and writes the run into
`outdir`: `events.tsv`, `state.tsv`, `plasmo.trees` and `params.tsv`,
with `model.slim`, `cmdline.txt`, `run.log` and `DONE` (see *The run
directory*). `params.tsv` records every parameter, floats to the 15
significant digits passed to SLiM. The launch is shell-free (processx),
so it handles the SLiM quoting pitfalls for you.

## Usage

``` r
run_sim(
  params = vivax_params(),
  seed = 1L,
  outdir,
  model = slim_model(),
  verbose = NULL,
  timeout = Inf,
  echo = FALSE,
  overwrite = FALSE,
  skip_existing = FALSE
)
```

## Arguments

- params:

  a
  [`vivax_params()`](https://shaziaruybal.github.io/malariafwd/reference/vivax_params.md)
  object, or a named list of overrides. Either is checked by
  [`vivax_params()`](https://shaziaruybal.github.io/malariafwd/reference/vivax_params.md)
  again before the run, so an object edited after it was built is
  checked too.

- seed:

  integer RNG seed. SLiM starts from it (`-s`) and the model sets it as
  `SEED`, so the outputs and the tree-sequence provenance both record
  it.

- outdir:

  the run's directory, created when the run has finished (its parent is
  created if needed). A directory that already holds a run is not
  replaced unless `overwrite = TRUE`, and a directory holding other
  files is never used; these refusals are errors of class
  `malariafwd_outdir_in_use`.

- model:

  path to the SLiM model (defaults to the bundled template).

- verbose:

  if `TRUE`, run the model in VERBOSE mode (full per-event narration).
  `NULL` (default) leaves the param as set (quiet).

- timeout:

  seconds before the run is stopped, with an error saying so (`Inf` = no
  limit).

- echo:

  if `TRUE`, stream SLiM output live to the console.

- overwrite:

  if `TRUE`, a run already in `outdir` is replaced, once the new run has
  finished. Only run files are ever deleted.

- skip_existing:

  if `TRUE` and `outdir` already holds this run (the same parameters,
  seed and model), it is returned without running SLiM.

## Value

invisibly, a `vivax_run` object: the run directory, its parameter set,
seed, outcome, versions, genome and `complete` (see
[`load_run()`](https://shaziaruybal.github.io/malariafwd/reference/load_run.md)).
A run that goes extinct, with no parasite left and no importation still
to come, ends on that tick and gives a message; so does a founder-free
run whose imports never established a parasite.

## The run directory

The run is written into a new directory beside `outdir` and moved to
`outdir` only when it has finished, so `outdir` holds a whole run or
what it held before, never a mix. `DONE`, written last, records when the
run finished, how long it took, its seed, a key of its parameters and
seed, the md5 of the model, the versions, and the size of every file;
[`load_run()`](https://shaziaruybal.github.io/malariafwd/reference/load_run.md)
reports a directory whose `DONE` is missing or whose files have changed
since. SLiM runs in the run's directory (`OUTDIR` is `'.'`) on the copy
of the model it holds, `model.slim`, so the directory is self-contained:
`cmdline.txt` gives the arguments, and running `slim` with them in a
directory holding `model.slim` repeats the run. A custom model that
reads other files must use absolute paths. A run that fails, times out
or is interrupted leaves its files and `run.log` in `<outdir>.failed`;
if R itself is killed, a `<outdir>.partial-*` directory is left beside
`outdir`, to delete by hand.

## Examples

``` r
if (FALSE) { # slim_available()
run <- run_sim(vivax_params(SIM_DAYS = 100), seed = 1, outdir = tempfile("demo"))
run
summarise_run(run)
}
```
