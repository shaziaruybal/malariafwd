# Run a grid of parameter sets x seeds in parallel

Runs every (parameter set, seed) cell, each in its own directory
`<set>_seed<seed>` under `outdir_base`, in parallel with the current
`future` plan. Set `workers` for a quick local multisession; for a
cluster, set `future::plan(...)` yourself before calling and leave
`workers = NULL`, so that the grid runs under your plan.

## Usage

``` r
run_grid(
  param_sets,
  seeds = 1L,
  outdir_base,
  model = slim_model(),
  workers = NULL,
  verbose = NULL,
  skip_existing = TRUE,
  overwrite = FALSE,
  ...
)
```

## Arguments

- param_sets:

  a
  [`vivax_params()`](https://shaziaruybal.github.io/malariafwd/reference/vivax_params.md)
  object, or a list of them (or of lists of overrides). Its names name
  the cells' directories: they must differ, also ignoring case, and be
  usable as a directory name; blank names become `set01`, `set02`, ...
  by position.

- seeds:

  whole-number seeds to run for each set, none repeated.

- outdir_base:

  the grid's directory, created if needed.

- model:

  path to the SLiM model.

- workers:

  if set, run under `future::plan(multisession, workers)` for this call
  (the previous plan is restored afterward). Set it to the cores you can
  spare. If `NULL`, use the current plan.

- verbose:

  passed to
  [`run_sim()`](https://shaziaruybal.github.io/malariafwd/reference/run_sim.md).

- skip_existing:

  if `TRUE` (the default), a cell that already holds its run is not run
  again.

- overwrite:

  if `TRUE`, a cell that holds another run, or one that did not finish,
  is run again and its run replaced once the new one has finished. Files
  that are not part of a run are never replaced.

- ...:

  passed to
  [`run_sim()`](https://shaziaruybal.github.io/malariafwd/reference/run_sim.md).

## Value

a `data.table` with one row per cell, by set then seed: `set`, `seed`,
`outdir`, `ok` (the cell holds its run, run now or before), `skipped`
(it was not run again), `error` (the failure message, or `NA`),
`extinct_tick` (the tick on which the run went extinct and ended, `NA`
if it did not) and `never_started` (`TRUE` for a founder-free run whose
imports never established a parasite; see
[`load_run()`](https://shaziaruybal.github.io/malariafwd/reference/load_run.md)),
both `NA` for a failed cell, `elapsed_s`, `finished` (UTC),
`slim_version` and `malariafwd_version`, from the run's `DONE`, and then
one column per parameter whose value differs between the sets, of the
parameter's type. A cell that never started is `ok`: not starting is an
outcome of the grid, which `never_started` records. Its run holds no
parasite, so its event and genetic summaries are empty or `NA`: keep the
cells with `ok & !never_started` for those, and count the cells that
never started as an outcome. The cells that went extinct, those that
never started, and those not run again are counted in one message each.
A failed cell does not abort the grid – its directory is left as it was
and its files and log are in `<cell>.failed` – but it warns, and the
reason travels with the result. A warning given in the cells (such as a
parameter set's own warnings, which
[`run_sim()`](https://shaziaruybal.github.io/malariafwd/reference/run_sim.md)
repeats for each run) is given once, after the warning that lists the
failed cells, with the number of cells it came from. If the grid stops
part-way with an error, that error's message lists the warnings of the
cells whose results had reached this session; with workers, a cell's
warnings travel with its result. A cell whose worker dies is run again
on a new worker, its first attempt's directory left beside it as
`<cell>.partial-...`; if that worker dies too, the grid stops with an
error that names the cell. Under `options(warn = 2)` a warning in a cell
fails that cell, in any plan and also inside
[`suppressWarnings()`](https://rdrr.io/r/base/warning.html), and the
warning that lists the failed cells is then an error, so nothing is
returned; the cells that ran are on disk. The caller's random-number
stream is left as it was.

## Details

Every set is checked, and every cell's directory, before any cell
starts: a set
[`vivax_params()`](https://shaziaruybal.github.io/malariafwd/reference/vivax_params.md)
refuses, or a cell whose directory
[`run_sim()`](https://shaziaruybal.github.io/malariafwd/reference/run_sim.md)
would refuse, stops the grid with all of them listed. A cell that
already holds its run (the same parameters, seed and model) is not run
again, so calling `run_grid()` again with the same sets and seeds
resumes a grid that stopped part-way, and adding seeds or sets runs only
the new cells. The cells of an edited set hold another run and stop the
grid unless `overwrite = TRUE`.

`grid.tsv` in `outdir_base`, written when the grid finishes (also when
cells failed), records each cell: its outcome, versions, the key of its
run and every parameter.
[`load_grid()`](https://shaziaruybal.github.io/malariafwd/reference/load_grid.md)
reads it. Cells run by other means, such as one
[`run_sim()`](https://shaziaruybal.github.io/malariafwd/reference/run_sim.md)
per job of a job array, are recorded by calling `run_grid()` with the
same sets and seeds afterwards, which runs nothing. Two grids writing
into one `outdir_base` at the same time race on `grid.tsv`; a later call
rebuilds it.

## Examples

``` r
if (FALSE) { # slim_available()
sets <- list(low = vivax_params(BITING_RATE = 0.05, SIM_DAYS = 100),
             high = vivax_params(BITING_RATE = 0.12, SIM_DAYS = 100))
out <- tempfile("sweep")
g <- run_grid(sets, seeds = 1:2, outdir_base = out)   # workers = 2 for two at a time
g
load_grid(out)
}
```
