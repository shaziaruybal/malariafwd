# Load a completed run as a `vivax_run` object

Wraps a run's output directory; use
[`load_events()`](https://shaziaruybal.github.io/malariafwd/reference/load_events.md),
[`load_state()`](https://shaziaruybal.github.io/malariafwd/reference/load_state.md),
[`trees_path()`](https://shaziaruybal.github.io/malariafwd/reference/trees_path.md)
and
[`load_params()`](https://shaziaruybal.github.io/malariafwd/reference/load_params.md)
to access its contents.

## Usage

``` r
load_run(outdir)
```

## Arguments

- outdir:

  a run output directory, or a `vivax_run`, whose directory is read
  again.

## Value

a `vivax_run` object, a list of:

- `outdir`, the run directory (character);

- `params`, the parameter set, as
  [`load_params()`](https://shaziaruybal.github.io/malariafwd/reference/load_params.md)
  gives it;

- `seed`, the seed SLiM ran with (integer; a double beyond the 32-bit
  range, which only a run started with `slim` directly, with `-d SEED`,
  can have);

- `extinct_tick`, the tick on which the run went extinct and ended, or
  `NA` (integer);

- `never_started`, `TRUE` for a run that never started (logical);

- `slim_version` and `malariafwd_version` (character;
  `malariafwd_version` is `"unknown"` for a run started with `slim`
  directly);

- `genome_length`, the simulated genome in bp, and `n_chromosomes`
  (integers), and `recomb_within`, the crossover rate within a
  chromosome per bp per meiosis (double);

- `complete`, `TRUE` when the directory holds the whole run its `DONE`
  file describes (logical).

A field the run's `params.tsv` does not record is `NA`.

## Details

A run goes extinct on the first tick on which no parasite is left in any
host or mosquito, hypnozoites included, and no importation is still to
come. It ends on that tick: its `events.tsv` and `state.tsv` stop there,
and its `params.tsv` ends with the tick as an `EXTINCT_TICK` row, also
when the tick is the run's last day.

A founder-free run whose imports all fail to establish a parasite never
starts, whether by chance or because its settings let no import take
(`P_PRIMARY = 0` with a hypnozoite batch mean of 0, or an
`IMPORT_RATE_PER_TICK` window whose rate is 0 throughout, as under
`IMPORT_RATE_BOOST = 0` from a trigger before the window). It ends on
the first tick after its last import day, but not before tick 10, or on
its last day if that comes first, and its `params.tsv` ends with a
`NEVER_STARTED` row. It holds no parasite: its tree sequence has no
nodes and its events are at most imports that did not take. Every
function that reads a run takes such a run: its tables have their usual
columns and no rows, its counts are 0, and a quantity that needs a
parasite is `NA`
([`arg_summary()`](https://shaziaruybal.github.io/malariafwd/reference/arg_summary.md)
reports no samples).

The functions that summarise a window of a run take its length from
`SIM_DAYS`, so the days after an extinction count as days with no
parasite.

A directory without `params.tsv` or `events.tsv` is not a run and is
refused. One whose `DONE` file is missing, or whose files have changed
since the run finished (see
[`run_sim()`](https://shaziaruybal.github.io/malariafwd/reference/run_sim.md)),
loads with `complete = FALSE` and a warning of class
`malariafwd_incomplete_run`, given once per directory per session by
every function that reads the directory. A run started with `slim`
directly has no `DONE`.
`suppressWarnings(..., classes = "malariafwd_incomplete_run")` silences
it.

## Examples

``` r
run <- load_run(example_run()$outdir)   # a run directory, read again
run$seed
#> [1] 2
run$params$N_HOSTS
#> [1] 60
```
