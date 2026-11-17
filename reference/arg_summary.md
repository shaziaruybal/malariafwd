# Structure of the recorded ancestral recombination graph (ARG)

Whole-genealogy summary: how many local trees and edges, how deep it
goes, and – the key diagnostic for a forward run – how many ROOTS each
local tree has. A fully coalesced genealogy has one root; this model
seeds independent founder lineages, so the expected result is a forest
with roughly one root per surviving founder lineage.

## Usage

``` r
arg_summary(run, ids = NULL)
```

## Arguments

- run:

  a run: a `vivax_run` (from
  [`run_sim()`](https://shaziaruybal.github.io/malariafwd/reference/run_sim.md)
  or
  [`load_run()`](https://shaziaruybal.github.io/malariafwd/reference/load_run.md))
  or the path of a run directory.

- ids:

  optional pedigree ids to restrict to (default: all recorded genomes);
  ids that are not genomes in the tree sequence are left out, and none
  being one is an error.

## Value

a one-row data.table: `num_trees`, `num_nodes`, `num_edges`,
`num_samples`, `seq_length_mb`, `roots_min`/`roots_median`/`roots_max`,
`pct_trees_coalesced`, `max_root_time` (ticks before the end of the
run). A run that recorded no genome, such as one that never started,
gives the tree sequence's counts, no samples, 0 roots, and `NA` for
`pct_trees_coalesced` and `max_root_time`.

## Python

The IBD and ARG functions read a run's tree sequence with Python's
tskit, through the reticulate package. Either of two routes works.
*Automatic*: with no Python chosen, reticulate sets up its own
environment with tskit (\>= 1.0) and numpy the first time one of these
functions runs; this needs internet once. *Your own Python*: set
`RETICULATE_PYTHON` to a Python that has those modules (in `~/.Renviron`
to keep it) before the first such call of a session. A function that
cannot reach them stops with an error of class `malariafwd_no_python`
that says what to do, and
[`check_setup()`](https://shaziaruybal.github.io/malariafwd/reference/check_setup.md)
reports what is found.

## Examples

``` r
arg_summary(example_run())
#>    num_trees num_nodes num_edges num_samples seq_length_mb roots_min
#>        <int>     <int>     <int>       <int>         <num>     <int>
#> 1:      4124      1679     12850        1679      24.21467        33
#>    roots_median roots_max pct_trees_coalesced max_root_time
#>           <num>     <int>               <num>         <num>
#> 1:           33        33                   0           399
```
