# Pairs of genomes with a known true relationship, and their IBD

Assembles pairs whose relationship is known from the event log – a
genome with itself, clonal founder copies, a relapse against the
hypnozoite batch it activated from, the same lineage in two hosts,
siblings co-transmitted from one oocyst, and distinct founder lineages –
and attaches the ground-truth IBD of each. Feed it to
[`plot_relatedness_ladder()`](https://shaziaruybal.github.io/malariafwd/reference/plot_relatedness_ladder.md)
to check the genealogy against the textbook values.

## Usage

``` r
ibd_truth_pairs(run, max_per_class = 60)
```

## Arguments

- run:

  a run: a `vivax_run` (from
  [`run_sim()`](https://shaziaruybal.github.io/malariafwd/reference/run_sim.md)
  or
  [`load_run()`](https://shaziaruybal.github.io/malariafwd/reference/load_run.md))
  or the path of a run directory.

- max_per_class:

  the most pairs per relationship class, the first the event log gives.

## Value

a data.table of `cls` (a factor of the six classes), `a` and `b`
(integer pedigree ids) and `ibd` (double), one row per pair. A class
with no pair in the run has no rows, and a run with none, such as one
that never started, gives a table with no rows.

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
tp <- ibd_truth_pairs(example_run())
tp[, .(pairs = .N, mean_ibd = mean(ibd)), by = cls]
#>                      cls pairs  mean_ibd
#>                   <fctr> <int>     <num>
#> 1:    self (same genome)    20 1.0000000
#> 2:         clonal copies     4 1.0000000
#> 3: relapse vs its origin    60 1.0000000
#> 4: same lineage, 2 hosts    12 1.0000000
#> 5:      same-oocyst sibs     8 0.3680145
#> 6:     distinct founders    45 0.0000000
```
