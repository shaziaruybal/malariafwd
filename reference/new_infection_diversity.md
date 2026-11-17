# Genetic diversity of new infections, over the run or by window

The relatedness of the new infections of a window of days, those from
bites and imports: their mean pairwise IBD, the share of their pairs
that are clonal and how many lineages they carry. Relapses and
recrudescences are left out: their genomes come from parasites already
in the host. As transmission falls, the clonal share rises and the
lineages fall. With `bin`, each window of the run's grid gets a row.

## Usage

``` r
new_infection_diversity(
  run,
  burn_frac = 0.5,
  from = NULL,
  to = NULL,
  bin = NULL,
  cap = 250,
  clonal = 0.5,
  min_pairs = 15,
  seed = 1
)
```

## Arguments

- run:

  a run: a `vivax_run` (from
  [`run_sim()`](https://shaziaruybal.github.io/malariafwd/reference/run_sim.md)
  or
  [`load_run()`](https://shaziaruybal.github.io/malariafwd/reference/load_run.md))
  or the path of a run directory.

- burn_frac, from, to:

  the window, as
  [`summarise_run()`](https://shaziaruybal.github.io/malariafwd/reference/summarise_run.md)
  takes it: the days after the first `burn_frac` of `SIM_DAYS`, or days
  `from` to `to` (`to = NULL` is `SIM_DAYS`).

- bin:

  `NULL` for one row over the window, or a whole number of days, at
  least 1, for a row for every window of the run's grid (days `0` to
  `bin - 1`, `bin` to `2 * bin - 1`, ..., as in
  [`force_of_infection()`](https://shaziaruybal.github.io/malariafwd/reference/force_of_infection.md))
  that meets it.

- cap:

  the most new infections drawn from a window.

- clonal:

  the IBD at or above which a pair is clonal.

- min_pairs:

  the fewest pairs of genomes a window needs for its measures, a whole
  number, at least 1; 15 pairs are 6 genomes.

- seed:

  seed for the draw when a window holds more than `cap`, the `k`-th
  window drawing with `seed + k - 1`, so repeat calls on the same run
  agree. The caller's RNG stream is left untouched.

## Value

a data.table, one row per window, by window: `win` (its first day on the
grid with `bin`; without, the window's first day), `days` (its days
inside the window asked for), `n` (the genomes compared: those drawn
that the tree sequence holds, or all those drawn when it is not read),
`n_pairs` (`n * (n - 1) / 2`), `mean_ibd`, `clonal_frac` (the share of
the pairs with IBD at least `clonal`, from 0 to 1) and `n_lineages`, the
last three `NA` below `min_pairs`. A window with no new infection has
`n = 0`.

## Details

Each window draws at most `cap` of its new infections. The tree sequence
is read once for every window, and only the pairs within a window are
compared; it is not read at all when no window has the `min_pairs` pairs
its measures need.

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
run <- example_run()
new_infection_diversity(run)                        # the second half of the run
#>      win  days     n n_pairs  mean_ibd clonal_frac n_lineages
#>    <int> <int> <int>   <int>     <num>       <num>      <int>
#> 1:   201   200    96    4560 0.1263772   0.1243421         43
new_infection_diversity(run, bin = 100, cap = 40)   # by 100 days
#>      win  days     n n_pairs  mean_ibd clonal_frac n_lineages
#>    <int> <int> <int>   <int>     <num>       <num>      <int>
#> 1:   200    99    22     231 0.1084360  0.11688312         13
#> 2:   300   100    40     780 0.1023746  0.09102564         26
#> 3:   400     1     0       0        NA          NA         NA
```
