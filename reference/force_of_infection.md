# Incidence and force of blood-stage infection, by origin

The cases of each window of the run, by origin, as counts (`n`, the
incidence) and per host-year (`rate`). Among the new clones,
`imported` + `local reinfection` give the mosquito-borne **force of
infection (FOI)**; adding `relapse` gives the **molecular force of
blood-stage infection (molFOB)** – a compound measure that is elevated
in *P. vivax* because relapses contribute new blood-stage clones
(Koepfli et al. 2013). With `detected_only = TRUE` only the cases the
observation layer reported are counted: the incidence, FOI and molFOB of
a line list. Relapses are split into `relapse (imported)` and
`relapse (local)` – a relapse of an imported genotype is a distinct
truth category. Set `split_imported = FALSE` for a single `relapse`
class. Either way `origin` comes back as a factor with **all** levels
present, so the shape of the result depends on the argument and never on
what the run happened to contain: pooling a grid cannot split one
quantity across three columns.

## Usage

``` r
force_of_infection(
  run,
  bin = 30,
  split_imported = TRUE,
  exclude_masked = FALSE,
  detected_only = FALSE
)
```

## Arguments

- run:

  a run: a `vivax_run` (from
  [`run_sim()`](https://shaziaruybal.github.io/malariafwd/reference/run_sim.md)
  or
  [`load_run()`](https://shaziaruybal.github.io/malariafwd/reference/load_run.md))
  or the path of a run directory.

- bin:

  window width in days.

- split_imported:

  split relapses by the provenance of the genotype. `TRUE` (default)
  gives five origin levels, `FALSE` gives four.

- exclude_masked:

  if `TRUE`, leave out the masked relapses, into a patent clone of their
  own lineage (see
  [`classify_cases()`](https://shaziaruybal.github.io/malariafwd/reference/classify_cases.md)):
  molFOB as a genotyping study measures it.

- detected_only:

  if `TRUE`, count only the reported cases (`detected == 1`). Detection
  is drawn once per bite, import, activation or recrudescence
  (`P_DETECT_PRIMARY`, `P_DETECT_RELAPSE`, `P_DETECT_RECRUD`).

## Value

a data.table with one row for every window of the run and every origin,
by window: `win` (the window's first day, an integer: 0, `bin`,
`2 * bin`, ...), `origin`, `n` (the cases, clone onsets as rows of
[`classify_cases()`](https://shaziaruybal.github.io/malariafwd/reference/classify_cases.md),
0 where there are none;
[`episodes()`](https://shaziaruybal.github.io/malariafwd/reference/episodes.md)
groups them by host and day), `days` (the days the window covers within
the run, which is `SIM_DAYS` days long whether or not it went extinct:
`bin - 1` for the first window, which starts on day 1, and fewer for the
last when the run ends inside it), `rate` (per host-year, scaled by
`days`) and `is_acquisition`. Sum only the `is_acquisition` rows for FOI
or molFOB: a recrudescence is returned as its own row so incidence stays
complete, but it is not a newly acquired clone.

## Details

A **recrudescence** is returned as its own origin with
`is_acquisition = FALSE`. It is a new patent episode but not a new
clone, so including it would inflate molFOB – which counts clone
acquisitions. Worth noting that a real molFOB, measured by comparing
genotypes between study visits, would often count a recrudescence as a
new clone, because the clone was undetectable at the previous visit; the
gap between the two sums here estimates that bias.

## References

Koepfli C, et al. (2013) A high force of *Plasmodium vivax* blood-stage
infection drives the rapid acquisition of immunity in Papua New Guinean
children. PLoS Negl Trop Dis 7(9):e2403.

## See also

[`plot_force_of_infection()`](https://shaziaruybal.github.io/malariafwd/reference/plot_force_of_infection.md)

## Examples

``` r
force_of_infection(example_run(), bin = 100)
#> Key: <win>
#>       win             origin     n  days      rate is_acquisition
#>     <int>             <fctr> <int> <int>     <num>         <lgcl>
#>  1:     0           imported     7    99 0.4301347           TRUE
#>  2:     0  local reinfection    11    99 0.6759259           TRUE
#>  3:     0 relapse (imported)     2    99 0.1228956           TRUE
#>  4:     0    relapse (local)     5    99 0.3072391           TRUE
#>  5:     0      recrudescence     0    99 0.0000000          FALSE
#>  6:   100           imported     3   100 0.1825000           TRUE
#>  7:   100  local reinfection    18   100 1.0950000           TRUE
#>  8:   100 relapse (imported)     3   100 0.1825000           TRUE
#>  9:   100    relapse (local)    23   100 1.3991667           TRUE
#> 10:   100      recrudescence     3   100 0.1825000          FALSE
#> 11:   200           imported     6   100 0.3650000           TRUE
#> 12:   200  local reinfection    17   100 1.0341667           TRUE
#> 13:   200 relapse (imported)     8   100 0.4866667           TRUE
#> 14:   200    relapse (local)    32   100 1.9466667           TRUE
#> 15:   200      recrudescence     0   100 0.0000000          FALSE
#> 16:   300           imported     7   100 0.4258333           TRUE
#> 17:   300  local reinfection    67   100 4.0758333           TRUE
#> 18:   300 relapse (imported)    13   100 0.7908333           TRUE
#> 19:   300    relapse (local)    72   100 4.3800000           TRUE
#> 20:   300      recrudescence     2   100 0.1216667          FALSE
#> 21:   400           imported     0     1 0.0000000           TRUE
#> 22:   400  local reinfection     0     1 0.0000000           TRUE
#> 23:   400 relapse (imported)     0     1 0.0000000           TRUE
#> 24:   400    relapse (local)     0     1 0.0000000           TRUE
#> 25:   400      recrudescence     0     1 0.0000000          FALSE
#>       win             origin     n  days      rate is_acquisition
#>     <int>             <fctr> <int> <int>     <num>         <lgcl>
```
