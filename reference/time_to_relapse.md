# Time from an infection to the first relapse of its hypnozoite batch

One row per hypnozoite batch: the hypnozoites one genotype of a bite or
an import laid down in a host, so a bite that delivers two genotypes can
lay down two. A batch's clock starts on the day it was established, the
day of the bite, or for an import that also brought a blood-stage
infection, the day after it. Founder hosts carry no batch until a bite
or an import gives them one.

## Usage

``` r
time_to_relapse(run, emptied = c("censor", "follow"), exclude_masked = FALSE)
```

## Arguments

- run:

  a run: a `vivax_run` (from
  [`run_sim()`](https://shaziaruybal.github.io/malariafwd/reference/run_sim.md)
  or
  [`load_run()`](https://shaziaruybal.github.io/malariafwd/reference/load_run.md))
  or the path of a run directory.

- emptied:

  `"censor"` (stop the clock on the day a batch empties without
  relapsing) or `"follow"` (follow it to `SIM_DAYS`).

- exclude_masked:

  if `TRUE`, the first relapse is the first that is not masked.

## Value

a data.table, one row per batch, by `batch_tick` and
`batch_pedigree_id`: `batch_pedigree_id`, `human_id`, `lineage_id`,
`imported` (an imported batch), `batch_tick` (the day it was
established), `n_hypnozoites` (its size then), `time` (days from
`batch_tick` to the first relapse or to censoring), `event` (1 relapse,
0 censored) and `end_reason` (`"relapse"`, `"emptied"` or `"run_end"`).
A run with no batch gives no row.

## Details

The clock stops at the batch's first relapse (`event = 1`). Without one
(`event = 0`) it stops on the day the batch's last hypnozoite left it,
by death or radical cure (`emptied = "censor"`), or at the run's last
day, `SIM_DAYS`, for a batch that still held hypnozoites then, and with
`emptied = "follow"` for one that emptied without relapsing too.
Censoring at emptying gives the time to relapse of the hypnozoites that
live; following to the end gives the curve a trial reports, where
radical cure shows as relapses prevented.

With `exclude_masked = TRUE` a masked relapse, into a patent clone of
its own lineage (see
[`classify_cases()`](https://shaziaruybal.github.io/malariafwd/reference/classify_cases.md)),
does not stop the clock, although its hypnozoites leave the batch: the
first relapse is the first a study could see as a recurrence.

## See also

[`plot_relapse_km()`](https://shaziaruybal.github.io/malariafwd/reference/plot_relapse_km.md),
[`dormancy_durations()`](https://shaziaruybal.github.io/malariafwd/reference/dormancy_durations.md)

## Examples

``` r
tt <- time_to_relapse(example_run())
head(tt)
#>    batch_pedigree_id human_id lineage_id imported batch_tick n_hypnozoites
#>                <int>    <int>     <char>   <lgcl>      <int>         <int>
#> 1:                52       10          7    FALSE         13             5
#> 2:                99       40          4    FALSE         28             7
#> 3:               128       10        127     TRUE         34             1
#> 4:               154       50          8    FALSE         43             1
#> 5:               173       14        172     TRUE         54             1
#> 6:               185       27        184     TRUE         63             5
#>     time event end_reason
#>    <int> <int>     <char>
#> 1:    20     1    relapse
#> 2:    22     1    relapse
#> 3:   128     0    emptied
#> 4:     4     0    emptied
#> 5:    84     1    relapse
#> 6:    20     1    relapse
```
