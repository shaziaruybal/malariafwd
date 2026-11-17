# Load the event log of a run

Reads `events.tsv`. Every run gives the same 21 columns, in the order
and classes the model writes them
([events.tsv](https://shaziaruybal.github.io/malariafwd/reference/malariafwd-outputs.md)
describes them), whatever the run holds: a run with no event gives a
table with no rows. A column the file lacks is `NA`, a column it adds
follows the 21, and a value that does not fit its column stops with the
column named.

## Usage

``` r
load_events(run)
```

## Arguments

- run:

  a run: a `vivax_run` (from
  [`run_sim()`](https://shaziaruybal.github.io/malariafwd/reference/run_sim.md)
  or
  [`load_run()`](https://shaziaruybal.github.io/malariafwd/reference/load_run.md))
  or the path of a run directory.

## Value

a `data.table` with one row per event, in the order they were logged.

## Examples

``` r
ev <- load_events(example_run())
ev[, .N, by = event_type]
#>                   event_type     N
#>                       <char> <int>
#>  1:      vh_inoculation_seed    16
#>  2:         hv_uptake_clonal   797
#>  3:               mosq_death  1206
#>  4:          human_clearance   140
#>  5:           vh_inoculation   206
#>  6:    hv_uptake_recombinant   114
#>  7:    vh_inoculation_import    41
#>  8:          bloodcure_given    83
#>  9:    case_bloodstage_clear   138
#> 10:    hypnozoite_activation   158
#> 11:            radcure_given    22
#> 12: radcure_hypnozoite_clear    29
#> 13:          radcure_failure     5
#> 14:        bloodcure_persist    11
#> 15:      bloodcure_subpatent     5
#> 16:     hypnozoite_clearance    47
#> 17:            recrudescence     5
```
