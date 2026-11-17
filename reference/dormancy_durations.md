# Hypnozoite dormancy durations, from a batch to each of its relapses

One row per relapse (a `hypnozoite_activation`) whose batch's row is in
the log: the days from the day its hypnozoite batch was laid down to the
day it activated, latency included. An import's batch, laid down a day
after its primary, counts from its own row; founder seeding lays down no
batch. These are the durations of the hypnozoites that activated before
the run ended, so their mean is shorter than the mean time to
activation;
[`time_to_relapse()`](https://shaziaruybal.github.io/malariafwd/reference/time_to_relapse.md)
gives the censored view, one row per batch.

## Usage

``` r
dormancy_durations(events)
```

## Arguments

- events:

  an events table from
  [`load_events()`](https://shaziaruybal.github.io/malariafwd/reference/load_events.md),
  or a run (a `vivax_run` or its directory).

## Value

a data.table, one row per relapse, by `batch_tick`, `batch_pedigree_id`
and `tick`, every column an integer: `human_id`, `batch_pedigree_id`
(the batch it woke from), `pedigree_id` (the relapse's genome),
`batch_tick` (the day the batch was laid down), `tick` (the day it
activated) and `duration` (`tick - batch_tick`, in days).

## See also

[`time_to_relapse()`](https://shaziaruybal.github.io/malariafwd/reference/time_to_relapse.md)

## Examples

``` r
d <- dormancy_durations(example_run())
summary(d$duration)   # days from a batch to each of its relapses
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>    1.00   19.25   42.50   58.35   78.75  282.00 
```
