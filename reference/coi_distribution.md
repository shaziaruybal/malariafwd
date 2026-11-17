# Distribution of complexity of infection over the infected host-days

The patent lineages (`moi_active`) of every infected host on each day
from `from` to `to`, as counts of host-days and their shares.

## Usage

``` r
coi_distribution(state, from = 0, to = NULL)
```

## Arguments

- state:

  a state table from
  [`load_state()`](https://shaziaruybal.github.io/malariafwd/reference/load_state.md),
  or a run (a `vivax_run` or its directory).

- from, to:

  the days, both included: `from = 0` is day 1, and `to = NULL` is the
  state table's last day.

## Value

a data.table, one row per COI, by COI: `moi_active` (integer), `n` (the
infected host-days, integer) and `frac` (their share, from 0 to 1,
summing to 1).

## Examples

``` r
coi_distribution(example_run(), from = 200)
#>    moi_active     n        frac
#>         <int> <int>       <num>
#> 1:          1  2241 0.744023904
#> 2:          2   589 0.195551129
#> 3:          3   126 0.041832669
#> 4:          4    47 0.015604250
#> 5:          5     9 0.002988048
```
