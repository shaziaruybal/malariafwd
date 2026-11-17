# Host-state prevalence and the hypnozoite reservoir over time

The share of hosts in each of five states on each day of the state
table: no parasite; patent blood-stage clones, without or with
hypnozoites; subpatent clones and no patent one, with or without
hypnozoites; and hypnozoites alone, a reservoir no blood survey can see.
Each day also carries its mean number of hypnozoites per host.

## Usage

``` r
host_state_prevalence(state)
```

## Arguments

- state:

  a state table from
  [`load_state()`](https://shaziaruybal.github.io/malariafwd/reference/load_state.md),
  or a run (a `vivax_run` or its directory).

## Value

a data.table, one row per day of the state table and state, by day and
state: `tick`, `cls` (a factor with levels `"uninfected"`,
`"active only"`, `"active + reservoir"`, `"subpatent"` and
`"reservoir only (invisible)"`), `n` (hosts, 0 when none), `frac` (`n`
over the day's hosts) and `hypnozoites_per_host` (the day's mean
`hyp_reservoir_size` over its hosts, the same on each of its rows, so
`unique(x[, .(tick, hypnozoites_per_host)])` is the reservoir's
trajectory). A state table without `moi_subpatent` counts no subpatent
host.

## See also

[`plot_host_state_prevalence()`](https://shaziaruybal.github.io/malariafwd/reference/plot_host_state_prevalence.md)

## Examples

``` r
hs <- host_state_prevalence(example_run())
hs[tick == 400]   # the last day
#>     tick                        cls     n       frac hypnozoites_per_host
#>    <int>                     <fctr> <int>      <num>                <num>
#> 1:   400                 uninfected    24 0.40000000                  2.3
#> 2:   400                active only     1 0.01666667                  2.3
#> 3:   400         active + reservoir    16 0.26666667                  2.3
#> 4:   400                  subpatent     0 0.00000000                  2.3
#> 5:   400 reservoir only (invisible)    19 0.31666667                  2.3
```
