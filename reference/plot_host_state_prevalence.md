# Host-state prevalence over time (stacked area)

The share of hosts in each state, day by day, stacked; subpatent
carriers and the hosts with hypnozoites alone are bands of their own.

## Usage

``` r
plot_host_state_prevalence(x)
```

## Arguments

- x:

  a run, a state table from
  [`load_state()`](https://shaziaruybal.github.io/malariafwd/reference/load_state.md),
  or the table
  [`host_state_prevalence()`](https://shaziaruybal.github.io/malariafwd/reference/host_state_prevalence.md)
  returns.

## Value

a ggplot object.

## See also

[`host_state_prevalence()`](https://shaziaruybal.github.io/malariafwd/reference/host_state_prevalence.md)

## Examples

``` r
plot_host_state_prevalence(example_run())
```
