# Kaplan-Meier curve of the time to first relapse, per hypnozoite batch

The fraction of hypnozoite batches not yet relapsed against the days
since the infection that laid each down, from
[`time_to_relapse()`](https://shaziaruybal.github.io/malariafwd/reference/time_to_relapse.md),
as the Kaplan-Meier estimate.

## Usage

``` r
plot_relapse_km(x, title = NULL)
```

## Arguments

- x:

  a run, or the table
  [`time_to_relapse()`](https://shaziaruybal.github.io/malariafwd/reference/time_to_relapse.md)
  returns.

- title:

  plot title.

## Value

a ggplot object.

## See also

[`time_to_relapse()`](https://shaziaruybal.github.io/malariafwd/reference/time_to_relapse.md)

## Examples

``` r
plot_relapse_km(time_to_relapse(example_run()))
```
