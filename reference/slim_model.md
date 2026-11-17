# Path to the bundled SLiM model template

The single source of truth for the model;
[`run_sim()`](https://shaziaruybal.github.io/malariafwd/reference/run_sim.md)
snapshots it into each run's output directory for provenance.

## Usage

``` r
slim_model()
```

## Value

file path to `inst/slim/malariafwd_vivax.slim`.

## Examples

``` r
slim_model()
#> [1] "/home/runner/work/_temp/Library/malariafwd/slim/malariafwd_vivax.slim"
file.exists(slim_model())
#> [1] TRUE
```
