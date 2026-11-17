# Path to a run's tree-sequence file

Path to a run's tree-sequence file

## Usage

``` r
trees_path(run)
```

## Arguments

- run:

  a run: a `vivax_run` (from
  [`run_sim()`](https://shaziaruybal.github.io/malariafwd/reference/run_sim.md)
  or
  [`load_run()`](https://shaziaruybal.github.io/malariafwd/reference/load_run.md))
  or the path of a run directory.

## Value

the path to `plasmo.trees`.

## See also

[plasmo.trees](https://shaziaruybal.github.io/malariafwd/reference/malariafwd-outputs.md)

## Examples

``` r
trees_path(example_run())
#> [1] "/tmp/RtmpjMBTAb/malariafwd-example-run-60b2e786/plasmo.trees"
file.exists(trees_path(example_run()))
#> [1] TRUE
```
