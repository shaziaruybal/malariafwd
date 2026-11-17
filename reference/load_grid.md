# Load the record of a grid

Reads the `grid.tsv` that
[`run_grid()`](https://shaziaruybal.github.io/malariafwd/reference/run_grid.md)
writes in its `outdir_base`: one row per cell, with
[`run_grid()`](https://shaziaruybal.github.io/malariafwd/reference/run_grid.md)'s
columns. `outdir` is resolved against `outdir_base`, so a grid that was
moved or copied still loads. It warns, naming up to five, when a cell
recorded as `ok` no longer holds its whole run, and says so when the
cells were run by different versions of SLiM or malariafwd.

## Usage

``` r
load_grid(outdir_base, params = c("varying", "all"))
```

## Arguments

- outdir_base:

  the grid's directory.

- params:

  `"varying"` (the default) for the parameters whose values differ
  between the cells, `"all"` for every parameter, each of its type.

## Value

a `data.table`, as
[`run_grid()`](https://shaziaruybal.github.io/malariafwd/reference/run_grid.md)
returns.

## Examples

``` r
if (FALSE) { # slim_available()
out <- tempfile("sweep")
run_grid(vivax_params(SIM_DAYS = 100), seeds = 1:2, outdir_base = out)
load_grid(out)
}
```
