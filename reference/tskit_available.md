# Can the tree-sequence functions reach Python with tskit?

`TRUE` when Python with tskit (\>= 1.0) and numpy can be used through
reticulate. It never errors or warns, so it can guard examples, vignette
chunks and tests. It may start Python, and with no Python chosen,
reticulate then sets up its own environment.

## Usage

``` r
tskit_available()
```

## Value

`TRUE` or `FALSE`.

## Python

The IBD and ARG functions read a run's tree sequence with Python's
tskit, through the reticulate package. Either of two routes works.
*Automatic*: with no Python chosen, reticulate sets up its own
environment with tskit (\>= 1.0) and numpy the first time one of these
functions runs; this needs internet once. *Your own Python*: set
`RETICULATE_PYTHON` to a Python that has those modules (in `~/.Renviron`
to keep it) before the first such call of a session. A function that
cannot reach them stops with an error of class `malariafwd_no_python`
that says what to do, and
[`check_setup()`](https://shaziaruybal.github.io/malariafwd/reference/check_setup.md)
reports what is found.

## Examples

``` r
tskit_available()
#> [1] TRUE
```
