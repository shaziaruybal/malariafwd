# Is a SLiM binary available?

The binary is the one the `MALARIAFWD_SLIM` environment variable names
if it is set, otherwise `slim` on the PATH.
[`check_setup()`](https://shaziaruybal.github.io/malariafwd/reference/check_setup.md)
shows the path it found.

## Usage

``` r
slim_available()
```

## Value

`TRUE` if that binary is found.

## Examples

``` r
slim_available()
#> [1] FALSE
```
