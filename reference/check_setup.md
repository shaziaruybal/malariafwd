# Check what malariafwd needs: R, SLiM, the model, and Python with tskit

Reports each component the package uses, with its version and where it
was found, and what to do about any that is missing. It never errors.
SLiM and the model are needed to run; reticulate, Python, tskit and
numpy for the IBD and ARG functions (see
[`tskit_available()`](https://shaziaruybal.github.io/malariafwd/reference/tskit_available.md));
ggplot2 and patchwork for the plots.

## Usage

``` r
check_setup(python = TRUE)
```

## Arguments

- python:

  if `FALSE`, Python is not started and its rows are not checked.

## Value

a data.frame of class `malariafwd_setup`, one row per component, with
columns `component`, `ok` (`NA` when not checked), `version`, `detail`
and `advice` (`NA` when there is nothing to do).

## Examples

``` r
check_setup(python = FALSE)
#>   [ok] malariafwd  0.1.0      /home/runner/work/_temp/Library
#>   [ok] R           4.6.2      /opt/R/4.6.2/lib/R
#>   [--] SLiM                   slim
#>        install SLiM >= 5.2 and put it on PATH, or set MALARIAFWD_SLIM=/path/to/slim in ~/.Renviron (RStudio on macOS may not see the shell's PATH)
#>   [ok] model                  /home/runner/work/_temp/Library/malariafwd/slim/malariafwd_vivax.slim
#>   [ok] reticulate  1.47.0     for the IBD and ARG functions
#>   [  ] Python                 not checked
#>   [  ] tskit                  not checked
#>   [  ] numpy                  not checked
#>   [ok] ggplot2     4.0.4      optional, for the plots
#>   [ok] patchwork   1.3.2      optional, for the plots
#> 1 problem: see the advice above.
```
