# Convert a target realised mean (and variance) into lambda and k

The inverse of
[`uptake_moments()`](https://shaziaruybal.github.io/malariafwd/reference/uptake_moments.md):
give it the mean – and optionally the variance – you want the model to
*produce* (e.g. a mean oocyst count per infected mosquito reported in a
feeding study) and it returns the `LAMBDA_*` / `K_OOCYSTS` values to
set.

## Usage

``` r
uptake_calibrate(mean, var = NULL)
```

## Arguments

- mean:

  target realised (zero-truncated) mean; must be \> 1 and at most 1000,
  the largest mean
  [`vivax_params()`](https://shaziaruybal.github.io/malariafwd/reference/vivax_params.md)
  accepts.

- var:

  optional target realised variance.

## Value

list with `lambda`, `k`, and the `achieved` moments from
[`uptake_moments()`](https://shaziaruybal.github.io/malariafwd/reference/uptake_moments.md).

## Details

If `var` is omitted, or is not larger than the Poisson variance at that
mean, a Poisson (`k = 0`) solution is returned – the negative binomial
can only *add* dispersion. The mean is always matched. The
zero-truncated negative binomial the model draws can be only so
dispersed at a given mean: its variance is largest as `k` approaches 0
(about 11.1 at a mean of 3, 271 at a mean of 10), and `k = 0.001`, the
smallest `k` the model draws, comes within 0.5% of that limit. A target
variance that `k = 0.001` falls more than 1% short of returns
`k = 0.001` with a warning.

## Examples

``` r
uptake_calibrate(2.5)                 # Poisson: lambda giving a realised mean of 2.5
#> $lambda
#> [1] 2.231612
#> 
#> $k
#> [1] 0
#> 
#> $achieved
#>    mean     var 
#> 2.50000 1.82903 
#> 
uptake_calibrate(2.5, var = 6)        # overdispersed: lambda + k
#> $lambda
#> [1] 0.2391864
#> 
#> $k
#> [1] 0.06533696
#> 
#> $achieved
#> mean  var 
#>  2.5  6.0 
#> 
```
