# Realised mean and variance of a zero-truncated count distribution

`LAMBDA_OOCYSTS` / `LAMBDA_SPOROZOITES` are **pre-truncation** rates, so
they are *not* the mean count you would read off a paper: conditioning
away the zeros pulls the realised mean up. This returns what the model
actually produces for a given `lambda` (and dispersion `k`), so you can
check a parameter choice against a reported figure.

## Usage

``` r
uptake_moments(lambda, k = 0)
```

## Arguments

- lambda:

  pre-truncation mean (\> 0).

- k:

  dispersion; `0` (default) = Poisson, `> 0` = negative binomial.

## Value

named numeric of `mean` and `var` (the realised, zero-truncated values).

## Details

With base mean \\\lambda\\, dispersion \\k\\ and zero-probability
\\p_0\\, the truncated moments are \\E\[X\] = \lambda/(1-p_0)\\ and
\\Var\[X\] = E\[X^2\]/(1-p_0) - E\[X\]^2\\. `k <= 0` is Poisson (\\p_0 =
e^{-\lambda}\\); `k > 0` is negative binomial with \\Var = \lambda +
\lambda^2/k\\ and \\p_0 = (k/(k+\lambda))^k\\.

## See also

[`uptake_calibrate()`](https://shaziaruybal.github.io/malariafwd/reference/uptake_calibrate.md)
for the inverse (target moments -\> lambda, k).

## Examples

``` r
uptake_moments(2)          # Poisson lambda = 2 -> realised mean ~2.31
#>     mean      var 
#> 2.313035 1.588974 
uptake_moments(2, k = 1)   # negative binomial, same lambda -> realised mean ~3.0, more spread
#> mean  var 
#>    3    6 
```
