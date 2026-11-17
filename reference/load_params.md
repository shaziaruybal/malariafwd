# The parameter set a run was made from

Reads `params.tsv`
([params.tsv](https://shaziaruybal.github.io/malariafwd/reference/malariafwd-outputs.md))
into a
[`vivax_params()`](https://shaziaruybal.github.io/malariafwd/reference/vivax_params.md)
object holding every parameter, in order and of its type:
`run_sim(load_params(run), seed = run$seed, outdir = tempfile())` gives
the same run again. Values are not checked here
([`run_sim()`](https://shaziaruybal.github.io/malariafwd/reference/run_sim.md)
checks them). The seed, the versions, the derived genome and the outcome
rows are on the
[`load_run()`](https://shaziaruybal.github.io/malariafwd/reference/load_run.md)
object instead.

## Usage

``` r
load_params(run)
```

## Arguments

- run:

  a run: a `vivax_run` (from
  [`run_sim()`](https://shaziaruybal.github.io/malariafwd/reference/run_sim.md)
  or
  [`load_run()`](https://shaziaruybal.github.io/malariafwd/reference/load_run.md))
  or the path of a run directory.

## Value

a
[`vivax_params()`](https://shaziaruybal.github.io/malariafwd/reference/vivax_params.md)
object.

## Details

A file from another version of the package can lack some parameters,
which take their defaults, or hold rows that are not parameters, which
are left out. Either warns, once per directory per session, with class
`malariafwd_params_mismatch`.

## Examples

``` r
p <- load_params(example_run())
p
#> <vivax_params> 63 parameters
#>   non-default:
#>     SIM_DAYS                 = 400
#>     N_HOSTS                  = 60
#>     N_FOUNDER_HOSTS          = 16
#>     N_FOUNDER_LINEAGES       = 10
#>     N_MOSQ                   = 150
#>     LAMBDA_OOCYSTS           = 1
#>     LAMBDA_SPOROZOITES       = 1.5
#>     RADCURE_COVERAGE         = 0.3
#>     BLOODCURE_COVERAGE       = 0.5
#>     P_BLOODCURE_SURVIVAL     = 0.1
#>     P_RECRUD_GIVEN_FAIL      = 0.5
#>     P_DETECT_RELAPSE         = 0.5
#>     IMPORT_START             = 1
#>     IMPORT_RATE_PER_TICK     = 0.05
if (FALSE) { # slim_available()
run <- run_sim(p, seed = 3, outdir = tempfile("seed3"))   # the same settings, another seed
}
```
