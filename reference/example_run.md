# A small run that ships with the package

A run of the package's own model, small enough to bundle: 60 hosts and
150 mosquitoes for 400 days, with imports, case management,
recrudescence and co-transmission switched on, so that every kind of
case, recurrence and transmission link the package names occurs in it.
Every function that reads a run takes it, so the examples, the vignettes
and a first look at the package need no SLiM; the relatedness functions
still need Python with tskit
([`tskit_available()`](https://shaziaruybal.github.io/malariafwd/reference/tskit_available.md)).

## Usage

``` r
example_run()
```

## Value

a `vivax_run` (see
[`load_run()`](https://shaziaruybal.github.io/malariafwd/reference/load_run.md)).

## Details

The first call in a session unpacks the run into a temporary directory,
and later calls return the same run, unpacking it again if a file in its
directory is missing or has changed. Its settings are
`load_params(example_run())` and its seed `example_run()$seed`;
`data-raw/example_run.R`, in the package's source repository, makes it.

## See also

[`load_run()`](https://shaziaruybal.github.io/malariafwd/reference/load_run.md),
[`run_sim()`](https://shaziaruybal.github.io/malariafwd/reference/run_sim.md)

## Examples

``` r
run <- example_run()
run
#> <vivax_run> seed 2 
#>   outdir: /tmp/RtmpMJLewE/malariafwd-example-run-60b2e786 
#>   400 days, 60 hosts, 150 mosquitoes
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
#>   SLiM 5.2, malariafwd 0.1.0
table(classify_cases(run)$origin)
#> 
#>          imported local reinfection     recrudescence           relapse 
#>                23               113                 5               158 
```
