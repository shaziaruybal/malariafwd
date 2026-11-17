# What each parameter means, its unit and the values it may take

The model's parameters as documented on
[`vivax_params()`](https://shaziaruybal.github.io/malariafwd/reference/vivax_params.md):
one row per parameter, by group.

## Usage

``` r
param_info(params = NULL)
```

## Arguments

- params:

  `NULL` for every parameter, or the names of the ones to give, in that
  order.

## Value

a `data.table` of character columns: `name`, `group`, `default` (as
documented), `type` (`"whole number"`, `"number"`, `"TRUE or FALSE"` or
`"text"`), `unit`, `allowed`, `meaning`, `source` (where the default
came from) and `setting` (how to choose one, for the parameters that are
easy to get wrong); `source` and `setting` are `NA` when none is given.

## Examples

``` r
param_info(c("N_HOSTS", "BITING_RATE"))
#>           name        group default         type                 unit
#>         <char>       <char>  <char>       <char>               <char>
#> 1:     N_HOSTS   Population      40 whole number                hosts
#> 2: BITING_RATE Transmission     0.1       number per mosquito per day
#>                                                                                                                                                          allowed
#>                                                                                                                                                           <char>
#> 1:                                                                                                                                             whole number >= 1
#> 2: number > 0; a mosquito bites at most once a day, so it warns when the peak over the run, with its seasonal, intervention and resurgence multipliers, passes 1
#>                                                                                                                 meaning
#>                                                                                                                  <char>
#> 1:                                                                                      Humans; the number is constant.
#> 2: Daily probability that a mosquito bites. Bites per human per day are BITING_RATE x N_MOSQ / N_HOSTS, not this value.
#>    source
#>    <char>
#> 1:   <NA>
#> 2:   <NA>
#>                                                                                                                                                                                                                      setting
#>                                                                                                                                                                                                                       <char>
#> 1:                                                      prevalence and COI are scale-dependent, so a transmission rate fitted at one population size does not transfer: calibrate at the size you will publish (roadmap S10)
#> 2: the main transmission knob, but prevalence and complexity of infection are scale-dependent, so a value fitted on a small population does not transfer: calibrate it at the population size you will publish (roadmap S10)
```
