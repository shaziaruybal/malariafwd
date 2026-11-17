# Force of blood-stage infection: molFOB vs FOI

Plots the compound force of blood-stage infection (**molFOB** – every
new blood-stage clone, relapses included, recrudescences excluded)
against the force of infection (**FOI** – the new clones from bites and
imports). The distance between the lines is the relapse contribution,
which in *P. vivax* is typically the larger share and is invisible to
any method that equates new episodes with new transmission.

## Usage

``` r
plot_force_of_infection(x, smooth = 0.4)
```

## Arguments

- x:

  a run, or the table
  [`force_of_infection()`](https://shaziaruybal.github.io/malariafwd/reference/force_of_infection.md)
  returns (for a run, with its defaults).

- smooth:

  span for a loess trend through each series; `NA` draws the raw series,
  as does a span that loess cannot fit to the run's windows without a
  warning (too few windows for it).

## Value

a ggplot object.

## Details

Naming follows Koepfli et al. (2013): molFOB counts new blood-stage
clones however they arise, FOI those a bite or an import brought. A
recrudescence is a new episode but not a new clone, so it is in neither
series. Both series are drawn in every window of the run, as zero where
there was no case.

## See also

[`force_of_infection()`](https://shaziaruybal.github.io/malariafwd/reference/force_of_infection.md)

## Examples

``` r
run <- example_run()
plot_force_of_infection(run)

plot_force_of_infection(force_of_infection(run, bin = 15))

plot_force_of_infection(force_of_infection(run, bin = 60))   # 7 windows: no trend
```
