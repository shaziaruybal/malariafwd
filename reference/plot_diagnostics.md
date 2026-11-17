# One-glance run-health dashboard

A six-panel diagnostic for a single run: prevalence, incidence by
origin, the force of blood-stage infection by origin (FOI vs molFOB; see
[`force_of_infection()`](https://shaziaruybal.github.io/malariafwd/reference/force_of_infection.md)),
mean COI, the hypnozoite reservoir, and circulating-lineage richness –
all over time, with the burn-in marked. Surfaces the common failure
modes at a glance: extinction (lines crash to zero), saturation
(prevalence pinned near 100%), and non-stationarity (metrics still
trending at the burn-in). The two "by origin" panels share a palette.
The force panel stacks the new clones: imported and local reinfection
(FOI), and relapse, which molFOB adds (Koepfli et al. 2013). A
recrudescence, a new episode but not a new clone, appears in the
incidence panel only. The subtitle gives
[`summarise_run()`](https://shaziaruybal.github.io/malariafwd/reference/summarise_run.md)'s
numbers for the days after the burn-in.

## Usage

``` r
plot_diagnostics(run, bin = 30, burn_frac = 0.5)
```

## Arguments

- run:

  a run: a `vivax_run` (from
  [`run_sim()`](https://shaziaruybal.github.io/malariafwd/reference/run_sim.md)
  or
  [`load_run()`](https://shaziaruybal.github.io/malariafwd/reference/load_run.md))
  or the path of a run directory.

- bin:

  window width in days.

- burn_frac:

  the share of the run's `SIM_DAYS` days marked as burn-in; the subtitle
  gives
  [`summarise_run()`](https://shaziaruybal.github.io/malariafwd/reference/summarise_run.md)'s
  numbers for the rest.

## Value

a patchwork object (combined ggplot).

## Details

The panels run to the end of the run's last whole window of `bin` days:
a part-window after it, which the panels that count per window would
draw as a dip, is left out. A run that goes extinct shows lines that
drop to zero and stay there, and the subtitle says on which day it went
extinct. A window with no case is drawn as zero.

## References

Koepfli C, et al. (2013) PLoS Negl Trop Dis 7(9):e2403 (molFOB).

## Examples

``` r
plot_diagnostics(example_run())
```
