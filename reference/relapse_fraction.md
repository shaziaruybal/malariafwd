# Relapse fraction, by episodes and by acquisitions

Relapses as a share of a run's blood-stage cases, in two readings, both
counted per clone onset (a row of
[`classify_cases()`](https://shaziaruybal.github.io/malariafwd/reference/classify_cases.md),
so a bite that delivers two genotypes counts twice):

## Usage

``` r
relapse_fraction(events, exclude_masked = FALSE)
```

## Arguments

- events:

  an events table from
  [`load_events()`](https://shaziaruybal.github.io/malariafwd/reference/load_events.md),
  or a run (a `vivax_run` or its directory).

- exclude_masked:

  if `TRUE`, leave out the masked relapses.

## Value

a one-row data.table: `n_relapse`, `n_new_infection` (imported and local
reinfection), `n_recrudescence`, `relapse_frac_episodes` and
`relapse_frac_acquisitions`, each fraction `NA` when its denominator is
0 (a run with no case).

## Details

- `relapse_frac_episodes` – relapses over **every** case: imported,
  local reinfection, relapse and recrudescence.

- `relapse_frac_acquisitions` – relapses over the **acquisitions**:
  imported, local reinfection and relapse. A recrudescence is a new
  clinical episode but not a newly acquired parasite (its genotype was
  already in the host), so it is left out. This is the fraction the
  calibration in section 13 of the model description compares with
  published studies.

The two are equal when a run has no recrudescence, as at the default
`P_RECRUD_GIVEN_FAIL = 0`. Both count each host's first case, which is
not a recurrence;
[`recurrences()`](https://shaziaruybal.github.io/malariafwd/reference/recurrences.md)
gives the recurrences alone.

Masked relapses, into a patent clone of their own lineage (see
[`classify_cases()`](https://shaziaruybal.github.io/malariafwd/reference/classify_cases.md)),
are counted unless `exclude_masked = TRUE` leaves them out of the
numerator and both denominators, as a study that sees only the relapses
that bring a new clone would. `NU_HYPNOZOITE_BATCH` was fitted with
every relapse counted (0.63 at its shipped value of 3, inside the
calibration target), and leaving the masked relapses out gives 0.37,
below it (section 13 of the [model
description](https://shaziaruybal.github.io/malariafwd/MODEL_DESCRIPTION.html)).

## See also

[`classify_cases()`](https://shaziaruybal.github.io/malariafwd/reference/classify_cases.md),
[`summarise_run()`](https://shaziaruybal.github.io/malariafwd/reference/summarise_run.md)

## Examples

``` r
relapse_fraction(example_run())
#>    n_relapse n_new_infection n_recrudescence relapse_frac_episodes
#>        <int>           <int>           <int>                 <num>
#> 1:       158             136               5             0.5284281
#>    relapse_frac_acquisitions
#>                        <num>
#> 1:                  0.537415
```
