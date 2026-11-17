# A run's epidemiology over one window of days

One row that summarises a run over a window of its days, by default the
second half. Every column covers the same days, so the rows of different
runs, or of a grid, mean the same thing.

## Usage

``` r
summarise_run(
  run,
  burn_frac = 0.5,
  from = NULL,
  to = NULL,
  exclude_masked = FALSE
)
```

## Arguments

- run:

  a run: a `vivax_run` (from
  [`run_sim()`](https://shaziaruybal.github.io/malariafwd/reference/run_sim.md)
  or
  [`load_run()`](https://shaziaruybal.github.io/malariafwd/reference/load_run.md))
  or the path of a run directory.

- burn_frac:

  the share of `SIM_DAYS` left out as burn-in, at least 0 and below 1:
  the window starts on day `round(SIM_DAYS * burn_frac) + 1`. Not with
  `from`.

- from, to:

  the window's first and last days, both included: `from = 0` is day 1,
  and `to = NULL` is `SIM_DAYS`.

- exclude_masked:

  if `TRUE`, leave masked relapses out of every relapse column.

## Value

a one-row data.table, nothing rounded:

- the window and the run: `from_tick`, `to_tick` and `n_days` (integer),
  `exclude_masked`, `extinct` and `never_started` (logical), and
  `extinct_tick` (integer, `NA` for a run that did not go extinct);

- host state, over the window's host-days: `prevalence_frac` (the share
  with a patent clone), `coi_mean` and `coi_median` (patent lineages per
  infected host-day) and `polyclonal_frac` (the share of infected
  host-days with two or more), each `NA` with no infected host-day, and
  `hypnozoites_per_host`;

- rates per person-year: `eir_per_py` (infectious bites; `NA` for a
  state table without `infectious_bites`), `foi_per_py` (new clones from
  bites and imports), `molfob_clones_per_py` (and from relapses),
  `molfob_episodes_per_py` (and from recrudescences) and
  `relapses_per_py_at_risk` (`NA` with no time at risk);

- fractions: `relapse_frac_episodes` and `relapse_frac_acquisitions`, as
  [`relapse_fraction()`](https://shaziaruybal.github.io/malariafwd/reference/relapse_fraction.md)
  gives them, and `imported_frac` (imports among the new clones from
  bites and imports), each `NA` with nothing to divide;

- counts (integer): `n_imported`, `n_local_reinfection`, `n_relapse`,
  `n_relapse_masked`, `n_recrudescence` and `n_bloodstage_courses`
  (blood-stage courses, by case management or MDA);

- denominators: `person_years` and `person_years_at_risk`.

## Details

**The window.** The first `burn_frac` of the run's `SIM_DAYS` days are
burn-in, so the window is days `round(SIM_DAYS * burn_frac) + 1` to
`SIM_DAYS`; `from` and `to` set it directly. A run that went extinct, or
never started, is summarised over the same days, the days after its last
state row counting as days with no parasite. A run that stopped early
without either record is refused unless `to` is at most its last
recorded day.

**What is counted.** Cases are the clone onsets in the window, as
[`classify_cases()`](https://shaziaruybal.github.io/malariafwd/reference/classify_cases.md)
gives them: a bite that delivers two genotypes is two new clones, though
one infectious bite, and founder seeding is not an onset. A
recrudescence is an episode but not an acquisition, so it is in
`molfob_episodes_per_py` and `relapse_frac_episodes` only. Masked
relapses, into a patent clone of their own lineage, are counted unless
`exclude_masked = TRUE` leaves them out of every relapse column;
`n_relapse_masked` counts them either way. `NU_HYPNOZOITE_BATCH` was
fitted with every relapse counted (see
[`relapse_fraction()`](https://shaziaruybal.github.io/malariafwd/reference/relapse_fraction.md)).

**Person-years.** The EIR, FOI and molFOB divide by the person-years of
every host, `N_HOSTS` times the window's days over 365. The relapse rate
divides by the person-years at risk: each host counts from the day after
its first infection by a bite or an import, a blood-stage primary or a
hypnozoite batch, to the window's end; founder seeding does not put a
host at risk.

**Pooling.** To pool runs, sum the counts and the person-years and
divide, or weight each rate by its person-years: an average of rates
weights a short or small run as much as a long one.

## See also

[`force_of_infection()`](https://shaziaruybal.github.io/malariafwd/reference/force_of_infection.md),
[`relapse_fraction()`](https://shaziaruybal.github.io/malariafwd/reference/relapse_fraction.md),
[`classify_cases()`](https://shaziaruybal.github.io/malariafwd/reference/classify_cases.md),
[`plot_diagnostics()`](https://shaziaruybal.github.io/malariafwd/reference/plot_diagnostics.md)

## Examples

``` r
run <- example_run()
summarise_run(run)                        # the second half of the run
#>    from_tick to_tick n_days exclude_masked extinct extinct_tick never_started
#>        <int>   <int>  <int>         <lgcl>  <lgcl>        <int>        <lgcl>
#> 1:       201     400    200          FALSE   FALSE           NA         FALSE
#>    prevalence_frac coi_mean coi_median polyclonal_frac hypnozoites_per_host
#>              <num>    <num>      <num>           <num>                <num>
#> 1:       0.2503333 1.338549          1       0.2563249             1.705667
#>    eir_per_py foi_per_py molfob_clones_per_py molfob_episodes_per_py
#>         <num>      <num>                <num>                  <num>
#> 1:   3.680417       2.92             6.722083               6.782917
#>    relapses_per_py_at_risk relapse_frac_episodes relapse_frac_acquisitions
#>                      <num>                 <num>                     <num>
#> 1:                 5.55657             0.5605381                 0.5656109
#>    imported_frac n_imported n_local_reinfection n_relapse n_relapse_masked
#>            <num>      <int>               <int>     <int>            <int>
#> 1:     0.1354167         13                  83       125               69
#>    n_recrudescence n_bloodstage_courses person_years person_years_at_risk
#>              <int>                <int>        <num>                <num>
#> 1:               2                   58     32.87671             22.49589
summarise_run(run, from = 101, to = 200)  # its second hundred days
#>    from_tick to_tick n_days exclude_masked extinct extinct_tick never_started
#>        <int>   <int>  <int>         <lgcl>  <lgcl>        <int>        <lgcl>
#> 1:       101     200    100          FALSE   FALSE           NA         FALSE
#>    prevalence_frac coi_mean coi_median polyclonal_frac hypnozoites_per_host
#>              <num>    <num>      <num>           <num>                <num>
#> 1:       0.1376667 1.179177          1       0.1791768            0.6933333
#>    eir_per_py foi_per_py molfob_clones_per_py molfob_episodes_per_py
#>         <num>      <num>                <num>                  <num>
#> 1:   1.885833   1.338333                 2.92                 3.1025
#>    relapses_per_py_at_risk relapse_frac_episodes relapse_frac_acquisitions
#>                      <num>                 <num>                     <num>
#> 1:                4.186149             0.5098039                 0.5416667
#>    imported_frac n_imported n_local_reinfection n_relapse n_relapse_masked
#>            <num>      <int>               <int>     <int>            <int>
#> 1:     0.1363636          3                  19        26               15
#>    n_recrudescence n_bloodstage_courses person_years person_years_at_risk
#>              <int>                <int>        <num>                <num>
#> 1:               3                   17     16.43836             6.210959
```
