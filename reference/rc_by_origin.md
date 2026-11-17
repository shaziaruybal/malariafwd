# Case reproduction number (Rc) by true origin

The mean number of secondary cases each case caused, by the case's true
origin. A case is a blood-stage clone onset, a row of
[`classify_cases()`](https://shaziaruybal.github.io/malariafwd/reference/classify_cases.md);
its secondary cases are the targets of its links in
[`transmission_links()`](https://shaziaruybal.github.io/malariafwd/reference/transmission_links.md),
infected by a bite's blood stage or by a relapse of the batch the bite
left.

## Usage

``` r
rc_by_origin(
  run,
  unit = c("episode", "acquisition"),
  kind = c("all", "direct", "liver"),
  from = 0,
  followup = 0
)
```

## Arguments

- run:

  a run: a `vivax_run` (from
  [`run_sim()`](https://shaziaruybal.github.io/malariafwd/reference/run_sim.md)
  or
  [`load_run()`](https://shaziaruybal.github.io/malariafwd/reference/load_run.md))
  or the path of a run directory.

- unit:

  `"episode"` or `"acquisition"`: what a case is (see Details).

- kind:

  `"all"`, `"direct"` or `"liver"`: the links counted, as in
  [`transmission_links()`](https://shaziaruybal.github.io/malariafwd/reference/transmission_links.md).

- from:

  count the cases that began on day `from` or later.

- followup:

  count the cases that began on day `SIM_DAYS - followup` or earlier.

## Value

a data.table, one row per origin, every origin present: `origin` (a
factor of the origins
[`classify_cases()`](https://shaziaruybal.github.io/malariafwd/reference/classify_cases.md)
gives, without `"recrudescence"` for `unit = "acquisition"`), `n_cases`
(integer), `n_secondary` (the sum of the cases' link weights), `Rc`
(`n_secondary / n_cases`) and `onward_frac` (the share of cases with a
secondary case), the last two `NA` with no case. To pool runs, divide
the summed `n_secondary` by the summed `n_cases`.

## Details

**Cases.** With `unit = "episode"` every clone onset is a case, a
recrudescence included. With `unit = "acquisition"` a recrudescence is
not a case: its onward transmission is credited to the clone it
re-expanded from, through repeated recrudescences, so each acquisition
has one count. A subpatent holder's transmission is credited to the case
that began it. Founders infect, but are not cases.

**Secondary cases.** Each secondary case counts once: one from an
outcrossed oocyst credits 1/2 to each parent. A bite that delivers two
genotypes gives two secondary cases. Masked relapses count like any
other; join
[`classify_cases()`](https://shaziaruybal.github.io/malariafwd/reference/classify_cases.md)
on the links to leave them out. `kind = "direct"` and `kind = "liver"`
add up to `kind = "all"`.

**The end of the run.** A case near the end has had less time to
transmit, above all through relapses: a hypnozoite stays dormant
`1 / (ALPHA_ACTIVATE + MU_HYPNOZOITE)` days on average, about four
months at the defaults, and more than a year for 1 in 20. `followup`
counts only the cases that began at least that many days before
`SIM_DAYS`. A run that went extinct is not censored.

**The phase of the epidemic.** Rc reflects the phase the counted cases
began in: above 1 while an epidemic grows from its founders, near 1 at a
self-sustaining steady state, and lower where imports sustain
transmission. To compare runs, choose `from` and `followup` inside a
stretch that
[`plot_diagnostics()`](https://shaziaruybal.github.io/malariafwd/reference/plot_diagnostics.md)
shows to be steady; `followup` alone keeps a growing run's early cases.

## References

Fraser C (2007) Estimating individual and household reproduction numbers
in an emerging epidemic. PLoS ONE 2(8): e758.
[doi:10.1371/journal.pone.0000758](https://doi.org/10.1371/journal.pone.0000758)

## See also

[`transmission_links()`](https://shaziaruybal.github.io/malariafwd/reference/transmission_links.md),
[`classify_cases()`](https://shaziaruybal.github.io/malariafwd/reference/classify_cases.md)

## Examples

``` r
rc_by_origin(example_run(), from = 100, followup = 100)
#>               origin n_cases n_secondary       Rc onward_frac
#>               <fctr>   <int>       <num>    <num>       <num>
#> 1:          imported       9         0.0 0.000000   0.0000000
#> 2: local reinfection      35        50.5 1.442857   0.3714286
#> 3:           relapse      67        81.0 1.208955   0.2985075
#> 4:     recrudescence       3         0.0 0.000000   0.0000000
```
