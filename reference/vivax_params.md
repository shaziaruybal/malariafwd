# Build a validated parameter set for the vivax model

Returns a named list of all model parameters, starting from the defaults
(which mirror the SLiM model) and overriding any you name. `SEED` and
`OUTDIR` are set by
[`run_sim()`](https://shaziaruybal.github.io/malariafwd/reference/run_sim.md),
not here.

## Usage

``` r
vivax_params(...)

# S3 method for class 'vivax_params'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- ...:

  named parameter overrides, e.g. `N_HOSTS = 300, BITING_RATE = 0.09`;
  ignored by
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html).

- x:

  a `vivax_params` object.

- row.names, optional:

  ignored.

## Value

a named list of every model parameter (class `vivax_params`): the
defaults with your overrides, checked.
[`print()`](https://rdrr.io/r/base/print.html) shows the values that
differ from the defaults, and
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) gives the
table: `name`, `value` (as SLiM is given it), `default`, `changed`,
`group`, `unit` and `meaning`.

## Details

The defaults are an illustrative preview, 40 hosts for a year, not a
calibrated setting: at them the relapse rate is above its calibration
target even with `NU_HYPNOZOITE_BATCH = 3`. `NU_HYPNOZOITE_BATCH = 3` is
the value fitted to four calibration targets, ranges taken from
published *P. vivax* studies, at moderate transmission (`BITING_RATE`
0.08 with 2 mosquitoes per host and 400 hosts): section 13 of the [model
description](https://shaziaruybal.github.io/malariafwd/MODEL_DESCRIPTION.html)
gives the targets and their sources, and
`inst/validation/calibrate_nu.R` repeats the fit. Re-fit it when
transmission or scale change, since both move the relapse rate.
`BITING_RATE`, with `N_MOSQ / N_HOSTS`, sets the transmission intensity
and is for each study to choose.

Time is in days: one model tick is one day. A duration parameter `D`
(such as `MEAN_BLOOD_STAGE_DURATION` or `MOSQ_LIFESPAN`) is the mean of
a geometric duration: the event happens each day with probability `1/D`.
The treatment delay (`TREATMENT_DELAY`) and, under the default
`RECRUD_DELAY_DIST = "geometric"`, the recrudescence delay
(`MEAN_RECRUD_DELAY`) follow the same rule from the day after the case
or the suppression, so they are at least one day; `EIP_DAYS` and
`GAMETOCYTE_DELAY` are fixed delays. Rates (`BITING_RATE`,
`ALPHA_ACTIVATE`, `MU_HYPNOZOITE`, `DELTA_LATENCY`) are daily
probabilities; a mosquito bites at most once a day, so the realised
biting rate is capped at 1, and `vivax_params()` warns when the peak
over the run, with its seasonal, intervention and resurgence
multipliers, passes it.

Values the model cannot run are refused here, before SLiM starts, naming
the parameter; the *Parameters* section below gives what each may be,
and
[`param_info()`](https://shaziaruybal.github.io/malariafwd/reference/param_info.md)
the same as a table. Each value is a number, `TRUE`/`FALSE` or a string,
as its default is; numbers are finite, and strings are UTF-8 text
without control characters. The tables are also in section 10 of the
[model
description](https://shaziaruybal.github.io/malariafwd/MODEL_DESCRIPTION.html).

Several sizes are capped far above anything realistic, so a cap binds
only on an absurd value; each parameter's row below gives its own,
[`uptake_moments()`](https://shaziaruybal.github.io/malariafwd/reference/uptake_moments.md)
the oocyst one. A founder-free run whose imports never establish a
parasite, by chance or because its settings let none take, is not
refused here; it ends as never started (see
[`load_run()`](https://shaziaruybal.github.io/malariafwd/reference/load_run.md)).
Imports arrive on days 1 to `SIM_DAYS`: a pulse or an import window set
outside the run brings nothing, without a warning, which is how imports
are switched off in one arm of a grid. Each host takes at most one
import a day, so a day's imports beyond `N_HOSTS` never arrive. This
warns when a day's fixed imports inside the run exceed `N_HOSTS`
(`N_IMPORTS` on `IMPORT_TICK`, `N_IMPORTS_PER_TICK` in the window, or
both on the same day), and when `IMPORT_RATE_PER_TICK` at its peak over
the window, with its seasonal and resurgence multipliers and with
`N_IMPORTS_PER_TICK`, draws more than `N_HOSTS` imports on a day with
probability above 0.01.

Each value is stored as its default's type and as SLiM is given it (a
double to the 15 significant digits
[`run_sim()`](https://shaziaruybal.github.io/malariafwd/reference/run_sim.md)
passes), so a set rebuilt from a run's `params.tsv` is the set it used.
Common forms are converted, such as a named number or a whole number
carrying floating-point noise; a number of another class, such as a
`bit64::integer64`, is refused: give the number itself.

**Editing a set.** A set edited with `$<-`, `[[<-` or
[`modifyList()`](https://rdrr.io/r/utils/modifyList.html) keeps its
class and is not checked at the edit. It is checked in full when it is
used:
[`run_sim()`](https://shaziaruybal.github.io/malariafwd/reference/run_sim.md)
and
[`run_grid()`](https://shaziaruybal.github.io/malariafwd/reference/run_grid.md)
build every set again through `vivax_params()`, and
`p <- do.call(vivax_params, p)` checks it at once. A default taken from
another parameter is set only when a set is built, so editing
`N_FOUNDER_HOSTS` leaves `N_FOUNDER_LINEAGES` as it was.

**Choosing a value.** The parameters below give their bounds and, where
a value is easy to get wrong, a `Setting:` note on how to choose one.
[`vignette("parameters")`](https://shaziaruybal.github.io/malariafwd/articles/parameters.md)
collects those notes and the cross-parameter advice they depend on, such
as the scale dependence of transmission.

## References

The works the sources column cites (full list: the model description's
references):

Ayele T, et al. (2023) Curr Res Parasitol Vector Borne Dis 4:100137.
[doi:10.1016/j.crpvbd.2023.100137](https://doi.org/10.1016/j.crpvbd.2023.100137)

Jiang H, et al. (2011) Genome Biol 12:R33.
[doi:10.1186/gb-2011-12-4-r33](https://doi.org/10.1186/gb-2011-12-4-r33)

Koepfli C, et al. (2021) PLoS Negl Trop Dis 15:e0009672.
[doi:10.1371/journal.pntd.0009672](https://doi.org/10.1371/journal.pntd.0009672)

Mehra S, Taylor PG, McCaw JM, Flegg JA (2024a) J Math Biol 89:7.
[doi:10.1007/s00285-024-02088-7](https://doi.org/10.1007/s00285-024-02088-7)

Mehra S, Nosten F, Luxemburger C, White NJ, Watson JA (2024b) Proc Natl
Acad Sci USA 121:e2401024121.
[doi:10.1073/pnas.2401024121](https://doi.org/10.1073/pnas.2401024121)

Ross A, et al. (2025) PLoS Negl Trop Dis 19:e0013567.
[doi:10.1371/journal.pntd.0013567](https://doi.org/10.1371/journal.pntd.0013567)

Taylor AR, et al. (2019) Nat Commun 10:5595.
[doi:10.1038/s41467-019-13412-x](https://doi.org/10.1038/s41467-019-13412-x)

White MT, et al. (2014) eLife 3:e04692.
[doi:10.7554/eLife.04692](https://doi.org/10.7554/eLife.04692)

White MT, et al. (2018) Nat Commun 9:3300.
[doi:10.1038/s41467-018-05860-8](https://doi.org/10.1038/s41467-018-05860-8)

## See also

[`param_info()`](https://shaziaruybal.github.io/malariafwd/reference/param_info.md),
[`run_sim()`](https://shaziaruybal.github.io/malariafwd/reference/run_sim.md),
[`load_params()`](https://shaziaruybal.github.io/malariafwd/reference/load_params.md),
[`uptake_calibrate()`](https://shaziaruybal.github.io/malariafwd/reference/uptake_calibrate.md)

## Parameters

Every parameter, by group. Section 10 of the [model
description](https://shaziaruybal.github.io/malariafwd/MODEL_DESCRIPTION.html)
holds the same tables.

### Run

Length, label and what is recorded.

- `SIM_DAYS`:

  Length of the run; one tick is one day. A run that goes extinct, or
  never starts, ends sooner. Unit: days. Default 365; whole number in
  \[2, 1e9\]; a run that reaches day 99999 must move the days that
  default to 99999 past its end.

- `SIM_ID`:

  A label for the run, written into params.tsv. Keep it to ASCII on
  Windows: a non-ASCII label does not survive the way the run's
  arguments reach SLiM there, and the run fails. Default "default";
  text; UTF-8, without control characters.

- `VERBOSE`:

  Narrate every event in run.log. The outputs are the same either way.
  Default FALSE; TRUE or FALSE.

- `RETAIN_FULL_HISTORY`:

  Keep every genome in the tree sequence, which IBD to infections that
  have already cleared needs. With FALSE most earlier genomes leave the
  tree, and ibd_pairs() warns that the ids it was asked about are
  absent; they look up as NA. Default TRUE; TRUE or FALSE.

### Population

Hosts, mosquitoes and the founding infections.

- `N_HOSTS`:

  Humans; the number is constant. Unit: hosts. Default 40; whole number
  \>= 1. Setting: prevalence and COI are scale-dependent, so a
  transmission rate fitted at one population size does not transfer:
  calibrate at the size you will publish (roadmap S10).

- `N_MOSQ`:

  Mosquitoes; the number is constant, and every one starts uninfected.
  Unit: mosquitoes. Default 100; whole number \>= 1. Setting: the
  transmission rate depends only on BITING_RATE x N_MOSQ / N_HOSTS, so
  the two trade off, but the genetics do not: spreading the same uptakes
  over more mosquitoes lowers co-transmission. Hold it fixed and vary
  BITING_RATE where a result turns on co-transmission, and calibrate at
  publication scale (roadmap S10).

- `N_FOUNDER_HOSTS`:

  Hosts infected at the start: establishment, whether the epidemic
  sustains itself. Unit: hosts. Default 12; whole number \>= 0; at most
  N_HOSTS; 0 needs an import path inside the run. Setting: hold this
  fixed and vary N_FOUNDER_LINEAGES to change founding diversity without
  changing whether the epidemic establishes.

- `N_FOUNDER_LINEAGES`:

  Distinct founding genomes among the founder hosts; the others carry
  clonal copies, seeded on day 2. The genealogical diversity knob: fewer
  lineages give higher background relatedness at the same establishment.
  Unit: genomes. Default N_FOUNDER_HOSTS; whole number \>= 0; at most
  N_FOUNDER_HOSTS, and at least 1 when there are founder hosts. Setting:
  the clean genealogical diversity knob: lower it, at fixed
  N_FOUNDER_HOSTS, to raise background relatedness without changing
  establishment.

### Transmission

Bites, infection and the mosquito's incubation.

- `BITING_RATE`:

  Daily probability that a mosquito bites. Bites per human per day are
  BITING_RATE x N_MOSQ / N_HOSTS, not this value. Unit: per mosquito per
  day. Default 0.1; number \> 0; a mosquito bites at most once a day, so
  it warns when the peak over the run, with its seasonal, intervention
  and resurgence multipliers, passes 1. Setting: the main transmission
  knob, but prevalence and complexity of infection are scale-dependent,
  so a value fitted on a small population does not transfer: calibrate
  it at the population size you will publish (roadmap S10).

- `P_TRANS_VH`:

  b: an infectious bite infects the human. Unit: probability. Default
  0.5; number in \[0, 1\].

- `P_TRANS_HV`:

  c: a mosquito that bites an infectious human takes up parasites. Unit:
  probability. Default 0.5; number in \[0, 1\].

- `MOSQ_LIFESPAN`:

  Mean mosquito lifespan: a mosquito dies each day with probability
  1/MOSQ_LIFESPAN. Unit: days. Default 10; number \>= 1.

- `EIP_DAYS`:

  Extrinsic incubation period: days from uptake until a mosquito's
  products can be transmitted (fixed). Unit: days. Default 10; whole
  number \>= 0.

- `MAX_MOI`:

  A host already carrying this many blood-stage lineages, patent or
  subpatent, takes no new inoculation, counting the clones it gained
  earlier that day. Relapses bypass it, and a bite that delivers several
  genotypes is admitted whole. 1 means no superinfection; at 20 it is a
  guard that seldom binds. Unit: lineages. Default 20; whole number
  \>= 1. Setting: a runaway guard, not a quantity to vary: at 20 it
  essentially never binds, and set low it becomes what you are
  measuring. Confirm it never bound in a published run (roadmap S8).

- `MAX_LINEAGES_PER_MOSQ`:

  A mosquito carrying this many lineages takes up no more. One uptake
  can add up to four per oocyst, so the count can pass it. Unit:
  lineages. Default 20; whole number \>= 0.

### Blood stage

A clone's duration and its infectiousness to mosquitoes.

- `MEAN_BLOOD_STAGE_DURATION`:

  Mean duration of a blood-stage clone: it clears each day with
  probability 1/MEAN_BLOOD_STAGE_DURATION. Unit: days. Default 60;
  number \>= 1. Source: White et al. (2014). Setting: sets how much
  episodes overlap, and so polyclonality and the outcrossing that
  follows from it: sweep it when a result depends on COI or on
  relatedness within a host (roadmap S9).

- `P_PRIMARY`:

  A transmitted genotype causes a blood-stage primary; its hypnozoite
  batch is drawn either way. Unit: probability. Default 1; number in
  \[0, 1\]. Source: White et al. (2014), base model.

- `P_GAMETOCYTE`:

  An infection carries gametocytes, drawn once per infection. Vivax:
  about 0.9. Unit: probability. Default 1; number in \[0, 1\]. Setting:
  ships at 1, the identity: every infection is transmissible. Sweep a
  high but incomplete carriage of about 0.9 together with
  GAMETOCYTE_DELAY (roadmap S5).

- `GAMETOCYTE_DELAY`:

  Days from blood-stage onset until a clone can infect mosquitoes
  (fixed). 0 acts as 1: a clone is first taken up the day after it
  appears. Vivax: about 2 to 3 days. Unit: days. Default 0; whole number
  \>= 0. Setting: ships at 0, the identity: every blood-stage clone is
  transmissible at once. Vivax gametocytes appear early but not
  instantly, so sweep a short delay of about 3 days (roadmap S5).

### Mosquito uptake and co-transmission

How many oocysts an uptake forms, and how many genotypes a bite
delivers. Both lambdas are pre-truncation rates; 0 for both gives one
oocyst per uptake and one genotype per bite, and the defaults let a bite
co-transmit several.

- `LAMBDA_OOCYSTS`:

  Pre-truncation rate of the zero-truncated number of oocysts per
  successful uptake, each its own meiosis; 0 gives exactly one, and the
  default a mean of 2.31. See uptake_moments() and uptake_calibrate().
  Unit: oocysts. Default 2; number in \[0, 1000\]; with K_OOCYSTS, the
  mean count per uptake is at most 1000. Source: membrane feeds of wild
  An. arabiensis in Ethiopia give a median of 7.2 oocysts per midgut
  (Ayele et al. 2023). Setting: a rate before zero-truncation, not the
  mean oocysts a midgut carries: convert with uptake_calibrate() first.
  Vivax feeding counts are overdispersed and far higher than the
  default, so sweep it with K_OOCYSTS (roadmap S3).

- `K_OOCYSTS`:

  Negative-binomial dispersion of the oocyst count (variance lambda +
  lambda^2 / K); 0 is Poisson. It applies to oocysts only. Vivax oocyst
  counts are overdispersed. Default 0; number \>= 0; 0, or at least
  0.001. Setting: leave at 0 for a Poisson count; raise it to match the
  overdispersion of vivax feeding studies, which the default does not
  have (roadmap S3).

- `LAMBDA_SPOROZOITES`:

  Pre-truncation rate of the zero-truncated number of a mosquito's
  products an infectious bite samples, repeats of a lineage counted
  once, so a bite can deliver several genotypes; 0 gives one genotype
  per bite. No study has measured the genotypes a natural bite delivers,
  so the default is an assumption to sweep. Unit: products per bite.
  Default 2; number in \[0, 1000\]. Setting: a rate before
  zero-truncation, not the mean genotypes a feeding study reports:
  convert with uptake_calibrate() first. No measurement of genotypes per
  natural bite exists, so sweep it rather than trust a point value
  (roadmap S1).

### Hypnozoite reservoir

The liver-stage batches a bite seeds, and their activation and death.

- `NU_HYPNOZOITE_BATCH`:

  Mean of the geometric hypnozoite batch an infectious bite seeds, split
  over the genotypes it delivers. Fitted to four calibration targets at
  moderate transmission (section 13 of the model description); re-fit it
  per setting. Unit: hypnozoites per bite. Default 3; number in \[0,
  1e6\]; times its boost, at most 1e6 once the trigger falls inside the
  run. Source: White et al. (2014), section 2. Setting: calibrated to 3
  against relapse, molFOB, MOI and relapse-fraction targets for one
  setting, and it drives relapse rate, reservoir size and COI together:
  re-fit it per setting rather than carrying 3 over (roadmap S2).

- `ALPHA_ACTIVATE`:

  Daily probability that a hypnozoite activates: a relapse. Unit: per
  hypnozoite per day. Default 1/170; number in \[0, 1\]. Source: an
  average carriage of about 170 days (Mehra et al. 2024b); White et
  al. (2014) estimate 1/332. Setting: with K_LATENCY, this sets relapse
  periodicity and so how identifiable a relapse is: report tropical and
  temperate settings separately rather than averaging them (roadmap S7).

- `MU_HYPNOZOITE`:

  Daily probability that a hypnozoite dies, also while latent. Unit: per
  hypnozoite per day. Default 1/442; number in (0, 1\]. Source: as in
  Mehra et al. (2024a), who credit White et al. (2014), whose table
  gives 1/425.

- `K_LATENCY`:

  Latency stages a batch passes before it can activate: 0 is tropical, 2
  temperate. Unit: stages. Default 0; whole number \>= 0. Source: White
  et al. (2014). Setting: left at 0 the reservoir is tropical, with no
  delay before a hypnozoite can activate; 2 with DELTA_LATENCY gives the
  temperate pattern. Report the two separately (roadmap S7).

- `DELTA_LATENCY`:

  Daily probability that a latent batch moves on one stage. Unit: per
  day. Default 1/100; number in (0, 1\]; checked only when K_LATENCY \>
  0.

### Genome

The neutral PvP01 genome the model simulates; the mutation rate is 0.

- `GENOME_MODE`:

  "full": the 14 PvP01 chromosomes at full length (24.2 Mb); "core":
  each at 90% of its length, a placeholder for the core genome. Default
  "full"; "full" or "core".

- `CHROMOSOMES`:

  Which PvP01 chromosomes to simulate, such as "1-14" or "1-3,5".
  Default "1-14"; text; chromosome numbers 1 to 14, each once, separated
  by commas, with low-high ranges.

- `RECOMB_RATE`:

  Crossover rate within a chromosome; chromosomes assort independently.
  Run time grows steeply with it. Unit: per bp. Default -1; number \<=
  0.5; 0 or below means the default map rate, 7.4e-7. Source: Jiang et
  al. (2011).

### Seasonality and vector control

Factors on the biting rate; seasonality also acts on the Poisson import
rate.

- `SEASONAL_AMP`:

  Amplitude A of the seasonal factor 1 + A cos(2 pi (t - SEASONAL_PHASE)
  / 365), floored at 0, on biting and on the Poisson import rate; 0
  means no seasonality. Default 0; number \>= 0.

- `SEASONAL_PHASE`:

  Day of the year of the seasonal peak. Unit: day of the year. Default
  0; number.

- `INTERVENTION_START`:

  First day of vector control. Unit: day. Default 99999 (never); whole
  number.

- `INTERVENTION_END`:

  Vector control acts on days INTERVENTION_START to
  INTERVENTION_END - 1. Unit: day. Default 99999 (never); whole number.

- `INTERVENTION_EFFICACY`:

  Biting is multiplied by 1 - INTERVENTION_EFFICACY inside the window;
  importation and treatment are not affected. Unit: fraction. Default 0;
  number in \[0, 1\].

### Detection and case management

A case is detected, decided where it arises and treated later, on the
host; radical cure is nested in a blood-stage course. Off by default.

- `P_DETECT_PRIMARY`:

  A primary infection or an import is detected, drawn once per bite (for
  all the genotypes it delivers) or per import. Unit: probability.
  Default 1; number in \[0, 1\].

- `P_DETECT_RELAPSE`:

  A relapse is detected, drawn once per activation. Relapses are often
  subpatent, so a realistic value is lower. Unit: probability. Default
  1; number in \[0, 1\]. Setting: relapses are often subpatent, and
  every observational result depends on this: sweep it over 0.15 to 1.0
  rather than leaving it at the default of 1, which reports every
  relapse (roadmap S6).

- `BLOODCURE_COVERAGE`:

  A detected case is treated: a course clears each patent clone the host
  carries, unless the clone survives. 0 turns the whole treatment path
  off. Unit: probability. Default 0; number in \[0, 1\].

- `P_BLOODCURE_SURVIVAL`:

  A clone survives a course: a patent one, or under MDA a subpatent one
  too. Unit: probability. Default 0; number in \[0, 1\]. Source: about
  0.10 in the literature.

- `RADCURE_COVERAGE`:

  A treated case also gets radical cure: a probability given treatment,
  not a coverage of its own, so it does nothing unless
  BLOODCURE_COVERAGE \> 0. Unit: probability. Default 0; number in \[0,
  1\].

- `P_RAD_SURVIVAL`:

  Radical cure fails: the patient's whole reservoir survives
  (RADCURE_MODE "person"), about 71% cure; or each hypnozoite survives
  ("hypnozoite"), where the same number would mean 95% efficacy. Unit:
  probability. Default 0.29; number in \[0, 1\]. Source: near the 70%
  effectiveness White et al. (2018) assume for a 14-day primaquine
  course, given poor adherence.

- `RADCURE_MODE`:

  The unit of radical-cure failure: "person", or "hypnozoite", which
  thins the reservoir hypnozoite by hypnozoite and is the only mode that
  can clear part of it. Default "person"; "person" or "hypnozoite".

- `TREATMENT_DELAY`:

  Mean days from a case to its course: geometric, at least one day, so a
  course never lands on the day the case arises. Unit: days. Default 7;
  number in \[1, 1e6\]. Source: White et al. (2014) use 7 days.

- `TREATMENT_START`:

  First day on which a detected case can be prescribed a course; case
  management has its own window, separate from vector control. Unit:
  day. Default 0; whole number.

- `TREATMENT_END`:

  Cases detected on days TREATMENT_START to TREATMENT_END - 1 are
  treated. The window gates the decision, not the delivery: a course
  prescribed before the end still arrives. Unit: day. Default 99999
  (never); whole number.

### Mass drug administration

One round, through the same course as case management.

- `MDA_TICK`:

  Day of one MDA round, through the same course as case management:
  blood-stage drugs and radical cure, whatever RADCURE_COVERAGE is, for
  each host it reaches. It also reaches subpatent clones and clones that
  became patent that day, and it is not gated by the treatment window.
  Unit: day. Default 99999 (never); whole number; not 1 when a founder
  host carries a copy of another founder's genome.

- `MDA_COVERAGE`:

  The round reaches a host, drawn for each host. Unit: probability.
  Default 1; number in \[0, 1\].

### Recrudescence

A clone that survives a course can be suppressed below patency and
return. Off by default.

- `P_RECRUD_GIVEN_FAIL`:

  A clone that survives a course is suppressed below patency, to
  recrudesce later, rather than persisting. At 0 the third recurrence
  state is inert, which matches recrudescence being 0.2 to 0.3% of
  recurrences after chloroquine or artesunate. Unit: probability.
  Default 0; number in \[0, 1\]. Source: Taylor et al. (2019).

- `MEAN_RECRUD_DELAY`:

  Mean days from suppression to recrudescence, at least one; a
  literature assumption pending parameterisation. Geometric at the
  default: a median of 7 days, about 94% within 28 days. Unit: days.
  Default 10.5; number in \[1, 1e6\]. Source: close to the prior of
  Taylor et al. (2019). Setting: a literature assumption rather than an
  estimate from data, and it sets how many recrudescences a 28-day
  follow-up sees: sweep it, though it matters only when
  P_RECRUD_GIVEN_FAIL is above 0 (roadmap S11).

- `RECRUD_DELAY_DIST`:

  Shape of that delay: "geometric", the discrete form of an exponential
  delay; "poisson", 1 + Poisson(D - 1) with mean D, a tight-window
  sensitivity arm; or "fixed", round(D) days with halves rounded up (11
  days at 10.5). Every shape is at least one day. Default "geometric";
  "geometric" or "poisson" or "fixed".

- `P_DETECT_RECRUD`:

  A recrudescence is detected; it should usually sit well below
  P_DETECT_PRIMARY. Unit: probability. Default 0.5; number in \[0, 1\].
  Source: Ross et al. (2025) put recurrence detectability at 0.13 to
  0.49 of a primary's.

- `SUBPATENT_INFECTIVITY`:

  Infectiousness of a suppressed clone relative to a patent one. Neither
  0 nor 1 is defensible. Unit: relative. Default 0.15; number in \[0,
  1\]. Source: 0.035 to 0.10 in the model of White et al. (2018); 0.10
  to 0.20 from the feeding studies compiled by Koepfli et al. (2021, Fig
  5).

### Resurgence

Multipliers that switch on from one day, to start a resurgence.

- `RESURGENCE_TRIGGER_TICK`:

  Day from which the four boosts apply. Unit: day. Default 99999
  (never); whole number.

- `BITING_RATE_BOOST`:

  Multiplies BITING_RATE from the trigger. Unit: multiplier. Default 1;
  number \>= 0.

- `ALPHA_ACTIVATE_BOOST`:

  Multiplies ALPHA_ACTIVATE from the trigger, capped so that activation
  and death stay within one day's probability. Unit: multiplier. Default
  1; number \>= 0.

- `NU_HYPNOZOITE_BATCH_BOOST`:

  Multiplies NU_HYPNOZOITE_BATCH from the trigger. Unit: multiplier.
  Default 1; number \>= 0.

### Importation

Three routes, which add up: a pulse, and a fixed daily count and a
Poisson rate in a window. Each host takes at most one import a day, and
a host at MAX_MOI takes none.

- `IMPORT_TICK`:

  Day of one pulse of N_IMPORTS imports. Unit: day. Default 99999
  (never); whole number.

- `N_IMPORTS`:

  Size of the pulse. Unit: imports. Default 0; whole number \>= 0; warns
  above N_HOSTS.

- `IMPORT_START`:

  First day of the import window. Unit: day. Default 99999 (never);
  whole number.

- `IMPORT_END`:

  Imports in the window arrive on days max(1, IMPORT_START) to
  min(IMPORT_END - 1, SIM_DAYS). Unit: day. Default 99999 (never); whole
  number.

- `N_IMPORTS_PER_TICK`:

  A fixed number of imports on each day of the window. Unit: imports a
  day. Default 0; whole number \>= 0; warns when a day's fixed imports
  exceed N_HOSTS.

- `IMPORT_RATE_PER_TICK`:

  Mean of a Poisson number of imports on each day of the window; it
  follows the season. Unit: imports a day. Default 0; number in \[0,
  1e6\]; at its seasonal and boosted peak, at most 1e6 a day.

- `IMPORT_RATE_BOOST`:

  Multiplies IMPORT_RATE_PER_TICK from RESURGENCE_TRIGGER_TICK. Unit:
  multiplier. Default 1; number \>= 0.

## Examples

``` r
vivax_params(N_HOSTS = 300, BITING_RATE = 0.09, GENOME_MODE = "core")
#> <vivax_params> 63 parameters
#>   non-default:
#>     N_HOSTS                  = 300
#>     BITING_RATE              = 0.09
#>     GENOME_MODE              = core
p <- vivax_params(N_HOSTS = 300)
subset(as.data.frame(p), changed)
#>      name value default changed      group  unit
#> 5 N_HOSTS   300      40    TRUE Population hosts
#>                           meaning
#> 5 Humans; the number is constant.
```
