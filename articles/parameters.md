# Choosing parameters

A run takes 63 parameters, and
[`vivax_params()`](https://shaziaruybal.github.io/malariafwd/reference/vivax_params.md)
gives every one a default, so a first run needs none of them.
[`?vivax_params`](https://shaziaruybal.github.io/malariafwd/reference/vivax_params.md)
lists them all with their units and bounds, and
[`param_info()`](https://shaziaruybal.github.io/malariafwd/reference/param_info.md)
returns the same thing as a table. This page is about the smaller set
you have to think about: the ones where a plausible-looking value
quietly changes what you are measuring.

Five things catch people out, and two more are worth knowing before you
start. All seven come before the per-parameter table.

## The defaults are a transparent model, not a realistic one

Five parameters ship at the value that switches their mechanism **off**,
so a first run shows the model’s mechanics with no carriage or
observation layer on top of them:

``` r

off <- param_info(c("P_PRIMARY", "P_GAMETOCYTE", "GAMETOCYTE_DELAY",
                    "P_DETECT_PRIMARY", "P_DETECT_RELAPSE"))
knitr::kable(off[, c("name", "default", "meaning")],
             col.names = c("Parameter", "Default", "What it does"))
```

| Parameter | Default | What it does |
|:---|:---|:---|
| P_PRIMARY | 1 | A transmitted genotype causes a blood-stage primary; its hypnozoite batch is drawn either way. |
| P_GAMETOCYTE | 1 | An infection carries gametocytes, drawn once per infection. Vivax: about 0.9. |
| GAMETOCYTE_DELAY | 0 | Days from blood-stage onset until a clone can infect mosquitoes (fixed). 0 acts as 1: a clone is first taken up the day after it appears. Vivax: about 2 to 3 days. |
| P_DETECT_PRIMARY | 1 | A primary infection or an import is detected, drawn once per bite (for all the genotypes it delivers) or per import. |
| P_DETECT_RELAPSE | 1 | A relapse is detected, drawn once per activation. Relapses are often subpatent, so a realistic value is lower. |

At those defaults every transmitted genotype becomes a blood-stage
infection, is infectious to mosquitoes from the next day, and is
detected. That is deliberate: it means a discrepancy in your output is
the model’s, not that of a layer you did not know was on. It is also not
*P. vivax*. Gametocytes appear early but not instantly, carriage is high
but incomplete, and relapses are often sub-patent — so a result you
intend to publish should move these off their defaults and report the
sensitivity. The per-parameter table below gives the ranges to sweep.

## `BITING_RATE` is per mosquito, not per human

`BITING_RATE` is the daily probability that *one mosquito* bites. It is
not the human biting rate — bites received per person per day,
Ross-Macdonald’s `ma` — which is what “biting rate” usually means in the
malaria literature. The two differ by `N_MOSQ` / `N_HOSTS`:

``` r

p <- vivax_params()
unname(p$BITING_RATE * p$N_MOSQ / p$N_HOSTS)
#> [1] 0.25
```

so the same `BITING_RATE` means different exposure at a different
`N_MOSQ` / `N_HOSTS` ratio.

The *transmission rate* depends only on that product, so `BITING_RATE`
and `N_MOSQ` trade off against each other — but the genetics do not.
Halving `BITING_RATE` and doubling `N_MOSQ` leaves exposure unchanged
while spreading the same uptakes over twice as many mosquitoes, each of
which then carries fewer lineages and so co-transmits and outcrosses
less. For a transmission-rate sweep, vary one of the two; where a result
turns on co-transmission or on recombination in the midgut, hold
`N_MOSQ` fixed and vary `BITING_RATE`.

## A rate fitted at one population size does not transfer

Prevalence and complexity of infection are scale-dependent in their own
right: a small population loses diversity to drift faster, and saturates
differently. So **a transmission rate fitted on 40 hosts does not
transfer to 4,000** — calibrate at the population size you intend to
publish, not at the size that was quick to fit.

## The uptake lambdas are rates before truncation, not mean counts

`LAMBDA_OOCYSTS` and `LAMBDA_SPOROZOITES` parameterise *zero-truncated*
draws: the model has already decided that transmission happened, so the
count it needs is conditioned away from zero. Conditioning pulls the
realised mean up, and the gap is largest exactly where feeding studies
report their numbers. A lambda of 2 does not give a mean of 2:

``` r

uptake_moments(2)          # Poisson
#>     mean      var 
#> 2.313035 1.588974
uptake_moments(2, k = 1)   # negative binomial, same lambda, more spread
#> mean  var 
#>    3    6
```

So a mean oocyst count read off a paper cannot be set as a lambda.
Convert it first. Membrane feeds of wild *An. arabiensis* in Ethiopia
give a median of 7.2 oocysts per midgut (Ayele et al. 2023), well above
the default’s realised 2.31 — to simulate that, ask for it as a target:

``` r

uptake_calibrate(7.2)              # a target realised mean
#> $lambda
#> [1] 7.194595
#> 
#> $k
#> [1] 0
#> 
#> $achieved
#>     mean      var 
#> 7.200000 7.161087
uptake_calibrate(7.2, var = 60)    # and a target variance, if you have one
#> $lambda
#> [1] 5.314459
#> 
#> $k
#> [1] 0.576476
#> 
#> $achieved
#> mean  var 
#>  7.2 60.0
```

[`uptake_calibrate()`](https://shaziaruybal.github.io/malariafwd/reference/uptake_calibrate.md)
returns the `lambda` and `k` to set, with the moments they achieve.

## Two parameters are guards, not quantities

`MAX_MOI` and `MAX_LINEAGES_PER_MOSQ` exist to stop a run growing
without bound. They are not biology, and at their defaults of 20 they
should never bind. Set one low and it stops being a guard and becomes
the thing your run measures — a `MAX_MOI` of 3 does not model an immune
ceiling at 3, it truncates superinfection and silently stops topping up
the hypnozoite reservoir of the hosts it refuses.

Check that they did not bind. A mosquito’s `moi_active` is exactly the
count the cap tests:

``` r

mosq <- load_state(example_run(), entity = "mosquito")
max(mosq$moi_active)                                   # against MAX_LINEAGES_PER_MOSQ
#> [1] 17
```

In this small bundled run just 1 of the 143 mosquitoes ever infected
reached 17 of the 20 lineages it was allowed, on 14 of 4,623 infected
mosquito-days. The guard held, but with only 3 lineages to spare — which
is why it is worth looking rather than assuming.

For a host, `MAX_MOI` counts the blood-stage lineages it carries whether
or not they are patent, so sum the two columns (this omits clones gained
earlier the same day, making it a close lower bound):

``` r

host <- load_state(example_run())
max(host$moi_active + host$moi_subpatent)              # against MAX_MOI
#> [1] 5
```

## Changing diversity without changing the epidemic

Founding diversity and establishment are separable, and conflating them
is the commonest way to confound a relatedness result. `N_FOUNDER_HOSTS`
sets how many hosts start infected, and so whether the epidemic takes
off at all; `N_FOUNDER_LINEAGES` sets how many *distinct* genomes those
hosts carry, the rest being clonal copies. Hold the first fixed and
lower the second to raise background relatedness while leaving
establishment alone:

``` r

diverse <- vivax_params(N_FOUNDER_HOSTS = 12, N_FOUNDER_LINEAGES = 12)  # maximal diversity
related <- vivax_params(N_FOUNDER_HOSTS = 12, N_FOUNDER_LINEAGES = 3)   # same epidemic, more IBD
```

## Tropical or temperate, reported separately

`ALPHA_ACTIVATE` and `K_LATENCY` together set relapse periodicity, and
so how distinguishable a relapse is from a reinfection — the quantity
much of this package exists to measure. At `K_LATENCY = 0` the reservoir
is tropical: a hypnozoite can activate immediately. `K_LATENCY = 2` with
`DELTA_LATENCY` gives the temperate pattern, with a delay first.

These are different epidemiologies, not a range to average over. Report
them as separate scenarios.

## How to choose each of them

The table below is generated from the package’s parameter dictionary, so
it cannot drift from
[`?vivax_params`](https://shaziaruybal.github.io/malariafwd/reference/vivax_params.md).
It covers the parameters that carry explicit guidance; the rest are
documented in
[`?vivax_params`](https://shaziaruybal.github.io/malariafwd/reference/vivax_params.md)
and in section 10 of the [model
description](https://shaziaruybal.github.io/malariafwd/MODEL_DESCRIPTION.html).
The `S` numbers point to the sensitivity-analysis table in the
[roadmap](https://shaziaruybal.github.io/malariafwd/ROADMAP.html), which
gives the evidence behind each recommendation.

### Population

| Parameter | Default | How to choose one |
|:---|:---|:---|
| `N_HOSTS` | `40` | Prevalence and COI are scale-dependent, so a transmission rate fitted at one population size does not transfer: calibrate at the size you will publish (roadmap S10). |
| `N_MOSQ` | `100` | The transmission rate depends only on BITING_RATE x N_MOSQ / N_HOSTS, so the two trade off, but the genetics do not: spreading the same uptakes over more mosquitoes lowers co-transmission. Hold it fixed and vary BITING_RATE where a result turns on co-transmission, and calibrate at publication scale (roadmap S10). |
| `N_FOUNDER_HOSTS` | `12` | Hold this fixed and vary N_FOUNDER_LINEAGES to change founding diversity without changing whether the epidemic establishes. |
| `N_FOUNDER_LINEAGES` | `N_FOUNDER_HOSTS` | The clean genealogical diversity knob: lower it, at fixed N_FOUNDER_HOSTS, to raise background relatedness without changing establishment. |

### Transmission

| Parameter | Default | How to choose one |
|:---|:---|:---|
| `BITING_RATE` | `0.1` | The main transmission knob, but prevalence and complexity of infection are scale-dependent, so a value fitted on a small population does not transfer: calibrate it at the population size you will publish (roadmap S10). |
| `MAX_MOI` | `20` | A runaway guard, not a quantity to vary: at 20 it essentially never binds, and set low it becomes what you are measuring. Confirm it never bound in a published run (roadmap S8). |

### Blood stage

| Parameter | Default | How to choose one |
|:---|:---|:---|
| `MEAN_BLOOD_STAGE_DURATION` | `60` | Sets how much episodes overlap, and so polyclonality and the outcrossing that follows from it: sweep it when a result depends on COI or on relatedness within a host (roadmap S9). |
| `P_GAMETOCYTE` | `1` | Ships at 1, the identity: every infection is transmissible. Sweep a high but incomplete carriage of about 0.9 together with GAMETOCYTE_DELAY (roadmap S5). |
| `GAMETOCYTE_DELAY` | `0` | Ships at 0, the identity: every blood-stage clone is transmissible at once. Vivax gametocytes appear early but not instantly, so sweep a short delay of about 3 days (roadmap S5). |

### Mosquito uptake and co-transmission

| Parameter | Default | How to choose one |
|:---|:---|:---|
| `LAMBDA_OOCYSTS` | `2` | A rate before zero-truncation, not the mean oocysts a midgut carries: convert with uptake_calibrate() first. Vivax feeding counts are overdispersed and far higher than the default, so sweep it with K_OOCYSTS (roadmap S3). |
| `K_OOCYSTS` | `0` | Leave at 0 for a Poisson count; raise it to match the overdispersion of vivax feeding studies, which the default does not have (roadmap S3). |
| `LAMBDA_SPOROZOITES` | `2` | A rate before zero-truncation, not the mean genotypes a feeding study reports: convert with uptake_calibrate() first. No measurement of genotypes per natural bite exists, so sweep it rather than trust a point value (roadmap S1). |

### Hypnozoite reservoir

| Parameter | Default | How to choose one |
|:---|:---|:---|
| `NU_HYPNOZOITE_BATCH` | `3` | Calibrated to 3 against relapse, molFOB, MOI and relapse-fraction targets for one setting, and it drives relapse rate, reservoir size and COI together: re-fit it per setting rather than carrying 3 over (roadmap S2). |
| `ALPHA_ACTIVATE` | `1/170` | With K_LATENCY, this sets relapse periodicity and so how identifiable a relapse is: report tropical and temperate settings separately rather than averaging them (roadmap S7). |
| `K_LATENCY` | `0` | Left at 0 the reservoir is tropical, with no delay before a hypnozoite can activate; 2 with DELTA_LATENCY gives the temperate pattern. Report the two separately (roadmap S7). |

### Detection and case management

| Parameter | Default | How to choose one |
|:---|:---|:---|
| `P_DETECT_RELAPSE` | `1` | Relapses are often subpatent, and every observational result depends on this: sweep it over 0.15 to 1.0 rather than leaving it at the default of 1, which reports every relapse (roadmap S6). |

### Recrudescence

| Parameter | Default | How to choose one |
|:---|:---|:---|
| `MEAN_RECRUD_DELAY` | `10.5` | A literature assumption rather than an estimate from data, and it sets how many recrudescences a 28-day follow-up sees: sweep it, though it matters only when P_RECRUD_GIVEN_FAIL is above 0 (roadmap S11). |

## Sweeping rather than guessing

Where the table says *sweep*,
[`run_grid()`](https://shaziaruybal.github.io/malariafwd/reference/run_grid.md)
takes the parameter sets and the seeds and runs them in parallel:

``` r

grid <- run_grid(
  params = lapply(c(0.15, 0.4, 0.7, 1.0),
                  function(d) vivax_params(P_DETECT_RELAPSE = d, SIM_DAYS = 730)),
  seeds  = 1:20,
  outdir = file.path(tempdir(), "detect-sweep"))
```

Two habits make a sweep interpretable. Vary one mechanism at a time,
since several of these parameters move prevalence, COI and relatedness
together and a joint change cannot be attributed afterwards. And use
enough seeds: at a small `N_HOSTS` the seed-to-seed spread can exceed
the effect you are testing, so a difference between two single runs is
not a result.

## Where the defaults came from

`param_info()$source` gives the origin of each default in one line, and
[`?vivax_params`](https://shaziaruybal.github.io/malariafwd/reference/vivax_params.md)
shows it alongside the parameter. Section 13 of the [model
description](https://shaziaruybal.github.io/malariafwd/MODEL_DESCRIPTION.html)
records what has been calibrated and against which targets, and what
remains a literature assumption — `NU_HYPNOZOITE_BATCH` was fitted for
one setting and should be re-fitted for yours, and `MEAN_RECRUD_DELAY`
has no estimate from data at all. The roadmap’s sensitivity table lists
the eleven places where a choice most needs testing.
