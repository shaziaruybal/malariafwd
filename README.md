# malariafwd <img src="man/figures/logo.png" align="right" height="139" alt="malariafwd hex logo: a mosquito and a human either side of a fast-forward symbol, over a lineage tree" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/shaziaruybal/malariafwd/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/shaziaruybal/malariafwd/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://shaziaruybal.github.io/malariafwd/LICENSE.html)
<!-- badges: end -->

> **Pre-release (0.1.0).** The package is under active development and both the interface and
> the model may change before 1.0.0. 
>
> If you are interested in using *malariafwd* or in contributing, please get in touch by opening
> an [issue](https://github.com/shaziaruybal/malariafwd/issues) or contacting
> Dr Shazia Ruybal-Pesántez (shazia.ruybal@isglobal.org).

*malariafwd* is an individual-based forward simulator of *Plasmodium vivax* transmission and genetics that leverages the
[SLiM](https://messerlab.org/slim/) and [*tskit*](https://tskit.dev/) ecosystem. The model couples host and vector dynamics, co-transmission, meiotic recombination in the mosquito, and a hypnozoite reservoir, with whole-genome simulation, so every parasite genome, including latent hypnozoites, has a known origin and genealogy. Full parasite genealogies across the 14 chromosomes of the PvP01 reference genome are recorded as tree sequences for efficient storage. This also enables access to ancestral recombination graphs (ARG) and provides ground-truth identity-by-descent (IBD). 

*malariafwd* is built for questions whose answers are normally unobservable: for instance, is a recurrence a relapse, a reinfection or a recrudescence? How much of a case load comes from the hypnozoite reservoir? How robustly does an IBD or COI estimator recover the ground truth? How identifiable are IBD and other genetic metrics under given transmission and parasite diversity conditions? 

Further details on the model can be found in the [model documentation](https://shaziaruybal.github.io/malariafwd/MODEL_DESCRIPTION.html).

## Requirements

- **SLiM ≥ 5.2.** On any system, `conda install -c conda-forge slim` (builds for Linux, macOS on
  Intel or Apple silicon, and Windows). On macOS the installer from the
  [SLiM page](https://messerlab.org/slim/) works too; on Linux, conda or a build from source (the
  SLiM manual, chapter 2). On Debian and Ubuntu, `apt install slim` installs an unrelated login
  manager, not this SLiM. Windows is untested. If R cannot find SLiM (RStudio on macOS does not
  see the shell's `PATH`), add `MALARIAFWD_SLIM=/full/path/to/slim` to `~/.Renviron` and restart R.
- For the IBD and ARG functions and `simulate_genotypes()`: **Python
  with `tskit` (≥ 1.0) and `numpy`**, through the `reticulate` package. Either route works.
  *Automatic*: with no Python chosen, reticulate sets up its own environment with these modules the
  first time an IBD function runs; this needs internet once and a few hundred MB in the user cache.
  *Your own Python*: set `RETICULATE_PYTHON=/path/to/python` (in `~/.Renviron` to keep it) before
  the first IBD call of a session.
- On a cluster whose compute nodes have no internet, run `check_setup()` once on a node that has
  it; later sessions that share the same user cache may then reuse the environment. If that does
  not work, use your own Python.

`check_setup()` reports what is found.

## Install

```r
# install.packages("pak")
pak::pak("shaziaruybal/malariafwd", dependencies = TRUE)
```

`dependencies = TRUE` also installs the suggested packages that the plots and the IBD functions
use (ggplot2, patchwork and reticulate); without it you get only what the package needs to load,
and the first plot or IBD call fails. `remotes::install_github("shaziaruybal/malariafwd",
dependencies = TRUE)` does the same thing. The vignettes are on the
[website](https://shaziaruybal.github.io/malariafwd/); to install them as well, use remotes, which
can build them:

```r
remotes::install_github("shaziaruybal/malariafwd", dependencies = TRUE, build_vignettes = TRUE)
```

## Quick start

```r
library(malariafwd)
check_setup()   # SLiM, the model, and Python with tskit
                # the first call builds reticulate's Python environment: a few hundred MB, once

# one run: 100 hosts and 250 mosquitoes for 18 months, a few seconds on a modern laptop.
p   <- vivax_params(N_HOSTS = 100, N_MOSQ = 250, BITING_RATE = 0.1, SIM_DAYS = 548)
run <- run_sim(p, seed = 1, outdir = file.path(tempdir(), "demo"))
# every parameter, its unit and the values it may take: ?vivax_params or param_info()
# no SLiM yet? every call below also takes the run that ships with the package:
# run <- example_run()

# the headline numbers, over the run's second half: prevalence, complexity of infection (COI),
# entomological inoculation rate (EIR), force of infection (FOI), molecular force of blood-stage
# infection (molFOB), relapses
summarise_run(run)

# one call, several panels: did the run settle, saturate or go extinct?
plot_diagnostics(run)

# one host: its blood-stage clones coloured by lineage, over the hypnozoites it carries
ev   <- load_events(run)                 # one row per event; ?events.tsv describes the columns
host <- ev[event_type == "hypnozoite_activation", .N, by = human_id][order(-N), human_id[1]]
plot_host_trajectory(run, host)

# the population: who is infected, and who carries only an invisible reservoir
plot_host_state_prevalence(run)

# clones arriving (molFOB) against infections arriving (FOI)
plot_force_of_infection(run)

# ground-truth IBD by known relationship, read from the recorded genealogy
plot_relatedness_ladder(run, style = "bars")

# the same IBD as segments along the 14 PvP01 chromosomes
tp   <- ibd_truth_pairs(run)
pick <- tp[cls %in% c("relapse vs its origin", "same-oocyst sibs", "distinct founders"),
           .SD[1], by = cls]
plot_ibd_pair(run, pick$a, pick$b, labels = sprintf("%s\nIBD %.2f", pick$cls, pick$ibd))
```

## How it works

The model runs as a **subprocess** (`processx`, shell-free) with parameters passed as `-d`
overrides; it writes `events.tsv`, `state.tsv`, `plasmo.trees`, and a provenance snapshot into
each run directory (`?events.tsv` and `?state.tsv` describe the columns). `reticulate` + `tskit` are used only to compute IBD from the recorded tree
sequence. Multi-run sweeps parallelise via `future`. Set `VERBOSE = TRUE` to see the
full per-event narration; the default is quiet.

## Learn more

- [Getting started](https://shaziaruybal.github.io/malariafwd/articles/malariafwd.html): a run,
  from its parameters to its ground truth and IBD.
- [Choosing parameters](https://shaziaruybal.github.io/malariafwd/articles/parameters.html): the
  parameters that are easy to get wrong, and how to choose a value for each.
- [Reading the outputs](https://shaziaruybal.github.io/malariafwd/articles/anatomy.html): what
  each output holds, and how to read it.
- [Reference](https://shaziaruybal.github.io/malariafwd/reference/): every function, starting with
  the terms they share (`?malariafwd`).
- [Model description](https://shaziaruybal.github.io/malariafwd/MODEL_DESCRIPTION.html): the
  model in full, with its parameters and assumptions.
- [Design decisions](https://shaziaruybal.github.io/malariafwd/DECISIONS.html): why each
  modelling choice was made.
- [Validation](https://shaziaruybal.github.io/malariafwd/VALIDATION.html): how the model is
  checked, and what each check found.
- [Derived metrics](https://shaziaruybal.github.io/malariafwd/METRICS.html): what each number a run
  yields means, and how it is counted.
- [Ongoing and future work](https://shaziaruybal.github.io/malariafwd/ROADMAP.html): the sensitivity
  analyses still needed and the extensions planned.
- The validation scripts are installed with the package, with the results of their last run in
  `RESULTS.md`: list them with `dir(system.file("validation", package = "malariafwd"))` and run
  one with `source()`, for example
  `source(system.file("validation", "audit_events_and_ibd.R", package = "malariafwd"))`.

## Limitations

Worth knowing before drawing a conclusion from a run. The full list, with the reasoning behind each,
is in the [model description](https://shaziaruybal.github.io/malariafwd/MODEL_DESCRIPTION.html). Many of these features are part of ongoing or planned future work. 

* **No immunity.** No age-prevalence gradient, no decline in density or detectability with exposure, no transmission-blocking immunity. This is a large simplifying assumption but deliberate because immunity would confound the genetic signal. Clones also do not compete, and detectability does not depend on parasitaemia. This is a key priority item being worked on actively.
* **No liver-stage delay.** A blood-stage infection begins the day it is inoculated, so
  recurrence intervals run about 9 days early, which matters most when the interval is
  what distinguishes a relapse from a reinfection. This is currently a high priority item being added. 
* **Homogeneous mixing**, with no spatial structure or heterogeneous biting, and no
  post-treatment prophylaxis — so a treated host can be reinfected or relapse the next day. This would overstate recurrence wherever treatment coverage is appreciable. This is a planned item during improvement of the treatment features of the model.
* **The genome is neutral** (mutation rate 0): no selection, no drug resistance, no sequence
  divergence. IBD is exact because it is read from the genealogy rather than inferred. Mutations could be overlaid backwards in time using msprime. 
* **Genotypes are painted, not drawn from a real population.** `simulate_genotypes()` gives each
  founder template and import independent alleles from the frequencies you supply, so there is
  no background relatedness or linkage disequilibrium among them. This is under active development. Genetic outputs in the current model are focused on ground-truth COI and IBD. 
* **Imported parasites are unrelated to local ones.** Imports are drawn as fresh founders, so
  they are maximally distinguishable — any import-detection result is an **upper bound**.

## Citing and getting help

Please open an [issue](https://github.com/shaziaruybal/malariafwd/issues) for bugs and questions,
and get in touch if you are interested in using the package, since it is under active development. Cite malariafwd
with `citation("malariafwd")`; it also lists the papers for SLiM and tskit, which run the model and
record its genealogy as tree sequences. An archived DOI will be added at the first tagged release.

## Development and AI use

Model and package design, parameter choices and calibration were decided by the author, and
each design decision is recorded in the [design
decisions](https://shaziaruybal.github.io/malariafwd/DECISIONS.html). The first
iterations of the SLiM model and basic package functions were written by the author without AI assistance. Anthropic's
Claude (via Claude Code) was then used as a coding assistant to add features, build the
package infrastructure, and write tests and documentation. The author reviewed all changes and is
responsible for this software.

## Funding

This work was supported by the "la Caixa" Foundation (ID 100010434), Junior Leader Fellowship
LCF/BQ/PI25/12100035.
