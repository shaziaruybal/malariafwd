# malariafwd 0.1.0

First public release — a **pre-release**, under active development. The interface and the
model may both change before 1.0.0.

**malariafwd** is an R interface to a forward, individual-based, genotype-resolved simulator of
*Plasmodium vivax* transmission, written in SLiM 5 (Eidos) with tskit tree-sequence recording.
One SLiM individual is one parasite genome, so every relatedness quantity is read from the
recorded genealogy rather than inferred from markers: **IBD here is ground truth**.

## What it models

* **Transmission.** Individual-based human–mosquito transmission (Ross–Macdonald), simulated
  day by day, with superinfection. There is no intrinsic incubation — blood stage starts at
  inoculation — so hosts are not SEI.

* **The hypnozoite reservoir.** A White (2014) / Mehra batch model: each inoculation deposits a
  batch of hypnozoites, each with a daily activation (relapse) and death hazard. Temperate
  latency is optional (`K_LATENCY`, `DELTA_LATENCY`). The reservoir is a second infected
  compartment, invisible to any blood survey.

* **Three-state recurrence.** A recurrence can be a **relapse** (a hypnozoite activates), a
  **reinfection** (a new bite or an importation), or a **recrudescence** (a clone that survived
  a treatment course was suppressed below patency and later re-expanded). The third state is
  what a transmission model usually omits, and all three carry ground-truth labels — so a
  method that tries to tell them apart can be scored against the truth.

* **Case management.** Treatment is decided where a case arises and applied to the **host**
  `TREATMENT_DELAY` days later, clearing every patent blood-stage clone the patient carries.
  Radical cure is **nested inside** a blood-stage course (`RADCURE_COVERAGE` is
  P(radical cure | treated)) and is liver-stage only; its efficacy is per person by default,
  with `RADCURE_MODE = "hypnozoite"` restoring per-hypnozoite thinning — the only mode that can
  partially clear a reservoir. Case management has its own window (`TREATMENT_START` /
  `TREATMENT_END`), separate from vector control. Mass drug administration (`MDA_TICK`) gives one
  round of the same course to each host with probability `MDA_COVERAGE`, symptoms or not, so it
  also reaches subpatent clones and clones that became patent earlier that day.

* **Genetics.** The whole 14-chromosome PvP01 genome (24.2 Mb), neutral (mutation rate 0), with
  meiotic recombination in the mosquito oocyst. Each oocyst is its own meiosis producing a
  tetrad; selfing versus outcrossing **emerges** from drawing two gametes independently from
  the gametocyte pool rather than being imposed.

* **Mosquito uptake.** One blood meal forms `LAMBDA_OOCYSTS` oocysts (zero-truncated;
  `K_OOCYSTS > 0` makes it negative-binomial, since vivax membrane-feeding counts are
  overdispersed). One infectious bite delivers `LAMBDA_SPOROZOITES` products, collapsed to one
  genotype per distinct lineage, with a single per-bite hypnozoite batch split across them.
  `LAMBDA_* = 0` (the default) gives exactly one oocyst and one genotype, which makes
  co-transmission structurally impossible — set them above 0 if you need it.

* **Importation and interventions.** Importation (rate, fixed count or single pulse), vector
  control, seasonality, and resurgence triggers that switch on biting, activation, batch size
  or import rate at a chosen tick. All off by default, so a bare run is a clean closed system.

* **Truth labels.** `classify_cases()` labels every blood-stage clone onset by its true origin:
  `imported`, `local reinfection`, `relapse` or `recrudescence`. A separate `imported` flag
  marks import-derived genotypes, so an **imported relapse** is directly identifiable.
  `episodes()` gives every clone's span, and `recurrences()` pairs each recurrence with the
  earlier episodes of its host, with ground-truth IBD.

## Episodes and acquisitions are different units

Worth knowing before you compute a rate. A recrudescence is a new clinical **episode** but not
a newly **acquired** parasite — the genotype was already in the host. The package keeps the two
apart rather than choosing for you:

* `force_of_infection()` returns a recrudescence as its own row with `is_acquisition = FALSE`.
  Sum only the acquisition rows for FOI or molFOB.
* `relapse_fraction()` reports both: `relapse_fraction` over acquisitions (the quantity the
  published 0.55–0.85 vivax band refers to) and `relapse_frac_episodes` over every patent
  episode.
* `rc_by_origin(unit = "episode" | "acquisition")` credits onward transmission either to the
  recurrence itself or to the infection that originally brought the clone in. The second is the
  unit a reproduction number estimated from routine case counts refers to.

## The R interface

* `vivax_params()` builds a parameter set, kept in sync with the model and checked again before
  each run; `run_sim()` launches SLiM as a shell-free subprocess (`processx`) and
  returns a run handle; `run_grid()` runs parameter × seed sweeps in parallel via
  `future`/`furrr`, reporting any failed cell with its error.

* **Loaders and metrics.** `load_events()`, `load_state()`, `load_params()`; `summarise_run()`,
  `classify_cases()`, `episodes()`, `recurrences()`, `incidence_by_origin()`, `force_of_infection()`, `relapse_fraction()`,
  `rc_by_origin()`, `transmission_links()`, `transmission_chains()`, plus reservoir, relapse,
  COI and dormancy helpers.

* **Ground-truth IBD.** `ibd_pairs()` computes genome-wide IBD from the tree sequence and, with
  `segments = TRUE`, returns the shared blocks themselves — each tagged with the common ancestor
  it came from, keyed by pedigree id: `ibd_lookup()`, `as.matrix()` and `as.data.frame()` read it
  (`NA` for a genome the run's tree sequence does not hold); `ibd_segments()` and
  `plot_ibd_pair()` show *why* a pair's IBD is what it is, painted along the 14 chromosomes.

* **ARG-native summaries.** `arg_summary()`, `arg_tmrca()`, `arg_diversity()`,
  `arg_divergence()`. The model is neutral, so every *site* statistic is zero and tskit's
  **branch** mode is the right tool. Two properties of a forward-simulated genealogy are
  documented and guarded here: it is a **forest** — each founder lineage and each import is an
  independent root, so most pairs never coalesce — and on a forest **branch divergence is not a
  relatedness measure** (measured: pairs that do coalesce scored 29.1, pairs that never coalesce
  32.1, i.e. indistinguishable). Use `arg_tmrca()`, which reports a genealogical **depth**
  together with `coalesced`, the fraction of the genome over which a common ancestor exists.

* **Plots.** `plot_diagnostics()` is a one-call run-health dashboard (prevalence, incidence and
  force of infection by origin, COI, reservoir, lineage richness) that surfaces extinction,
  saturation and non-stationarity at a glance. `plot_host_trajectory()` follows one host's
  lineages as they appear, clear and return. Plots label rather than conclude: no title asserts
  a finding that may not hold for your run.

* **The event dictionary.** `event_types()` is the authority on the event types in
  `events.tsv`, with `ends_episode` and `is_onset` flags. Derive filters from it rather than
  restating event names — that is what the flags are for.

* **Calibrating the uptake draws.** `LAMBDA_*` is a **pre-truncation** rate, so it is not the
  mean count a feeding study reports — conditioning away the zeros pulls the realised mean up
  (λ = 2 gives 2.31 under Poisson, 3.0 at K = 1). `uptake_moments()` shows what a setting
  produces; `uptake_calibrate(mean, var)` converts a target realised mean and variance into the
  `lambda` and `K` to set.

## Validation

Three suites ship with the package and run against a live model:

* `inst/validation/audit_events_and_ibd.R` — 99 checks over the event log, state snapshots and
  the downstream ARG/IBD: event-type coverage, id uniqueness, hypnozoite and blood-stage
  conservation against `state.tsv`, oocyst and tetrad structure, tree retention, and the IBD
  ladder from clonal pairs down to unrelated founders.
* `inst/validation/check_parameters.R` — 40 checks over 34 configurations: structural
  invariants, directional parameter response, and seed reproducibility.
* `audit_arg_recombination.R` and `audit_arg_junctions.R` verify that recombination breakpoints
  follow the map and that the 13 chromosome junctions assort independently.

## Status and limitations

This is a pre-release. The epidemiological and genetic core is validated and calibrated against
published vivax bands. The treatment and recurrence machinery is implemented and tested, and
calibration of it is ongoing.

* **Everything beyond the core is off by default.** Treatment, MDA, importation, seasonality and
  resurgence are all disabled, so a bare run is an untreated closed system. In particular
  `P_RECRUD_GIVEN_FAIL = 0`, so the third recurrence state is inert unless you switch it on —
  which matches the evidence that recrudescence is a fraction of a percent of recurrences after
  chloroquine or artesunate, and material only where resistance is present.
* **Imports are unrelated to local parasites**, which makes them maximally distinguishable and
  so gives an upper bound on import-detection performance. A structured source population is
  planned.
* **No within-host densities**, so no competition between clones and no density-dependent
  detectability. Infectiousness is a per-clone property, not a function of parasitaemia.
* **No intrinsic incubation** and a fixed EIP; **homogeneous mixing**, with no spatial structure
  and no heterogeneous biting.
* **Detection is not origin-dependent** beyond the primary / relapse / recrudescence split, so
  active and passive surveillance cannot currently be distinguished.
* `NU_HYPNOZOITE_BATCH` (3) is fitted to one setting — re-fit it per setting against relapse rate
  and relapse fraction.

### In development

Work in progress for the next release:

* a **structured source population** for imported infections, so that imports carry real
  ancestry and can be related to one another (they currently enter as unrelated founders);
* **post-treatment prophylaxis**, following the malariasimulation protection curves;
* **more extensive calibration of the treatment and recurrence machinery**, including a
  chloroquine-resistant setting in which recrudescence is material rather than negligible.

`MODEL_DESCRIPTION.md` gives the full model with its assumptions, parameter tables and
validation status; `DECISIONS.md` records why each modelling choice was made.

## Requirements

SLiM (≥ 5.2) on `PATH`; for the IBD and ARG functions, Python with `tskit` (≥ 1.0), `pandas` and
`numpy`, through `reticulate`. `check_setup()` reports both. Everything else is R.
