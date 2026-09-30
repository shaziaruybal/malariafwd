# malariafwd 0.1.0

First public release of **malariafwd** — an R interface to a SLiM (Eidos) + tskit
individual-based forward simulator of *Plasmodium vivax* transmission and genetics.
(Versions 0.1–0.3 during development were internal, pre-release iterations; the public
version series starts here.)

## Model

* **Epidemiology:** SEIS individual-based human–mosquito transmission (Ross–Macdonald),
  simulated day by day; blood-stage infection with exponential clearance; superinfection.
* **Hypnozoite reservoir & relapses:** White (2014) / Mehra batch model — each inoculation
  deposits a batch of hypnozoites, each with a daily activation (relapse) and death hazard;
  optional temperate latency (`K_LATENCY`, `DELTA_LATENCY`).
* **Genetics:** whole 14-chromosome PvP01 genome (24.2 Mb), neutral, with forward mutation
  and meiotic recombination in the mosquito oocyst; a tskit tree sequence is recorded so
  identity-by-descent is **ground truth**, not marker-inferred.
* **Mosquito uptake & inoculation (WS-F):** one blood meal forms `LAMBDA_OOCYSTS` oocysts
  (zero-truncated; `K_OOCYSTS > 0` makes it negative-binomial, var = λ + λ²/K, since vivax
  membrane-feeding oocyst counts are overdispersed), each drawing its two gametes
  independently and uniformly from the gametocyte pool — so P(outcross) = 1 − Σpᵢ² *and*
  the pairing P(i,j) ∝ pᵢpⱼ emerge from the draws rather than being imposed. One infectious
  bite delivers `LAMBDA_SPOROZOITES` products (zero-truncated, weighted so a selfed oocyst
  counts as its four identical sporozoites), collapsed to one genotype per distinct lineage,
  with a single per-bite hypnozoite batch split across the co-inoculated genotypes.
  Zero-*truncation* (not inflation) is used because `P_TRANS_HV`/`P_TRANS_VH` already model
  whether the event succeeds. λ = 0 (the default) gives exactly one oocyst and one genotype.
  `N_OOCYSTS_PER_BITE` is retired and now stops the run if set to anything but 1.
* **Parasite diversity axis:** founder establishment and founder diversity are decoupled
  (`N_FOUNDER_HOSTS` × `N_FOUNDER_LINEAGES`); both default to `N_FOUNDERS`.
* **Interventions & importation:** case-based **radical cure** (liver-stage) and **blood-stage
  treatment** (schizonticidal, `BLOODCURE_COVERAGE` / `P_BLOODCURE_SURVIVAL`), mass drug
  administration, importation, seasonality, and resurgence triggers — all off by default, so
  a bare run is a clean closed system.
* **Truth labels:** a relapse now carries the imported origin of its source hypnozoite, so
  **imported relapses** (a relapse of an imported genotype) are a distinct, directly-labelled
  truth category; `classify_cases()` exposes the four classes via `origin` + `imported`
  (imported new infection / local reinfection / imported relapse / local relapse).

## R interface

* `vivax_params()` builds a validated parameter set; `run_sim()` launches SLiM as a
  shell-free subprocess (`processx`) and returns a run handle; `run_grid()` runs parameter ×
  seed sweeps in parallel via `future`/`furrr`.
* Tidy loaders (`load_events()`, `load_state()`), epidemiological + genetic metrics
  (`summarise_run()`, `classify_cases()`, `incidence_by_origin()`, `rc_by_origin()`,
  reservoir/relapse/COI helpers), ground-truth IBD via tskit (`ibd_pairs()`, `ibd_index()`),
  and plotting helpers.
* `event_types()` documents the `events.tsv` event-type dictionary.
* `uptake_moments()` / `uptake_calibrate()` for the uptake draws. `LAMBDA_*` is a **pre-truncation**
  rate, so it is *not* the mean count reported in a feeding study -- conditioning away the
  zeros pulls the realised mean up (lambda 2 -> 2.31 Poisson, 3.0 at K = 1). `uptake_moments()`
  shows what a setting produces; `uptake_calibrate(mean, var)` converts a target realised mean
  (and variance) into the `lambda` / `K` to set.
* **ARG-native summaries** (`arg_summary()`, `arg_tmrca()`, `arg_diversity()`,
  `arg_divergence()`). The model is neutral (mutation rate 0), so every *site* statistic is
  zero and tskit's **branch** mode is the right tool. Two properties of a forward-simulated
  genealogy are documented and guarded here: it is a **forest** (each founder lineage and
  each import is an independent genome root, so most pairs never coalesce), and on a forest
  **branch divergence is not a relatedness measure** — measured on a real run, pairs that do
  coalesce scored 29.1 and pairs that never coalesce scored 32.1, i.e. indistinguishable,
  because the statistic accumulates branch length up to each sample's own root. Use
  `arg_tmrca()`, which returns the time to the common ancestor *and* `coalesced`, the
  fraction of the genome over which one exists (relapse vs its origin: coalesced 1.00;
  distinct founders: 0.00). Every function reports `coalesced` so the dilution is visible.
* `plot_host_trajectory()` gains `layout = "lineages"` (the default): the host's blood-stage
  lineages stacked in their own colours, with the hypnozoite reservoir on a second axis, so
  individual genotypes can be followed as they appear, clear and return. `layout = "summary"`
  keeps the previous COI-and-reservoir line view.
* Plot functions label rather than conclude. Titles that asserted a result (for example
  "Polyclonal infections predict more relapses") have been removed -- a figure should not
  state a finding that may not hold for the user's own run.
* `plot_diagnostics()` — a one-call run-health dashboard (prevalence, incidence + force of
  blood-stage infection by origin, COI, hypnozoite reservoir, circulating-lineage richness),
  which surfaces extinction / saturation / non-stationarity at a glance; and
  `force_of_infection()` returning FOI (mosquito-borne) and molFOB (compound, incl. relapses;
  Mueller et al. 2013) by origin.

## Correctness fixes (WS-G)

Four of these stem from one Eidos scoping rule: a plain assignment to a global inside a
callback or function creates a *local*, so the write was never visible.

* **Hypnozoite batches never died.** Batches that were exhausted or hit the age cap kept
  activating with their old count — in one check run **64% of batches (74/116) drew more
  hypnozoites than they ever held** (636 excess), inflating relapses. Batches are now retired
  with `batch_count = 0` and killed in the `early()` cleanup. A regression test asserts no
  batch over-draws. **Relapse-related calibration from before this fix needs re-running.**
* **`inoculation_id` restarted every tick** (max id 3 across a whole run); now unique.
* **`event_id` was always `<tick>_1`** (6745 rows → 1073 ids); now unique.
* **Relapse-triggered radical cure was mostly undone** — it thinned the batch, then the
  relapse bookkeeping overwrote it. It now runs after the host's batches are updated.
* **MDA errored** when a treated host was left with no infections; that host is now skipped.
* `ibd_index()`'s lookup returned 0 for a genome paired with **itself** (tskit only stores
  distinct pairs, so self-pairs fell through to the below-detection default) and carried
  node-map names on its return value. It now returns 1 for self-pairs and unnamed values.

An integrity audit of the event log, state snapshots and downstream ARG/IBD lives in
`inst/validation/audit_events_and_ibd.R` (62 checks: event-type coverage, id
uniqueness, hypnozoite and blood-stage conservation against `state.tsv`, oocyst/tetrad
structure, tree retention, the IBD ladder from clonal to unrelated founders, and the
ARG/TMRCA checks). `audit_arg_recombination.R` additionally verifies that breakpoints follow
the recombination map (within-chromosome density CV 0.013; KS vs uniform p = 0.29) and
`audit_arg_junctions.R` that the 13 chromosome junctions assort independently (50.7% same
IBD state either side, phi = -0.002) and explains the tetrad structure behind IBD = 0
between same-oocyst sibs. `explain_ibd.R` renders the segment and relatedness-ladder figures, and `check_parameters.R` runs 34 configurations (35 checks) covering structural invariants, directional parameter response and seed reproducibility. See `inst/validation/README.md`.

## API naming

Names follow the tidyverse conventions: one prefix per family (so autocompleting `ibd_`,
`arg_`, `load_`, `plot_`, `slim_` or `uptake_` finds everything in that family), verbs for
actions and nouns for objects. `run_` now means only "execute" -- the accessors that used it
possessively have moved to the families they belong to (`run_params` -> `load_params`,
`run_diagnostics` -> `plot_diagnostics`, `run_metrics` -> `summarise_run`). Other changes:
`calc_ibd` -> `ibd_pairs`, `relatedness_trajectory` -> `ibd_trajectory`,
`relapse_origin_ibd` -> `ibd_relapse_origin`, `who_classify` -> `classify_who` (matching
`classify_cases`), `compute_time_to_relapse` -> `time_to_relapse`,
`newinfection_diversity` -> `new_infection_diversity`, `malariafwd_events` -> `event_types`,
`pick_interesting_host` -> `pick_host`, `default_model` -> `slim_model`,
`trees_file` -> `trees_path`, `model_generation` -> `model_version`,
`chain_sizes` -> `transmission_chains`, `zt_moments`/`zt_calibrate` ->
`uptake_moments`/`uptake_calibrate`.

## Requirements

SLiM (≥ 5.0) on `PATH`; Python 3 with `tskit` + `pandas` (via `reticulate`) for IBD.
