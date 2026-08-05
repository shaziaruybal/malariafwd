# malariafwd 0.2.0

First packaged version of the *P. vivax* SLiM + tskit simulator (model iteration "v2").

## Model

* **Full 14-chromosome PvP01 genome** (24.2 Mb), modelled as one concatenated linkage group
  with independent assortment at chromosome junctions. `GENOME_MODE` (`"full"` / `"core"`)
  selects the coordinate extent; `CHROMOSOMES` selects a subset. Uniform recombination rate
  7.4e-7/bp. (Previously a single 100 kb chromosome.)
* Biology/logic review fixes: gametocyte gating wired in (`GAMETOCYTE_DELAY`, `P_GAMETOCYTE`;
  default identity); MDA remembers individuals before killing (IBD-safe); seasonal peak now
  lands on `SEASONAL_PHASE`; one relapse episode per activation event; activation-rate clamp.
* `VERBOSE` parameter (default `FALSE`): quiet by default, `TRUE` for full per-event narration.
  Does not affect results or data files.

## R interface

* `vivax_params()` — validated parameter object mirroring the model's defaults.
* `run_sim()` — run one simulation via `processx` (shell-free, robust quoting, provenance).
* `run_grid()` — parallel multi-seed / multi-parameter sweeps via `future`/`furrr`.
* `load_run()`, `load_events()`, `load_state()`, `trees_file()`, `run_params()`.
* `calc_ibd()` / `ibd_index()` — ground-truth IBD from the tree sequence (tskit), with a
  `max_time` option for recent-transmission IBD; vacant null haplosome nodes excluded.
* Metrics: `run_metrics()`, `classify_cases()`, `coi_distribution()`, `relapse_fraction()`,
  `relapses_per_person_year()`, `dormancy_durations()`, `rc_by_origin()`, `chain_sizes()`,
  `relapse_origin_ibd()`, `newinfection_diversity()`, and more.
* Plots: `plot_relapse_km()`, `plot_host_state_prevalence()`, `plot_host_trajectory()`,
  `plot_moi_vs_relapse()`.

See `DECISIONS.md` for the design rationale behind these choices.
