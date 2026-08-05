# malariafwd

Forward simulation of *Plasmodium vivax* transmission and genetics — an R interface to a
[SLiM](https://messerlab.org/slim/) (Eidos) + [tskit](https://tskit.dev/) individual-based
model. The model simulates blood-stage infection, the **hypnozoite reservoir and relapses**,
meiotic recombination in the mosquito, and importation/interventions on a neutral
14-chromosome PvP01 genome, recording a tree sequence for **ground-truth identity-by-descent**.

## Requirements

- **SLiM ≥ 5.0** on your `PATH` (or set `MALARIAFWD_SLIM`). Check with `slim_available()`.
- For IBD: **Python 3 with `tskit` + `pandas`** (via the `reticulate` package).

## Install

```r
# install.packages("remotes")
remotes::install_local("malariafwd")   # or install_github("<user>/malariafwd") once ported
```

## Quick start

```r
library(malariafwd)
slim_version()                                    # e.g. "5.2"

# one run
p   <- vivax_params(N_HOSTS = 300, BITING_RATE = 0.09)   # 52 validated params; full genome by default
run <- run_sim(p, seed = 1, outdir = "runs/demo")

# tidy outputs
ev  <- load_events(run)
st  <- load_state(run)
run_metrics(run)                                  # prevalence, COI, relapse fraction, ...

# ground-truth IBD (tskit)
cs  <- classify_cases(ev)
idx <- ibd_index(calc_ibd(run, ids = cs$pedigree_id))
relapse_origin_ibd(ev, idx)                       # relapses are clonal to their origin (IBD = 1)
newinfection_diversity(run)                       # relatedness of the transmitting population

# parallel sweep (see ?run_grid; use future::plan() for a cluster)
res <- run_grid(
  list(low = vivax_params(BITING_RATE = 0.06),
       hi  = vivax_params(BITING_RATE = 0.12)),
  seeds = 1:10, outdir_base = "runs/sweep", workers = 8)
```

## How it works

The model runs as a **subprocess** (`processx`, shell-free) with parameters passed as `-d`
overrides; it writes `events.tsv`, `state.tsv`, `plasmo.trees`, and a provenance snapshot into
each run directory. `reticulate` + `tskit` are used only to compute IBD from the recorded tree
sequence. Multi-run sweeps parallelise via `future`/`furrr`. Set `VERBOSE = TRUE` to see the
full per-event narration; the default is quiet.

## Design history

The genome model, biology fixes, and interface choices are recorded as an append-only decision
log in [`DECISIONS.md`](../DECISIONS.md) (the project root).
