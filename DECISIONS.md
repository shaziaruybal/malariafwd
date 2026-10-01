# Decision log — malariafwd (P. vivax SLiM + tskit simulator)

Internal history of the design decisions behind the model + R package, and *why* we made
them. This is deliberately separate from the standard R-package machinery:

- **`NEWS.md`** (later, in the package) = terse user-facing changelog ("what changed").
- **roxygen / vignettes** = API reference + how-to.
- **`DECISIONS.md`** (this file) = the reasoning history ("*why* we chose X over Y, and what
  we rejected"). This is the part normally lost.

**Conventions.** One entry per decision, newest at the top. Entries are *append-only*: when we
change our mind we add a new entry that marks the old one `Superseded by DNN` — we never delete,
so the trail of what we tried survives. Each entry: Context → Decision → Rationale → Alternatives
rejected → Consequences. Cross-reference model code and the review finding IDs (F1–F15).

Model file: `inst/slim/malariafwd_vivax.slim` (SLiM 5.2). This log starts at the first
full-genome prototype milestone (2026-08); entries written before the package was extracted
refer to paths from the prototype that preceded it.

---

## D16 — Verbosity is a PARAM (VERBOSE), one model file, not two · 2026-08-05 · Accepted
**Context.** The model emits ~30k lines/run of per-event narration. Considered a separate "quiet"
model file vs a runtime switch.
**Decision.** A single model with `defineDefault("VERBOSE", F)`; all per-event/per-tick narration goes
through a `vcat()` helper (`if (VERBOSE) catn(msg)`); GENOME/WARNING/SIMULATION-FINISHED stay on `catn()`
(always shown). 50 of 57 `catn()` converted to `vcat()`.
**Rationale.** One source of truth (two files would drift). Verified: VERBOSE=F → 28 log lines, VERBOSE=T
→ 9226, and events/state/trees are byte-identical between modes (narration only, never results). Not a
speed win (logging isn't the bottleneck) — the point is log-disk noise at sweep scale + a debuggable mode.
**Alternatives rejected.** Two model files (drift); deleting the narration (loses the debug mode).
**Consequences.** Default quiet; set VERBOSE=T to debug a single run.

## D15 — R↔SLiM interface: subprocess + processx + future, tskit via reticulate · 2026-08-05 · Accepted
**Context.** Deciding how the `malariafwd` package runs SLiM from R, esp. for multi-seed/multi-param
sweeps. Prior ad-hoc use of `slimr` felt slow/awkward.
**Decision.** (1) Run `slim` as a **subprocess** (static model template + `-d` param overrides + file
output) — NOT slimr. (2) Launch via **processx** (`run_sim()`): no shell, so no quote-stripping; creates
OUTDIR first; captures exit code/stdout; supports timeouts. (3) Parallelise `run_grid()` via
**future/furrr** (backend-agnostic: local multicore now, SLURM/HPC later by switching `plan()`).
(4) Use **reticulate + tskit ONLY for tree-sequence analysis** (IBD), never to drive SLiM. (5) HPC:
design run_grid on the pluggable future backend but don't build/test cluster tooling yet.
**Rationale (benchmarked, Apple M5).** R launch overhead is negligible (~0.2 s; R-launched run 5.2 s ≈
native). Runtime is dominated by the epidemiological compute, NOT the genome (1-chr ≈ 14-chr ≈ 5.4 s),
NOT logging, NOT R. This IS the slendr architecture (subprocess + tskit-via-reticulate); slimr adds a
stdout-parsing layer that explains the old slowness. Parallel scaling is HARDWARE-capped on Apple
Silicon (4 perf + 6 efficiency cores + macOS scheduling forks onto efficiency cores → ~2.7x ceiling at
8 workers); a homogeneous Linux server/cluster scales near-linearly → the venue for big sweeps.
**Pitfalls the launcher must own (both hit during benchmarking).** (a) SLiM's writeFile won't create
OUTDIR → mkdir first. (b) `system2`/`system` route through /bin/sh which strips Eidos string quotes
(`OUTDIR='...'` → malformed) → processx (no shell) or shQuote.
**Alternatives rejected.** slimr (overhead); parallel::mclapply (Unix-only, no cluster path); base
system2 (shell-quoting fragility); building HPC tooling now (premature).
**Consequences.** Deps: processx, future, furrr (Imports); reticulate (Suggests, IBD only). Also planned:
a VERBOSE flag to cut the ~30k-line/run debug log at scale (disk noise, not a runtime cost).

## D14 — Gametocyte maturation delay is a FIXED value, not a distribution · 2026-08-05 · Accepted
**Context.** `GAMETOCYTE_DELAY` sets the lag from blood-stage start to transmissibility (F1/D6).
Considered drawing a per-infection delay from a distribution over ~0–3 days.
**Decision.** Keep it a single fixed value (default 0; vivax-realistic preset ≈3).
**Rationale.** The delay is short relative to the ~60-day blood-stage infectious period, so its
*variance* has negligible dynamical effect — adding a per-infection draw is stochastic structure that
doesn't move results (parsimony). Vivax's early gametocytes also make the realistic value close to the
identity default.
**Alternatives rejected.** Poisson(mean) / uniform{0..3} per-infection draw — trivial (~2 lines) but
unjustified given the low impact; revisit only if a question specifically needs the delay's tail (e.g.
"early transmitter" fraction, treatment-timing effects).
**Consequences.** One fewer stochastic knob; upgrade path noted if needed.

## D13 — Single uniform recombination rate 7.4e-7 for both genome modes · 2026-08-05 · Accepted
**Context.** D4 paired `full` with the genome-wide rate 1.0e-6 and `core` with 7.4e-7, arguing each
reproduces its scaffold's empirical map length. On reflection this confounds the mode switch with a
rate change.
**Decision.** Use `within_r = 7.4e-7` /bp (euchromatic core rate) for BOTH modes. `GENOME_MODE`
becomes a pure coordinate/extent switch (which sequence we keep). `RECOMB_RATE` still overrides.
**Rationale.** The genome is neutral and the map is UNIFORM (not spatially varying). The genome-wide
1.0e-6 is only higher because it averages in hyper-recombinogenic subtelomeres we are not modelling
spatially; spreading that uniformly would overestimate recombination in the euchromatic core (where
all analyzable genome and IBD signal live). The core/euchromatic rate is the correct representative
uniform rate everywhere. It also keeps `full` and `core` directly comparable (identical per-bp
recombination; differ only in extent), so the coordinate choice can be isolated.
**Alternatives rejected.** Keep two rates (D4) — confounds extent with recombination.
**Supersedes.** The rate half of D4 (the chromosome lengths in D4 stand).
**Consequences.** full ≈ 17.9 Morgans, core ≈ 16.1 (both at 7.4e-7); ~1 crossover/chr/meiosis + free
assortment — still fixes the old ~1-crossover/genome limitation.

## D12 — Keep an ADR-style decision log (this file) · 2026-08-05 · Accepted
**Context.** Building the first full-genome prototype; want to remember what we tried and why,
beyond code comments and the R changelog.
**Decision.** Maintain `DECISIONS.md` as an append-only ADR log at the project root.
**Rationale.** Design rationale is normally lost; ADRs are the standard lightweight way to
preserve it. Keeps `NEWS.md`/roxygen for their standard purposes.
**Alternatives rejected.** (a) Rely on git commit messages — too granular, not narrative.
(b) A `docs/decisions/NNN-*.md` folder — heavier than needed for an internal single-maintainer log.
**Consequences.** New non-trivial decisions get an entry here in the same change.

## D11 — MDA remembers individuals before killing them (F3) · 2026-08-05 · Accepted
**Context.** MDA cleared blood-stage + hypnozoites via `killIndividuals` without a preceding
`treeSeqRememberIndividuals`, unlike every other clearance path.
**Decision.** Add `treeSeqRememberIndividuals(..., permanent=RETAIN_FULL_HISTORY)` before both
MDA kills.
**Rationale.** Under `RETAIN_FULL_HISTORY=T`, un-remembered cleared infections can be simplified
out of the tree, destroying the IBD-to-prior-infection ancestry that WS-E exists to preserve.
**Alternatives rejected.** Leave as-is — silently corrupts IBD in any MDA scenario.
**Consequences.** MDA scenarios now retain full genealogy like natural clearance.

## D10 — Seasonal peak lands on SEASONAL_PHASE (F4: sin→cos) · 2026-08-05 · Accepted
**Context.** Seasonal forcing used `1 + AMP·sin(2π(t−PHASE)/365)`, which peaks ~91 days *after*
`PHASE`, but the param is documented as the day-of-year of the peak.
**Decision.** Use `cos` so the maximum is exactly at `t = PHASE`.
**Rationale.** Least-surprise: the parameter should mean what it says. Anyone setting the peak
to a real transmission-season day was off by a quarter-year.
**Consequences.** Seasonal scenarios shift; document in NEWS as a behavior change.

## D9 — Clamp the hypnozoite activation probability (F7) · 2026-08-05 · Accepted
**Context.** Resurgence scenarios multiply the activation rate by `ALPHA_ACTIVATE_BOOST`; a large
boost can push the effective daily activation probability ≥ 1, breaking the competing-risks draw
(`rbinom` invalid; death denominator `MU/(1−α)` negative/explodes).
**Decision.** Clamp `effective_alpha_activate ≤ 1 − MU_HYPNOZOITE`.
**Rationale.** Keeps both competing-risks probabilities valid (α<1 and MU/(1−α)≤1). Biologically,
you cannot have >100% of a batch activating in one day; the guardrail just enforces that under
extreme reactivation-surge scenarios.
**Consequences.** Extreme relapse-resurgence scenarios saturate gracefully instead of erroring.

## D8 — One relapse EPISODE per activation event (F6) · 2026-08-05 · Accepted
**Context.** When N clonal hypnozoites from one batch activated on the same day, the code spawned
N separate blood-stage Individuals.
**Decision.** Spawn a single blood-stage relapse Individual per (batch, day) activation event;
the batch still loses N members (correct depletion).
**Rationale (biology).** A relapse is one recrudescent blood-stage episode. Co-activating members
of a batch are genetically identical clones (same inoculation) → they merge into one blood-stage
infection wave, not N simultaneous episodes. The old behavior over-counted relapse episodes,
delayed lineage clearance (all N copies had to clear), rolled detection/radical-cure N times
(inflating both), and N-weighted that lineage in the within-host mix (COI/outcrossing).
**Alternatives rejected.** Keep per-member spawning — biologically wrong and double-counts.
**Consequences.** Small in practice (n_activating is almost always 0–1 at ALPHA≈1/170/day), but
relapse counts, detection, and COI are now correct. Different *batches* (heterologous genotypes)
activating the same day still yield distinct episodes.

## D7 — Radical cure models the LIVER stage only (F5) · 2026-08-05 · Accepted
**Context.** `applyRadicalCure` thins the host's hypnozoite reservoir but leaves the blood-stage
infection, which then clears naturally over ~MEAN_BLOOD_STAGE_DURATION.
**Decision.** Keep it liver-stage-only (anti-relapse); add a clarifying comment. Blood-stage
(schizonticidal) treatment is a future addition.
**Rationale (biology).** "Radical cure" specifically = the 8-aminoquinolines (primaquine,
tafenoquine) that kill hypnozoites — the one thing that clears the dormant liver reservoir.
Blood schizonticides (chloroquine/ACT) are a *separate* drug class. Modeling radical cure as the
anti-hypnozoite lever is faithful and isolates the epidemiologically key mechanism (draining the
reservoir that sustains transmission). (MDA, by contrast, models a combined blood+liver regimen.)
**Consequences.** A treated case keeps transmitting from blood-stage for ~weeks — a real feature
of vivax control, not a bug. Add blood-stage treatment when treatment becomes a focus.

## D6 — Gametocyte gating wired in but defaults to identity (F1) · 2026-08-05 · Accepted
**Context.** `GAMETOCYTE_DELAY` and `P_GAMETOCYTE` were declared but never read: humans were
infectious the same tick they were infected, and every active infection was treated as gametocytemic.
**Decision.** Implement gating (`gametocytePool()` gates human→mosquito uptake on
`gametocyte_ready_tick` + a cached per-infection `Bernoulli(P_GAMETOCYTE)` carriage roll), but
**default to identity** (`GAMETOCYTE_DELAY=0`, `P_GAMETOCYTE=1.0`, plus a fast path that consumes
no RNG) so existing calibration is unchanged.
**Rationale (biology).** P. vivax gametocytes appear early (at/near symptom onset — no long
sequestration-maturation like P. falciparum's ~10–12 d), so infectiousness is ~concurrent with
blood-stage infection. The Imperial/White vivax models likewise key infectiousness off blood-stage
state with no long gametocyte lag. Wiring the params in makes them real for sensitivity analysis
while the identity default preserves behavior. Realistic preset: delay≈3, carriage≈0.9.
**Alternatives rejected.** (a) Delete the params — loses a useful sensitivity knob. (b) Bake in the
realistic preset as default — shifts the existing calibration for little gain.
**Consequences.** Empty gametocyte pool → uptake fails gracefully (mosquito not infected).

## D5 — Core-genome lengths use a 90% proxy pending exact Pv4 BED · 2026-08-05 · Accepted (interim)
**Context.** `GENOME_MODE="core"` needs per-chromosome subtelomere-masked lengths. The exact
MalariaGEN Pv4 core-region BED requires the `malariagen_data` package at a current version + Google
Cloud access, which was unavailable in the build environment.
**Decision.** Use `core = round(0.90 × full)` as a documented proxy (≈21.8 Mb, matching the
published ~22 Mb Pv core), structured as a one-line data swap.
**Rationale.** Unblocks the `core` option now; the approximation only affects the non-default path
(full is default, D3) and the exact coordinates drop in trivially later.
**Alternatives rejected.** Block on cloud access; hand-enter per-chromosome subtelomere boundaries
without a citable source.
**Consequences.** `PVP01_CORE` is marked PLACEHOLDER in code; replace with exact Pv4 coordinates
when available (would supersede this entry).

## D4 — Recombination rates from Jiang et al. 2011 · 2026-08-05 · Rate choice SUPERSEDED BY D13 (lengths stand)
**Context.** Need a per-bp crossover rate for the (uniform-within-chromosome) map.
**Decision.** `full` → 1.0×10⁻⁶ /bp (9.6 kb/cM, genome-wide incl. subtelomeres); `core` →
7.4×10⁻⁷ /bp (13.5 kb/cM, subtelomere-excluded). P. falciparum estimates (no comparable published
vivax genome-wide map; falciparum and vivax recombine at similar order).
**Rationale.** Each rate is the self-consistent uniform rate that reproduces the empirical total map
length over its own scaffold extent (full incl. hot subtelomeres → higher rate; core → lower).
Yields ~24 Morgans (full) / ~16 (core) — realistic, ~1 crossover per chromosome per meiosis plus
free assortment, which is what makes length-based IBD time-discriminating (fixes the old
100 kb/1e-5 ≈ 1-crossover limitation).
**Alternatives rejected.** A single rate for both modes — inconsistent with the differing extents.
**Consequences.** `RECOMB_RATE > 0` overrides for sensitivity analysis.

## D3 — Default genome is FULL (complete neutral scaffold) · 2026-08-05 · Accepted
**Context.** `GENOME_MODE` can be `full` (24.2 Mb, exact PvP01 lengths) or `core` (subtelomere-masked).
Initially defaulted to `core` (surveillance realism). User revisited.
**Decision.** Default `full`.
**Rationale.** The model is fully NEUTRAL (mutation rate 0; all genetics from the ground-truth
genealogy). Subtelomere masking exists in real data because those regions are hard to genotype/align
— a concern that does not exist in a simulation with exact ancestry. So the complete physical scaffold
is the honest, complete representation; `core` is only for deliberately mimicking a masked panel.
**Alternatives rejected.** Default `core` — imports a real-data limitation the simulator doesn't have.
**Supersedes.** The initial `core` default (never released).
**Consequences.** Runs cover the full 24.2 Mb (~9 s at prototype scale — negligible cost).

## D2 — 14 chromosomes as ONE concatenated linkage group (Option B) · 2026-08-05 · Accepted
**Context.** Extend the single 100 kb chromosome to all 14 PvP01 chromosomes. SLiM 5.2 supports
explicit multiple chromosomes.
**Decision.** Model the genome as ONE concatenated chromosome with recombination rate 0.5 at each of
the 13 chromosome junctions (independent assortment) and the within-chromosome rate elsewhere.
**Rationale.** Genetically equivalent to explicit chromosomes for every relatedness/IBD metric we
compute (a shared IBD segment cannot cross an r=0.5 junction). Verified empirically: an r=0.5 junction
yields P(breakpoint)=0.5, so the *existing* `addRecombinant` in `seedTetrad` gets independent
assortment for free — the meiosis code and the entire tskit IBD pipeline stay UNCHANGED, and output
remains a single `.trees` file. Validated: clonal IBD=1.000; same-oocyst-sib IBD distribution tightens
toward the center under 14 chromosomes vs bimodal under one (assortment working).
**Alternatives rejected.** Explicit SLiM multiple chromosomes — writes a *directory* of
`chromosome_N.trees` (one per chromosome), forcing a rewrite of `calc_ibd_subset` to aggregate across
14 files and reworking `seedTetrad` into per-chromosome meiosis (review finding F2). Only unique payoff
(per-chromosome tree sequences for chromosome-level selection scans) is not on the roadmap.
**Consequences.** A comment marks the switch-point if explicit chromosomes are ever needed.

## D1 — Package identity & scope: `malariafwd`, engine + IBD toolkit · 2026-08-05 · Accepted
**Context.** Packaging the run+analysis engine as a prototype R package for a new GitHub repo.
**Decision.** Name `malariafwd` (matches the model's internal name; room to extend beyond vivax).
Scope = **engine + IBD toolkit**: params object, `run_sim`/`run_grid`, `load_run`, tskit IBD, the
general metrics + plots from `utility_functions.R`. The ~117 analysis scripts + framework work
(E1/classify/design) stay in the prototype repo as the paper's code. Scaffold as `./malariafwd/`,
port to a fresh repo later. Do the model review + 14-chromosome genome FIRST so the package ships
around a corrected model.
**Rationale.** Ships the reusable engine, not the research; smallest coherent surface for a v0.1.
**Alternatives rejected.** Engine-only (too thin — loses the tested IBD/metrics); engine+framework
(couples the package to in-progress research).
**Consequences.** Package SemVer starts pre-1.0 (0.x); model iteration "v2" ⇔ 0.2.0.
