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

## D23 — Metrics derive their vocabulary from classify_cases(), and output schemas do not depend on run content · 2026-10-22 · Accepted
**Context.** A verification pass over the four remaining downstream defects found them to be one
bug class in two forms. **Vocabulary drift:** `new_infection_diversity()` matched
`event_type == "vh_inoculation"` only, so imported genomes were excluded from the package's
headline genetic-surveillance metric -- at a 17% import share that overstates mean IBD and clonal
fraction by ~44% and understates lineage richness by 40%, rising to 4.3x on mean IBD at a 61%
share, with all three biases pointing the same way as a genuine transmission decline.
`relapse_fraction()` likewise excluded imported acquisitions (+6.4% relative), but that turned
out to be the smallest of three problems with it: three different definitions of "relapse
fraction" were shipping in one project (this function, `epi_summary()` in the analyses repo, and
`MODEL_DESCRIPTION.md` §16.3), and the function's denominator and its title disagreed.
`transmission_links()` filtered `parent1_pedigree_id > 0` where the no-parent sentinel is -1 and
pedigree 0 is a real genome -- the first seed founder -- costing 8 of 1769 edges across the audit
runs and, on `treat`, an entire four-case transmission chain (7 components reported against 8).
**Content-dependent schemas:** `classify_who()` returned 9 columns on import-free runs and 6 on
the import run, with `imported` -- the truth flag -- present exactly where it is trivially FALSE
and absent on the one run where it carries signal. `ibd_segments()` was worse: the empty-result
path returned the same columns in a *different order*, and `rbindlist`'s default
`use.names = "check"` binds by position with only a message, pouring `gstart` values into `chrom`
(1021664 against 1). `ibd_pairs()` returned logical columns when a run had no shared ancestry.
**Decision.** Metrics derive their case vocabulary from `classify_cases()` instead of restating
event types. `relapse_fraction()` becomes relapses over **all** classified blood-stage episodes,
matching `epi_summary()`; its second column is renamed `n_reinfection` -> `n_nonrelapse` to stop
the name asserting something the denominator no longer means. `transmission_links()` tests
`!= -1`. `classify_who()` groups by whatever `classify_cases()` returns. `ibd_segments()` and
`ibd_pairs()` return a fixed column order and fixed column classes on every path. Both capped
subsamples get a `seed` argument.
**Rationale.** Deriving from `classify_cases()` is what makes this the last time: when Step 4 adds
a `recrudescence` origin, every one of these metrics picks it up with no further edit, which is
the point given that Step 4 is in scope for the pre-release. Aligning `relapse_fraction()` to
`epi_summary()` matters because that is the definition the published 0.55-0.85 band was scored
against, so the package and the calibration evidence now agree. The schema fixes are about
pooling: a grid sweep rbinds per-run tables, and a schema that varies with run *content* rather
than with arguments either errors on one cell or -- in the `ibd_segments` case -- silently
produces numerically plausible garbage.
**Alternatives rejected.** *Strict recurrences-only `relapse_fraction()`* -- defensible, and it
matches the old title, but it would disagree with the analyses repo and with the band the model
was calibrated against; it is documented in the roxygen as a one-liner instead, with both numbers
(0.351 against 0.684 on a heavily treated run) so nobody compares them by accident.
*Deterministic head-of-list subsampling* -- reproducible but biased toward early genomes; a fixed
local seed keeps the draw unbiased and restores the caller's RNG stream.
*Also retire `force_of_infection(split_imported = "auto")`* -- it is the same content-dependent
class, but changing that default is user-visible behaviour rather than a schema fix, and
`plot_force_of_infection()` already handles both vocabularies after D21. Deferred deliberately.
*Also drop `is_dormant == 0` from `classify_who()`'s import predicate* -- a dormant import
establishment is the travel marker for a later imported relapse, and the function disagrees with
`classify_cases()` on 239 cases because of it. That is a behaviour change, not a schema fix, and
needs its own decision.
**Consequences.** `relapse_fraction()` changes value on every run (`imports` 0.745 -> 0.700) and
renames one column; `summarise_run()$relapse_frac` moves with it. `new_infection_diversity()`
changes value on any run with imports. `transmission_chains()` gains components on three of seven
runs. `rc_by_origin()` does **not** move -- it is structurally immune because `classify_cases()`
excludes seed founders, so a test asserting a changed Rc here would be wrong. Fixing
`new_infection_diversity()` grew the `imports` pool past `cap = 250`, which exposed the unseeded
`sample()` and turned a previously deterministic call random (0.0529 / 0.0544 / 0.0526 / 0.0523
across four identical calls) -- hence the seed argument, added to `ibd_trajectory()` too.
**How to verify.** `classify_who()` must return one identical schema across all seven audit runs
and pool with `rbindlist` (10 columns, 4572 rows, 0 NA in `imported`, 518 "Imported" on the import
run) with row counts equal to `classify_cases()`. `relapse_fraction()` must equal
`mean(classify_cases(ev)$origin == "relapse")` exactly on every run. `transmission_links()` must
return >= 445/24/187/430/299/61/323 edges with no -1 source. `ibd_segments()` must return the same
column order whether a pair shares everything or nothing, and a `rbindlist` of both must keep
`chrom` within 1-14. `ibd_pairs()$ibd` must be integer/integer/numeric. Repeat calls to
`new_infection_diversity()` and `ibd_trajectory()` must agree, and must leave `.Random.seed`
unchanged. Suites: `test_dir()` 71/0, `audit_events_and_ibd.R` 66/0, `check_parameters.R` 41/0,
`R CMD check --as-cran --no-manual` 2 NOTEs.

---

## D22 — IBD counts all ancestry; recency is the caller's filter over returned segments · 2026-10-21 · Accepted
**Context.** `ibd_pairs(max_time=)` was documented as counting "only ancestry coalescing within
this many time units (days)" -- a genealogical **depth**. It passed the value straight to tskit's
`ibd_segments(max_time=)`, which compares against the common ancestor's **absolute node time**,
and SLiM writes node times backwards from the final tick (`time_units` = "ticks", end tick 400,
sample node times spread 0-399 with only 0.3% at zero). So `max_time = 180` meant "keep only
ancestry from the last 181 ticks" -- a window pinned to the end of the run, unrelated to the pair
being asked about. Measured on 1101 relapse-origin pairs that are clonal by construction (true IBD
1.000 for all 1101): the rule `segment kept iff mrca_time <= max_time` reproduces the function's
output for 1101/1101 pairs at six cutoffs, max |diff| 0.00e+00. Results are perfectly bimodal --
every pair reports exactly 1.000 or exactly 0.000, never between -- so the reported "mean IBD" is
a keep-rate, which is why it swings 0.00-0.43 across the seven audit runs at the same setting. It
is also biased, not merely noisy: at `max_time = 180` the surviving pairs have median relapse
latency 38 days and the erased ones 87 days, because a long-dormant relapse has an old common
ancestor by definition. That inverts the comparison the measure exists to make -- scoring each
case by its max IBD to prior cases population-wide, the ordering flips from relapse > local
reinfection to local reinfection > relapse in all three runs tested, and AUC falls 0.622 -> 0.511,
0.547 -> 0.436, 0.540 -> 0.416 (two of three below 0.5, i.e. anti-predictive). Separately,
`arg_tmrca()` returned `tree.tmrca()` raw while promising depth, so 60 relapse-vs-origin pairs
with an identical clonal relationship spanned 43-391 ticks; its `a == b` branch already returned
0.0, i.e. the function reported depth for self-pairs and absolute time for everything else.
**Decision.** Remove `max_time` from `ibd_pairs()`. All shared ancestry is counted, always.
Replace it with `segments = FALSE`, which when `TRUE` returns every shared segment
(`node_a`, `node_b`, `left`, `right`, `mrca_node`, `mrca_time`) plus `sequence_length`, so a
caller who wants a recency definition writes it as an explicit filter on
`mrca_time - max(t_a, t_b)`. `arg_tmrca()` now returns `tmrca` as a depth from the older sample
and the raw node time alongside it as `tmrca_node_time`.
**Rationale.** The argument had no justified use left. Its performance rationale (tskit's own
docstring recommends a max time for large tree sequences) turned out to be misattributed: the
cost in this code was `store_segments=True`, which was hard-coded on while the result was only
ever summed and discarded. tskit exposes `total_span` per pair without it -- 5.1s against 38.1s
over 14.8M segments on a 5413-sample tree, a 7.5x saving for an identical answer -- and `ids=`
was always the bigger lever. Its scientific rationale was a depth concept the implementation did
not provide; the cited authority (`DECISIONS D10d`) does not exist in this file or in any commit
of it, and the actual record, the prototype's `GENOMIC_SURVEILLANCE_DESIGN_PLAN.md` §10d, says
"coalescing **within** 180 days" and "coalescence-recency" -- depth intent throughout. Making the
segments the unit means the time convention is visible in the analysis code that chooses it rather
than hidden behind an argument whose meaning no caller could see.
**Alternatives rejected.** *Rename to `max_node_time` and keep it* -- honest, but it preserves a
footgun with no demonstrated use; anyone who genuinely wants an end-anchored window can filter the
returned segments. *Add `max_depth` alongside* -- it presumes one depth convention when the data
do not support a single obvious one: 33.8% of segments join members sampled >= 100 ticks apart and
0.3% have the older member *as* the common ancestor, where depth-from-older is 0 and
depth-from-younger is not. *Return segments through a lazy accessor* -- would tie an R object's
validity to a live Python session, breaking `saveRDS()` and `run_grid()` workers, to avoid a
row count that only arises from a call (`segments = TRUE` with no `ids`) nobody makes.
*Fix the documentation only* -- leaves the inversion above in place.
**Consequences.** Breaking change to two exported functions, taken deliberately before any
release (no tag, no published version). Every `ibd_pairs()` call that did not pass `max_time` --
which is every call inside the package, and 100 of the 102 IBD call sites across the three repos --
returns bit-identical numbers and runs faster. The two that did: `vignettes/malariafwd.Rmd` is
rewritten to demonstrate `segments` and to state that `mrca_time` is not a depth, and
`vivax_slim_prototype/e1_ibddef.R` needs its `ibd_time` arm re-specified -- that is the analysis
this defect was corrupting, and is a separate decision. `segments = TRUE` costs 3-7x the default,
growing with sample count. Also closes the `extract_segments` bridge on
`resources/PACKAGE_UTILS_WISHLIST.md`, and removes the dangling `D10d` cross-reference.
**How to verify.** `ibd_pairs(run, ids)` must equal the pre-change no-cutoff output bit for bit
(1101 relapse-origin pairs, max |diff| 0.00e+00, all 1.000); `segments = TRUE` must reconstruct
`ibd_fraction` from `sum(right - left) / sequence_length` for every pair (25,193 pairs / 105,918
segments, max |diff| 5.6e-17); `arg_tmrca()` on one fixed relationship must return one number
regardless of sampling date (60 relapse-vs-origin pairs, `tmrca` range 0-0, `tmrca_node_time`
range 43-391) and `tmrca` must never be negative. Both are now regression tests in
`test-model-invariants.R`. Suites: `test_dir()` 71/0, `audit_events_and_ibd.R` 66/0,
`check_parameters.R` 41/0, `R CMD check --as-cran --no-manual` 2 NOTEs.

---

## D21 — Windowed series report the days they actually cover; onset filters derive from the event dictionary · 2026-10-20 · Accepted
**Context.** The downstream sweep found four defects in the summary plots and in
`force_of_infection()`, all of them "the function is internally correct but silently disagrees
with the event log". Measured before the change, on the audit run set: `plot_diagnostics()`
showed a stacked incidence top of 5 where the true count was 2 and a plateau of 3 across windows
whose true count was 1 (7 of 24 window x origin combinations had no row, and `geom_area` joins
consecutive ROWS, interpolating straight across the gap); `force_of_infection()` divided an
11-day trailing window by a full 30-day bin, understating the final point ~3x and ending every
series in a crash that was pure binning artefact; the `layout = "summary"` onset rug drew 57
marks against 63 real blood-stage onsets on one host, omitting the imported and seed routes; and
episode bands kept `end == start`, drawing 5 phantom one-day episodes out of 10 on a host where
treatment is same-tick.
**Decision.** (1) Zero-fill the `(window x origin)` grid in `plot_diagnostics()` before stacking.
(2) Give `force_of_infection()` a `days` column holding the days each window actually covers, and
divide by `days` rather than `bin`; take the run end from the **state table**, not the last case.
(3) Filter episode onsets on the full inoculation vocabulary with `is_dormant == 0` applied across
**all** routes. (4) Require `end > start` for an episode band.
**Rationale.** Each defect misleads hardest in the regime the panel exists to interrogate. The
interpolation only bites on sparse, suppressed runs — exactly when the question is "did
transmission stop?". The trailing-window crash lands on the last point of every series, where a
reader looks for a trend. Adding `days` rather than dropping the partial window keeps the
estimate (it is unbiased, just noisier) and makes the short window **visible** to callers instead
of silently absent; every full window stays bit-identical, because there `days == bin`. The run
end must come from the state table rather than `max(case tick)`: on the `hypmode` audit run the
last case is at tick 373 but the run reaches 400, so a case-derived end would discard a window
whose true value is **zero** — the most informative reading there is.
**Alternatives rejected.** *Drop the partial trailing window* (which is what `plot_diagnostics()`
does via its `keep()` helper) — consistent, but it throws away a real estimate and, with a
case-derived run end, throws away true trailing zeros. *Leave `force_of_infection()` alone and fix
only the plot* — the metric function is the one analyses call directly. *Keep the rug's dormant
exclusion keyed to `vh_inoculation`* — this was my first attempt and it reintroduced the same
vocabulary-drift bug in the other direction: dormant **import** inoculations leaked through as
onsets, 68 marks against 63 real.
**Consequences.** `force_of_infection()` gains a column; its existing columns and all full-window
rates are unchanged, so no analysis output moves except the final point of a series on a run whose
length is not a multiple of `bin`. The `@return` roxygen is updated. Also clears the long-standing
"no visible binding" NOTE by declaring `series`, `grp` and `band` in `zzz.R`, so the three
undefined globals in `plots.R` no longer appear in `R CMD check`.
**How to verify.** For any run: the `plot_diagnostics()` incidence panel's stacked top must equal
`classify_cases()` counted per window (0 difference on `treat`, `baseline`, `hypmode`); every
`force_of_infection()` row must satisfy `rate == n / N_HOSTS / (days/365)` with `days == bin` for
all but the last window; the summary rug's mark count must equal the host's `is_dormant == 0`
onsets over all four inoculation routes (63 of 63); and no episode band may have `end <= start`.
Suites: `audit_events_and_ibd.R` 66/0, `check_parameters.R` 41/0, `test_dir()` 62/0,
`R CMD check --as-cran --no-manual` 2 NOTEs (new submission + `MODEL_DESCRIPTION.md`).

---

## D20 — Treatment is decided where a case arises and applied later, to the host · 2026-10-20 · Accepted
**Context.** Six audit findings shared one cause: treatment was a flag set on a single new
lineage inside `reproduction()`. Measured on the pre-change model: case-based radical cure
never reached the hypnozoite batch seeded by the bite that caused the case (46 of 87 import
batches and 28 of 64 local batches still relapsed under `RADCURE_COVERAGE = 1`,
`P_RAD_SURVIVAL = 0`), because SLiM does not make a juvenile visible until the callback ends;
**0 of 180** detected relapses ever received a schizonticide, because the flag was only set in
the inoculation path; a course cleared one clone rather than the patient; detection was drawn
per genotype, so one bite could be "detected and treated" for one of its genotypes and not the
other; radical cure could be given without a schizonticide; and clearance events hard-coded
`is_imported = 0`.

**Decision.** Split treatment into DECIDE and APPLY. `decideTreatment(h)` runs where a case
arises (a bite, an import, a relapse), draws detection ONCE PER EPISODE, draws blood-stage
treatment and — nested inside it — radical cure, and schedules the course `TREATMENT_DELAY`
days out (geometric, mean 7, White 2014). `applyTreatment(h, give_liver, src)` runs in `early()`
and hits every blood-stage clone in the host plus, if given, the whole reservoir.
`RADCURE_COVERAGE` becomes P(radical cure | treated). `P_RAD_SURVIVAL` becomes per PERSON,
default 0.29 (~71% cure, White 2018 / malariasimulation); `RADCURE_MODE = "hypnozoite"` keeps
the per-hypnozoite thinning as an experimental mode. MDA routes through the same
`applyTreatment()`.

**Rationale.** `early()` is the first point at which juveniles created in `reproduction()` are
visible, so applying there is what lets a course reach its own batch. Treating the host rather
than a lineage is what every reference model does and what drugs do. Nesting matches clinical
practice — primaquine is not given without a schizonticide. Keeping per-hypnozoite behind a flag
preserves the only mechanism by which radical cure can PARTIALLY clear a reservoir, which is
what reshapes post-treatment relapse diversity and is a capability a genetics model should not
lose.

**Alternatives rejected.** *Keep the per-lineage flag and patch the radical-cure timing* —
leaves the one-clone and untreated-relapse bugs. *Delete per-hypnozoite thinning outright* (the
original plan) — discards the partial-clearance capability; keeping it behind a tested flag
costs ~20 lines. *Carry `P_RAD_SURVIVAL = 0.05` across to the new semantics* — the same number
means 77% efficacy per hypnozoite but 95% per person, silently strengthening radical cure and
moving AWAY from the reference models.

**Consequences.** BREAKING: `RADCURE_COVERAGE` is now conditional on treatment, and
`P_RAD_SURVIVAL`'s meaning and default both change. Treated cases gain an infectious period
they did not have, so case management now suppresses transmission harder — measured: a
permanently-open treatment window can produce FEWER courses than a late one, because it
prevents the cases it would have treated. The case-management window gates the DECISION, not
delivery: a case detected just before `TREATMENT_END` still completes its course afterwards.
MDA no longer emits `mda_hypnozoite_survived`; seven new event types are added.

**How to verify.** With `BLOODCURE_COVERAGE = RADCURE_COVERAGE = 1`, `P_RAD_SURVIVAL = 0`,
`TREATMENT_DELAY = 0`: **zero** `hypnozoite_activation` events. 29/29 detected relapses cleared
by `case_bloodstage_clear`. Under superinfection, 79% of courses clear more than one clone
(mean 3.84, max 21) and 0 of 2,388 treated hosts retain an active clone. `TREATMENT_DELAY` 0
and 7 give mean lags of 0.0 and 7.0 days. `RADCURE_COVERAGE = 1` with `BLOODCURE_COVERAGE = 0`
gives no radical cure and an R-side warning. Treatment-free behaviour is unchanged against the
Step 0 baseline across 5 seeds (prevalence +0.004, reservoir −0.14, relapse fraction −0.002).

---

## D19 — Case management gets its own window, separate from vector control · 2026-10-19 · Accepted
**Context.** Every case-based treatment was gated on `INTERVENTION_START`/`INTERVENTION_END`
(five sites), which is the SAME window that scales biting and the importation rate by
`INTERVENTION_EFFICACY`. Two consequences, both silent: setting `RADCURE_COVERAGE` or
`BLOODCURE_COVERAGE` without also opening that window produced zero treatments and no warning
(this is how the issue was found), and switching treatment on forced a change in biting. The
two are unrelated programmes — case management is drugs given to diagnosed patients; vector
control is bednets and spraying — and nothing should make one imply the other.

**Decision.** Add `TREATMENT_START` / `TREATMENT_END`, defaulting to the whole run, and gate all
case-based treatment on them via a single `caseManagementOn()` helper. `INTERVENTION_*` keeps
biting and importation and no longer touches treatment. `validate_params()` refuses a coverage
set against a closed treatment window.

**Rationale.** Coverage becomes the switch, which is what a reader expects: set
`BLOODCURE_COVERAGE` and cases are treated. Expressing the gate once means the five sites
cannot drift apart. Defaulting to "always on" changes behaviour only for runs that set a
coverage, which is a deliberate, declared break rather than a silent one.

**Alternatives rejected.** *Leave the window shared and warn in `validate_params()` when a
coverage is set without one* — keeps an entanglement that has no modelling justification and
still forces a biting change on anyone who wants treatment. *Reuse `INTERVENTION_*` but add a
separate efficacy for treatment* — does not fix the forced coupling of the two windows' timing.

**Consequences.** BREAKING for any run that set a coverage AND opened the intervention window:
treatment now applies for the whole run unless `TREATMENT_*` is set.
`analysis/11_china_scenarios.R` is the only such analysis and needs `TREATMENT_START = P1`,
`TREATMENT_END = END` added to keep its current meaning; otherwise radical cure would act over
the pre-control baseline too and drain the reservoir before the programme starts. Expect a
stronger radical-cure effect once the window is scoped cleanly, and watch the relapse-driven
resurgence arm, which the script's own comment notes becomes degenerate if radical cure is too
effective. NOT addressed here: `INTERVENTION_EFFICACY` also suppresses importation, so border
measures and bednets remain the same knob — a modelling question, not a bug, and out of scope.

**How to verify.** A coverage alone produces treatments; opening the intervention window alone
produces none while changing prevalence; a bounded `TREATMENT_*` window produces fewer
treatments than an unbounded one; a coverage against a closed window errors. All four in
`tests/testthat/test-treatment.R`.

---

## D17 — Imports seeded in `early()` via a scratch subpopulation; second compartment cloned one tick later · 2026-10-19 · Accepted
**Context.** Two defects, one root cause. (a) Importation could never fire into an EMPTY parasite
population: imports are scheduled in `early()` but seeded inside `reproduction()`, and SLiM calls
`reproduction()` once per existing individual, so with zero parasites the queue accumulated
forever — exit 0, no events. Measured: founders present + `IMPORT_TICK=5` → 10 imports; no
founders + the same setting → 0 events. This made the documented naive-population mode
(`IMPORT_TICK = 1`) dead since it was written, and meant **no run reaching elimination could ever
be re-seeded** — i.e. the decline/resurgence scenarios. (b) An imported inoculation built its
blood-stage infection and its hypnozoite batch with two separate `addEmpty(NULL)` calls, giving
two unrelated ARG roots under one `lineage_id` (`IMPORT_INOCULATION_BUG.md`). Measured on
`mosaic_demo_s1`: IBD(primary, own batch) = 1.000 for local bites (100% of 40) but **0.000 for
imports (100% of 40)**; a relapse of an imported infection reads IBD 0.000 to the infection that
caused it, i.e. it looks like a brand-new unrelated infection.

**Decision.** Seed imports in `early()` by minting the genome in a disposable scratch
subpopulation (`sim.addSubpop`), moving it into the host with `takeMigrants()`, and removing the
scratch in the same callback. Roll the compartment outcomes FIRST; the migrated genome takes
whichever compartment exists (blood-stage if `P_PRIMARY` succeeds, else the hypnozoite batch). If
both are needed, queue the second and clone it from the first in `reproduction()` the following
tick, when the genome is no longer a "new juvenile".

**Rationale.** `addEmpty()`/`addCloned()`/`addRecombinant()` are `reproduction()`-only — verified:
`addEmpty() may only be called from a reproduction() callback`. `addSubpop()` is the only way to
mint a genome outside it; `takeMigrants()` is legal anywhere and preserves `pedigreeID`, `tag` and
all `setValue()` state (verified). A freshly minted genome cannot be cloned the same tick
(`parent must be visible in a subpopulation`), and two `addSubpop` genomes are independent roots,
so cloning next tick is the only route to IBD 1. This is the idiom the model already uses for
founder lineages (`sim.addSubpop(s, 1)`).

**Alternatives rejected.** *Seed in `early()` with `addEmpty`* — impossible, hard SLiM error.
*A permanent sentinel individual* so `reproduction()` always fires — it is alive at the end of the
run and so becomes a tree-sequence sample with a `pedigree_id`, silently contaminating every
all-samples analysis, and needs guards through the hot `reproduction()` path. *A phantom
"inoculum" genome cloned into both compartments*, restoring exact sibling symmetry with a local
bite — costs a one-day delay on every imported blood-stage infection plus a transient ghost
individual, to buy a topology that arrives for free under D18.

**Consequences.** Event-log IDs are UNCHANGED: `parent1_pedigree_id` stays `-1` for both import
rows (neither has a source genome inside the simulation), and `inoculation_id` / `lineage_id`
already carry "same inoculation, same lineage". The ARG gains a parent→child edge from primary to
batch; this divergence from the event log is deliberate, so nothing reconstructing transmission
from `parent1_pedigree_id` infers an intra-host blood→liver transmission that never happened.
Sampled ancestors are not a new class — 10.7% of sample nodes in `mosaic_demo_s1` already have
children, and 1.7% of edges already span ≤ 1 tick. Subpopulation ids cannot be reused
(`subpopulation p99 has been used already`), so a monotonic counter burns one id per import-tick
(~730 over a 40-year run against a base of ~1,400). Imported genomes carry the scratch population
id in the tskit population table — nothing in either repo reads tskit population ids (grepped),
but it is a latent trap for any future analysis mapping node → population → host. The two
`vh_inoculation_import` rows for one inoculation now land on CONSECUTIVE ticks sharing one
`inoculation_id`. Figures on import-bearing runs need regenerating: figure 32/32b's "imported
founder (arrived day 1010)" is pedigree 48469, which has `is_dormant = 1` — one of the phantom
batch roots this decision removes.

**How to verify.** (1) A no-founder, import-only run produces `vh_inoculation_import` events and a
non-empty event log. (2) A run driven to extinction, then given an import window, re-seeds.
(3) `IBD(primary, its own batch) == 1` for imported inoculations, matching local ones — this
un-skips the regression test in `tests/testthat/test-model-invariants.R`. (4) A queued second
compartment whose primary dies before the next tick is purged and the queue does not grow.

---

## D18 — Imports are independent founders for v0.1.0; a source deme is deferred to v0.2.0 · 2026-10-19 · Accepted (interim)
**Context.** Every import enters as an ARG root with no ancestry — exactly the limit of an
infinitely large, infinitely diverse source population, in which two imports can never be related
to each other. The population-genetics convention for immigration is a source ("ghost") deme with
migration, under which imported genomes carry real ancestry and coalesce with one another in the
source.

**Decision.** Keep independent founders for the pre-release; implement a source deme in v0.2.0.

**Rationale.** A source deme is structural: it adds source diversity as a calibration axis and
changes what an import *is*. It is not required to fix D17's defects, and the pre-release scope is
already treatment + recrudescence.

**Alternatives rejected.** *Implement the source deme now* — scope, and it would land
simultaneously with the treatment rework, making D17's regression ambiguous.

**Consequences.** Imports are maximally distinguishable by construction, so any result reporting
how identifiable imported cases are is conditional on that assumption: the E-id import arm's
AUC 0.80, with import founders near-disjoint (overlap 0.03), would fall under a realistic source.
This must be stated in the README limitations and in the software note, not left implicit. Under a
source deme the inoculum becomes a real genome, so D17's primary→batch edge is superseded by the
sibling topology of an ordinary bite, as a by-product rather than new work.

**How to verify.** README limitations names it; the software note states the assumption alongside
any import-identifiability result.

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
