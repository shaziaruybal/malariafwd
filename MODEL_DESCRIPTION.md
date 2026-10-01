# malariafwd v0.1.0 — model description

A forward, individual-based, genotype-resolved simulator of *Plasmodium vivax* transmission with
an explicit hypnozoite reservoir, implemented in SLiM 5 (Eidos) with tskit tree-sequence recording
and an R interface.

Engine: `inst/slim/malariafwd_vivax.slim`. Everything below was read from that file; where a
default is quoted it is the shipped value.

![Model structure](man/figures/model-structure.svg)

---

## 1. What the model is for

**One SLiM individual = one parasite genome.** Humans and mosquitoes are SLiM subpopulations;
the individuals inside them are parasite lineages, not hosts. The whole genealogy is recorded in a
tree sequence, so **identity-by-descent is ground truth read from the ARG, never inferred from
markers**.

That combination — genotype resolution *and* a hypnozoite reservoir — is what the model exists
for. It is built to answer questions of the form *"what can genetic data recover about vivax
recurrence, and where does that break down?"*: relapse vs reinfection identifiability, the value of
sequencing effort, reservoir dynamics under control and resurgence.

**It is not** a tool for predicting case counts, fitting incidence time series, or costing
interventions. Compartmental vivax models (White, Mehra/Flegg, Champagne) do those far better and
faster. See §12.

---

## 2. Scale and the tick cycle

Time advances in **daily ticks**. Each tick runs in this order:

| phase | what happens |
|---|---|
| `reproduction()` | hypnozoite activation (relapses) decided and seeded; biting; human→mosquito uptake and meiosis; mosquito→human inoculation; importation |
| `early()` | hypnozoite cleanup (retired batches killed); blood-stage clearance; mosquito death and replacement; interventions; `state.tsv` logging |

Ordering matters for one subtlety: a batch exhausted or age-capped during `reproduction()` is
**retired** by setting `batch_count = 0`, and killed by the cleanup pass at the top of `early()` in
the **same** tick, so there is no activation lag. (This replaced a queue that an Eidos scoping slip
had made function-local — the fix is WS-G, §13.)

Population sizes are **fixed**: `N_HOSTS` humans, `N_MOSQ` mosquitoes. A mosquito that dies is
replaced by a susceptible one, so vector density is constant and there is no larval or
density-dependent vector ecology.

---

## 3. Human hosts

Hosts carry a set of parasite lineages and are summarised by four states (Panel A):

| state | meaning |
|---|---|
| S₀ | uninfected, no reservoir |
| I₀ | blood-stage infection, no reservoir |
| S_L | hypnozoites only — **invisible to blood surveys** |
| I_L | blood stage **and** hypnozoites |

**No intrinsic incubation.** Blood stage begins at inoculation, so hosts are *not* SEI. Liver-stage
development time (~9 days in reality) is not represented.

**No immunity, by design.** There is no acquired blood-stage, anti-parasite, pre-erythrocytic or
transmission-blocking immunity, and no age structure. This is a deliberate scope decision: immunity
would confound the genetic signals the model exists to measure, and vivax immunity parameters are
poorly identified. It is the single largest simplification in the model (§11, §14).

**Blood-stage duration** is exponential with mean `MEAN_BLOOD_STAGE_DURATION` = 60 days, with a
hard ceiling `MAX_BROOD_AGE` = 180 days. Clearance is per-lineage, so a polyclonal host loses
lineages independently.

**Superinfection** is allowed and gated by `MAX_MOI` = 20. The gate is checked once, in the
infected-mosquito-bites-infected-human branch: a host already carrying `MAX_MOI` active lineages
has the **whole** inoculation skipped — primary *and* hypnozoite batch, so the reservoir also stops
being topped up. **Relapses bypass the gate entirely**, so realised MOI can exceed it. At 20 it is
a runaway guard that essentially never binds (COI is typically 1–3); set it low and you are
measuring the cap, not the biology.

---

## 4. Mosquitoes

Standard S_v → E_v → I_v, with:

- **Biting:** `BITING_RATE` = 0.1 is the per-*mosquito* daily biting probability. The
  Ross–Macdonald per-human rate is `ma = BITING_RATE × N_MOSQ / N_HOSTS`, **not** this value
  directly — a frequent source of confusion when comparing to the literature.
- **Transmission probabilities:** `P_TRANS_VH` = 0.5 (`b`, bite infects human),
  `P_TRANS_HV` = 0.5 (`c`, mosquito infected by biting an infected human).
- **Lifespan:** exponential, mean `MOSQ_LIFESPAN` = 10 days.
- **EIP:** **fixed** at `EIP_DAYS` = 10 days — not temperature-dependent, not distributed.
- **Carrying capacity:** `MAX_BROODS_PER_MOSQ` = 20 lineages; uptake is skipped above it. No
  fitness cost of burden (mortality and transmissibility do not decline with parasite load).

Mosquito death removes all of its broods, so a mosquito is a hard transmission bottleneck.

---

## 5. The hypnozoite reservoir

Follows the **White et al. 2014** batch model.

- Each transmitted genotype may cause an **active primary infection** with probability
  `P_PRIMARY` = 1.0 (White's base model treats this as certain given successful inoculation), and
  **independently** establishes a batch of hypnozoites.
- **Batch size** is geometric with mean `NU_HYPNOZOITE_BATCH` (ν) = 5.
  **Since WS-F the batch is drawn once per *bite*** and split uniformly across the co-inoculated
  genotypes, matching White's per-bite parameterisation (§12). Both, either, or neither of primary
  and batch can occur for a given genotype.
- **Activation:** each hypnozoite activates at `ALPHA_ACTIVATE` = 1/170 per day (Mehra; White 2014
  would give 1/344). An activation produces **one blood-stage episode** representing however many
  clonal hypnozoites activated that tick.
- **Death:** `MU_HYPNOZOITE` = 1/442 per day (White 2014), competing with activation.
- **Latency:** `K_LATENCY` = 0 (tropical) or 2 (temperate), with `DELTA_LATENCY` = 1/100 governing
  progression when enabled.
- **Age ceiling:** `MAX_HYPNOZOITE_AGE` = 1500 days — a safety ceiling only, ~3.4× White's 442-day
  mean lifespan, so death is governed by the exponential and the long late-relapse tail survives.

**A relapse is a clone of its originating hypnozoite** (`addCloned`), so relapse-to-origin IBD is
exactly 1. Relative to a host's *prior observed* infection it can be clonal, sibling, half-sibling
or unrelated — which is the whole point of the identifiability work.

> **Structural consequence worth stating explicitly.** Because the batch is divided among the
> genotypes that establish primaries, and `P_PRIMARY` = 1, **reservoir genotypes are a subset of
> blood-stage genotypes**; at the default (one genotype per bite) both are the same single
> genotype. Heterologous relapse can therefore only arise across separate bites. See §13.

---

## 6. The mosquito stage: uptake, meiosis, sporozoites

This is where the genetics are generated, and where v0.1.0 (WS-F) changed most.

**Uptake.** A successful blood meal forms `n ~ zero-truncated Poisson(LAMBDA_OOCYSTS)` oocysts
(negative binomial if `K_OOCYSTS > 0`, var = λ + λ²/K). Each oocyst is an independent meiosis with
its own `oocyst_id`; all oocysts from one meal share a `brood_id`.

> **λ is the PRE-truncation rate, not the mean count.** Conditioning away the zeros pulls the
> realised mean up: λ = 1 → 1.58, λ = 2 → 2.31; with K = 1, λ = 2 → 3.0. Use `uptake_moments()` /
> `uptake_calibrate()` to convert between a literature target and λ. Zero-*truncation* (not
> inflation) is correct because `P_TRANS_HV` already decides whether the event succeeds — allowing
> zeros would count failure twice.

**Gamete draw.** Each oocyst draws **two gametes independently and uniformly** from the host's
gametocyte pool. Selfing vs outcrossing **emerges** from whether the two happen to share a lineage;
it is not rolled as a coin flip. This reproduces P(outcross) = 1 − Σpᵢ² as a marginal *and* gives
the correct pairing P(i,j) ∝ pᵢpⱼ (Wong et al. 2018), which the previous "roll outcross, then pick
two distinct lineages" did not.

> `computeOutcrossingProb()` (Simpson) still exists but is **logged only**, as the `p_outcross`
> column. It no longer drives anything.

**Meiosis.** An outcrossed oocyst yields a full **tetrad**: two reciprocal recombinants (R12, R21,
sharing breakpoints and therefore IBD 0 to each other) plus two parental clones. A selfed oocyst
yields four identical products, stored as **one** object with `spz_weight = 4` — a storage
optimisation that keeps object count from scaling with clonal output, corrected for at sampling.

**Gametocyte gating.** A blood-stage infection is transmissible only from
`gametocyte_ready_tick = start + GAMETOCYTE_DELAY` and only if it carries gametocytes (per-infection
Bernoulli, `P_GAMETOCYTE`). Defaults are the identity (0 and 1.0). Vivax gametocytes appear early —
at or near symptom onset, unlike ~10–12 d in falciparum — so a realistic preset is a short delay
plus high carriage (≈3 and ≈0.9).

**Inoculation.** A bite draws `n ~ ZT(LAMBDA_SPOROZOITES)` products from everything past the EIP,
**weighted by `spz_weight`**, then collapses repeat draws to one object per distinct lineage.

> **`LAMBDA_SPOROZOITES` defaults to 0, which `drawZTCount` maps to exactly 1.** So at the shipped
> default **every bite delivers exactly one genotype, regardless of what the mosquito carries** —
> co-transmission is structurally impossible. Measured: mosquitoes carry a mean of 4.5 distinct
> lineages while 1.00 is inoculated. This is the subject of a pending change (§13).

`N_OOCYSTS_PER_BITE` is **retired**; it capped lineages rather than oocysts and now stops the run
if set to anything but 1.

---

## 7. Genome and genetics

- **PvP01 reference** (Auburn et al. 2016; GCA_900093555), 14 nuclear chromosomes.
- Implemented as **one concatenated linkage group** with recombination rate **0.5 at each of the 13
  junctions** (independent assortment) and **7.4 × 10⁻⁷ /bp within chromosomes** (13.5 kb/cM,
  euchromatic core rate, Jiang 2011). Genetically equivalent to explicit multi-chromosome for all
  IBD metrics, and keeps a single `.trees` file.
- `GENOME_MODE` is a pure **coordinate** switch: `"full"` (24.2 Mb, default) or `"core"`
  (subtelomere-masked, ~21.8 Mb). Both use the same recombination rate.
- **`initializeMutationRate(0.0)` — the genome carries no mutation at all.** It is a neutral
  scaffold; all genetics come from the recorded genealogy, and alleles are painted on downstream
  (PLAF/panel simulation). This is a strength for IBD (ground truth, no ascertainment) and a
  limitation for anything needing sequence divergence or selection.

> **The genealogy is a FOREST, not a tree.** Founder lineages have independent roots and ~93% of
> pairs never coalesce. Consequently **tskit branch-mode divergence is not a relatedness measure
> here** — use the IBD functions, not `diversity()`/`divergence()`.

**Gotcha for anyone writing new tskit code:** SLiM haploid individuals get **two** tskit nodes — the
real haplosome (`is_vacant = [0]`) and a null placeholder (`is_vacant = [1]`, isolated, IBD 0 with
everything). About half of `ts.samples()` are vacant and **must be filtered**. The package's
`ibd_pairs()` / `ibd_index()` already do this.

---

## 8. Importation

Three mechanisms, additive:

| parameter | behaviour |
|---|---|
| `IMPORT_RATE_PER_TICK` | genuine Poisson arrival rate — the realistic "leaky border" |
| `N_IMPORTS_PER_TICK` | deterministic floor, between `IMPORT_START`/`IMPORT_END` |
| `IMPORT_TICK` / `N_IMPORTS` | single pulse event |

Imports use the **same bite outcomes** as local transmission (primary and/or batch) and are
labelled `is_imported`, which is **carried through to relapses** so an imported relapse is a
distinct truth category.

**v0.1.0 limitation:** imported genotypes are drawn as unrelated founders — there is no structured
source population, so imports are maximally distinguishable from local parasites (measured founder
overlap 0.03). A structured source pool is planned (§15).

---

## 9. Interventions and the observation layer

**Interventions** (all off by default):

- **Vector control / seasonality** scale biting. Seasonality is
  `BITING_RATE × (1 + SEASONAL_AMP·cos(2π(t − SEASONAL_PHASE)/365))`, with `SEASONAL_PHASE` the
  day-of-year of the **peak**.
- **Case-based radical cure** (`RADCURE_COVERAGE`, efficacy `1 − P_RAD_SURVIVAL`) is **liver-stage
  only** — primaquine/tafenoquine. It deliberately does **not** clear the current blood-stage
  infection.
- **Case-based blood-stage treatment** (`BLOODCURE_COVERAGE`, efficacy `1 − P_BLOODCURE_SURVIVAL`)
  clears the active infection only.
- **MDA** (`MDA_TICK`) clears both.
- **Resurgence triggers** from `RESURGENCE_TRIGGER_TICK`: multiplicative boosts on biting
  (`BITING_RATE_BOOST`), activation (`ALPHA_ACTIVATE_BOOST`), batch size
  (`NU_HYPNOZOITE_BATCH_BOOST`) and imports (`IMPORT_RATE_BOOST`) — so one scenario can run a quiet
  baseline and then switch on a single mechanism.

**Observation layer.** A true episode becomes a *detected* (reported, genotyped) case with
probability `P_DETECT_PRIMARY` / `P_DETECT_RELAPSE`. Relapses are often sub-patent, so a realistic
setting is much lower for relapses (e.g. 0.15). Detected cases form the observed line list and can
trigger case-based treatment.

> **Known gap:** detection is not origin-dependent beyond the primary/relapse split, so
> *active vs passive* surveillance modes cannot currently be distinguished. This made one arm of
> the survey-design work come out null.

---

## 10. Parameters

λ values are **pre-truncation rates**; `GENOME_LENGTH`, `N_CHROMOSOMES`, `RECOMB_WITHIN` are
derived and written to `params.tsv` for provenance.

### Population and epidemiology

| parameter | default | meaning / source |
|---|---|---|
| `N_HOSTS` | 40 | humans |
| `N_MOSQ` | 100 | mosquitoes (constant) |
| `N_FOUNDERS` | 12 | legacy; decomposed below |
| `N_FOUNDER_HOSTS` | `N_FOUNDERS` | infected hosts at t = 0 (**establishment**) |
| `N_FOUNDER_LINEAGES` | `N_FOUNDER_HOSTS` | distinct founding genomes (**diversity** knob — the clean genealogical axis) |
| `BITING_RATE` | 0.1 | per-**mosquito** daily biting probability |
| `P_TRANS_VH` | 0.5 | `b` |
| `P_TRANS_HV` | 0.5 | `c` |
| `MOSQ_LIFESPAN` | 10 d | exponential |
| `EIP_DAYS` | 10 d | fixed |
| `MAX_MOI` | 20 | superinfection gate (relapses bypass) |
| `MAX_BROODS_PER_MOSQ` | 20 | mosquito lineage cap |
| `MEAN_BLOOD_STAGE_DURATION` | 60 d | |
| `MAX_BROOD_AGE` | 180 d | ceiling |

### Hypnozoite reservoir

| parameter | default | source |
|---|---|---|
| `P_PRIMARY` | 1.0 | White 2014 base model |
| `NU_HYPNOZOITE_BATCH` | 5.0 | White 2014 §2, geometric. **Placeholder** — calibrate per setting |
| `ALPHA_ACTIVATE` | 1/170 /d | Mehra (White 2014 = 1/344) |
| `MU_HYPNOZOITE` | 1/442 /d | White 2014 |
| `K_LATENCY` | 0 | 0 tropical, 2 temperate |
| `DELTA_LATENCY` | 1/100 | used only if `K_LATENCY > 0` |
| `MAX_HYPNOZOITE_AGE` | 1500 d | safety ceiling, ~3.4× mean lifespan |

### Mosquito-stage genetics

| parameter | default | note |
|---|---|---|
| `LAMBDA_OOCYSTS` | 0.0 | → exactly 1 oocyst. Field vivax DMFA median is **7.2/midgut** (wild *An. arabiensis*, Ethiopia) |
| `K_OOCYSTS` | 0.0 | 0 = Poisson; > 0 = NB. Vivax DMFA counts are overdispersed |
| `LAMBDA_SPOROZOITES` | 0.0 | → exactly 1 genotype per bite. **Co-transmission impossible at the default** |
| `N_OOCYSTS_PER_BITE` | 1 | **retired**; tripwire |
| `GAMETOCYTE_DELAY` | 0 | vivax realistic ≈ 2–3 d |
| `P_GAMETOCYTE` | 1.0 | vivax realistic ≈ 0.9 |

### Genome

| parameter | default |
|---|---|
| `GENOME_MODE` | `"full"` (24.2 Mb) |
| `CHROMOSOMES` | `"1-14"` |
| `RECOMB_RATE` | −1 (auto → 7.4e-7/bp) |
| mutation rate | **0** (hard-coded) |

### Interventions, observation, importation

| parameter | default |
|---|---|
| `INTERVENTION_START` / `_END` / `_EFFICACY` | 99999 / 99999 / 0 |
| `SEASONAL_AMP` / `SEASONAL_PHASE` | 0 / 0 |
| `P_DETECT_PRIMARY` / `P_DETECT_RELAPSE` | 1.0 / 1.0 |
| `RADCURE_COVERAGE` / `P_RAD_SURVIVAL` | 0 / 0.05 |
| `BLOODCURE_COVERAGE` / `P_BLOODCURE_SURVIVAL` | 0 / 0 |
| `MDA_TICK` | 99999 |
| `RESURGENCE_TRIGGER_TICK` + 4 boosts | 99999, all 1.0 |
| `IMPORT_*` | off |
| `RETAIN_FULL_HISTORY` | `T` — needed for IBD to already-cleared infections |
| `VERBOSE` | `F` |

---

## 11. Assumptions and simplifications

Stated plainly, because several matter for interpretation.

| # | assumption | consequence |
|---|---|---|
| 1 | **No immunity of any kind** | no age-prevalence gradient, no density/detectability decline with exposure, no transmission-blocking immunity. Deliberate — immunity would confound the genetic signal. **Largest simplification in the model.** |
| 2 | **No liver-stage delay** | blood stage starts at inoculation; recurrence-interval distributions are shifted ~9 d early |
| 3 | **Neutral genome, zero mutation** | no selection, no drug resistance, no sequence divergence; IBD is exact |
| 4 | **Fixed vector population, fixed EIP** | no larval ecology, no temperature dependence, no vector fitness cost of burden |
| 5 | **One genotype per bite at the default** | co-transmission structurally impossible; reservoir ⊆ blood-stage genotypes |
| 6 | **Equal oocyst productivity** | every oocyst contributes weight 4; real oocysts vary (~1,500–5,000 sporozoites) |
| 7 | **Gametocyte pool = infection counts** | outcrossing draws uniformly over *Individuals*; real outcrossing depends on each clone's gametocyte **density**, governed by replication and immunity |
| 8 | **Within-host densities not modelled** | no competition between clones, no density-dependent detectability |
| 9 | **Imports unrelated to local parasites** | imports maximally distinguishable; an upper bound on import-detection performance |
| 10 | **Detection not origin-dependent** | active vs passive surveillance indistinguishable |
| 11 | **Homogeneous mixing** | no spatial structure, no household clustering, no heterogeneous biting |
| 12 | **`NU_HYPNOZOITE_BATCH` is a placeholder** | shipped at 5; calibrate per setting against relapse rate / fraction |

---

## 12. Relationship to existing models

Two literatures, which this model sits between.

### Vivax hypnozoite models — rich reservoir mathematics, no genotypes

| model | inoculation step | genotypes |
|---|---|---|
| White 2014 (eLife) | primary **+** Geom(ν) batch per bite | no |
| White 2016 / 2018 | same, + exposure heterogeneity | no |
| Mehra 2022 (superinfection), 2024 (hybrid), Mehra & Taylor 2025 | same, + MOB / immunity / branching process | no |
| Anwar 2022 multiscale | White-derived | no |
| Champagne 2022 (Swiss TPH, `VivaxModelR`) | extends White 2016; case management + importation; **inverts an endemic steady state** to get λ and Rc from routine incidence | no |

The Anwar 2024 scoping review covers **48** vivax transmission models (38 compartmental, 9
agent-based): **none represent genotypes**, and White 2018 and descendants "do not account for the
variability in hypnozoites across mosquito bites" — which our per-bite geometric draw does.

Notably, Mehra 2022 **states** the biology we do not model — "a *P. vivax* parasite (sporozoite)
injected into a human host has two possible fates: it either gives rise to a primary (blood-stage)
infection … or develops into a hypnozoite" — but implements primary **and** batch in parallel, as
we do. Without genotypes the two are observationally identical, so there was nothing to gain.

### Genotype-resolved vector genetics — no hypnozoites

| model | inoculation step | hypnozoites |
|---|---|---|
| Wong et al. 2018 | oocyst count **and infected hepatocyte count** drawn; sporozoites sampled abundance-weighted into hepatocytes; hepatocyte count = 1 excluded as degenerate | no (Pf) |
| IDM within-vector 2017 | ~3,000 spz/oocyst (Poisson), 20% reach glands, **10 sporozoites per bite drawn from the pool, duplicates retained** | no (Pf) |
| EMOD genetics module | genotype-resolved Pf: chromosome structure, 25 loci, trackable recombination | vivax module unverified |
| SIMPLEGEN / PlasmoSim | genetic-model docs still in development | — |
| recombuddy | no mosquito stage; downstream block-mosaic painting (we use it this way) | — |

**Our position:** we are, as far as this review found, the only model with **both** genotype
resolution and a hypnozoite reservoir. That is the contribution, and it also means several
inoculation-step decisions have no precedent to defer to.

### Conventions we follow

Per-bite (not per-genotype) hypnozoite batches: White 2014, malariasimulation (batch = one
infective bite), OpenMalaria ("brood" = primary + relapses from the same bite), Mehra/Anwar (i.i.d.
batches across bites). Independent, uniform gamete draws giving P(i,j) ∝ pᵢpⱼ: Wong et al. 2018.

---

## 13. Verification, fixed defects, and calibration status

### Defects found and fixed in v0.1.0 (WS-G)

All five were reproduced empirically before fixing. Four were instances of the same Eidos trap: a
plain assignment to a global inside a callback or function creates a **local**.

| id | defect | measured impact |
|---|---|---|
| WS-G1 | hypnozoite batches never retired — exhausted and age-capped batches kept activating | **64% of batches over-drew**; relapses/host-yr fell ~3× after the fix (4.7 → 1.55 at ν = 3) |
| WS-G2 | duplicate event ids (`INOCULATION_ID` never incremented) | all inoculations shared id 1 — grouping by `inoculation_id` is invalid in pre-v0.1.0 output |
| WS-G3–5 | deferred radical-cure flag, MDA empty-host guard, tick-event id | corrected |

**Consequence for older output:** group bites by `(tick, mosq_id, human_id)`, not
`inoculation_id`, in anything produced before v0.1.0.

### Validation suite

`inst/validation/` ships five scripts: `audit_events_and_ibd.R` (62 checks),
`check_parameters.R` (41 checks), `audit_arg_recombination.R`, `audit_arg_junctions.R`,
`explain_ibd.R`. Verified properties include: clonal IBD = 1.000 exactly; relapse-to-origin 100%
clonal; same-oocyst sibling IBD centring on ~1/3 on one chromosome and tightening toward 0.5 under
14-chromosome assortment; reciprocal recombinants IBD 0 to each other.

### What is calibrated and what is not

| quantity | status |
|---|---|
| `NU_HYPNOZOITE_BATCH` | re-fitted on the corrected model; 3 sits inside all published bands (relapses/py, molFOB, median MOI, relapse fraction). Shipped default is **5** and is a placeholder |
| `BITING_RATE` | **not** calibrated at shipped defaults. Transmission targets are scale-dependent (prevalence and COI depend on `N_HOSTS`/`N_MOSQ`), so calibrate at the scale you will publish |
| hypnozoite rates | literature values (White, Mehra), not fitted |
| genome, recombination | reference values, not fitted |

> **Identification warning.** COI is jointly determined by biting rate *and* the uptake/inoculation
> parameters. Fitting both to COI leaves them confounded. Set the parasite-side parameters from
> biology, then calibrate `BITING_RATE` to prevalence/EIR.

### Pending model change

`LAMBDA_SPOROZOITES = 0` makes co-transmission structurally impossible (§6), which suppresses
heterologous relapse — the characteristic genetic signature of vivax recurrence. Scoped in
[`LIVER_STAGE_RESTRUCTURE.md`](LIVER_STAGE_RESTRUCTURE.md), to land before any full-scale
recalibration.

---

## 14. Where sensitivity analysis is needed

Ordered by how much the answer is likely to move.

| # | target | why | status |
|---|---|---|---|
| S1 | **`LAMBDA_SPOROZOITES`** (genotypes per bite) | no measurement of genotypes per *natural* bite exists; the Saimiri anchor (2–4 genotypes) comes from ~20 bites' worth of sporozoites | **priority** — sweep, don't trust a point value |
| S2 | **`NU_HYPNOZOITE_BATCH`** | placeholder at 5; drives relapse rate, reservoir size and COI | sweep per setting |
| S3 | **`LAMBDA_OOCYSTS` / `K_OOCYSTS`** | vivax DMFA counts overdispersed; field median 7.2 vs our effective 1 | sweep; NB option exists |
| S4 | **Per-oocyst sporozoite productivity** | we assume equal (weight 4); IDM draws ~3,000/oocyst Poisson with 20% gland migration, so an early or productive oocyst would dominate | not implemented; test as sensitivity |
| S5 | **`GAMETOCYTE_DELAY` / `P_GAMETOCYTE`** | shipped at identity; vivax realistic ≈ 3 d and ≈ 0.9 | knobs exist for exactly this |
| S6 | **`P_DETECT_RELAPSE`** | relapses are often sub-patent; everything observational depends on it | sweep 0.15–1.0 |
| S7 | **`ALPHA_ACTIVATE` / `K_LATENCY`** | tropical vs temperate strongly changes relapse periodicity and identifiability | already explored; keep as a reported axis |
| S8 | **`MAX_MOI`** | a runaway guard that becomes the measured quantity if set low | confirm it never binds in published runs |
| S9 | **`MEAN_BLOOD_STAGE_DURATION`** | sets episode overlap and therefore polyclonality and outcrossing | sweep |
| S10 | **Scale (`N_HOSTS`, `N_MOSQ`)** | prevalence and COI are scale-dependent; a rate fitted at preview scale does not transfer | always calibrate at publication scale |

---

## 15. Planned extensions

| extension | motivation |
|---|---|
| Co-transmission by default + independent hypnozoite genotyping | §13; makes heterologous relapse emergent |
| Structured import source pool | currently imports are unrelated founders, so import detection is an upper bound |
| Burden-scaled inoculum | sporozoite inoculum tracks salivary-gland load (Andolina 2024); heavily-infected mosquitoes should seed proportionally more |
| Per-oocyst productivity | S4 |
| Origin-dependent detection | to separate active vs passive surveillance |
| Fate determination (`P_HYPNOZOITE`) | measured at 0.58–0.63 and strain-dependent; **deferred** — would make ν emergent but trades one well-identified parameter for two poorly-identified ones |
| Immunity | explicitly out of scope for now; would require a separate identifiability argument |
| Spatial structure / heterogeneous biting | currently homogeneous mixing |
| Vector fitness cost of burden | `MAX_BROODS_PER_MOSQ` is a hard cap with no graded effect |

---

## 16. Outputs, reproducibility and cost

**Files per run:** `events.tsv` (the event log — one row per transmission, activation, clearance,
treatment event, with lineage/pedigree/oocyst ids, `is_dormant`, `is_imported`, `n_hypnozoites`,
`detected`), `state.tsv` (per-tick host and vector state), `params.tsv` (**all resolved parameters
plus `ACTUAL_SEED`** — full provenance), `plasmo.trees` (tskit tree sequence).

**Identifier hierarchy** — worth internalising, these are not interchangeable:

| id | meaning |
|---|---|
| `pedigree_id` | one SLiM individual = one genome; many per lineage |
| `lineage_id` | genetic identity / clonal group; propagates through clonal, parental and selfed descent; recombinants mint a fresh one. **COI counts distinct `lineage_id`** |
| `oocyst_id` | one meiosis / tetrad, ≤ 4 lineages |
| `brood_id` | one blood meal's co-hatching group, ≥ 1 oocyst |
| `inoculation_id` | one bite (**broken before v0.1.0** — see WS-G2) |

**Reproducibility:** `SEED` is recorded as `ACTUAL_SEED` in `params.tsv`. Requires SLiM ≥ 5.0 on
`PATH` and Python with `tskit` for IBD. `RETAIN_FULL_HISTORY = T` is required for IBD to
already-cleared infections — without it, relapse-to-prior comparisons silently lose ancestors.

**Cost:** the dominant expense is **pairwise IBD**, not the simulation. Tree-sequence loading is
second. Simulation itself is cheap (a 150-host, 3-year run is tens of seconds). Budget accordingly
when scaling the number of sampled genomes, and prefer cross-matrices over n² where possible.

---

## 17. Citation and licence

See `DESCRIPTION`, `DECISIONS.md` (ADR log, D1–D16) and `NEWS.md`.

> **Documentation defect to fix:** `NEWS.md` states the genome has "forward mutation". It does not —
> `initializeMutationRate(0.0)`.
