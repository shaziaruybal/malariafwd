# Scoping: restructure the infectious bite around liver-stage biology

**Status: proposed, not implemented.** Target `malariafwd` 0.2.0 (breaking).
Written 2026-10-01. Becomes ADR D17 if agreed.

---

## 1. The problem

A bite currently makes **two independent draws**:

```
k          = size(broods)                                 // from LAMBDA_SPOROZOITES
bite_batch = drawHypnozoiteBatchSize(effectiveNuBatch())  // Geom(nu), drawn separately
```

Nothing links them, and `LAMBDA_SPOROZOITES` defaults to `0`, which `drawZTCount` maps to a
literal `1`. So at the shipped default **exactly one genotype is inoculated per bite,
deterministically, regardless of what the mosquito carries.**

That is not a configuration choice a user can be expected to notice, and it is not a mild one.
Measured at moderate transmission (`analysis 05_cotransmission_sweep.R`, 3 seeds):

| setting | mosquitoes carrying >1 lineage | mean lineages in mosquito | genotypes **inoculated** per bite |
|---|---|---|---|
| λ_ooc 0, λ_spz 0 (**default**) | 100% | 4.45 | **1.00** |
| λ_ooc 2, λ_spz 2 | 97% | 3.58 | 1.08 |
| λ_ooc 4, λ_spz 4 | 100% | 9.78 | 1.44 |

The mosquitoes are richly polyclonal. Uptake and meiosis generate that diversity correctly.
**The inoculation step discards it.** Raising λ barely helps, because the new parameter samples
*with replacement and then collapses duplicates*, whereas the retired `N_OOCYSTS_PER_BITE` was a
*cap* that kept every distinct ready lineage. Those are not the same scale — 4 draws with
replacement from 4 equally abundant lineages yield 2.74 distinct in expectation, fewer once
`spz_weight` skews the weights. This is also a silent migration hazard: the tripwire catches
`N_OOCYSTS_PER_BITE != 1`, but a user translating their old `= 4` into `LAMBDA_SPOROZOITES = 4`
gets a materially weaker model with no warning.

**Consequence for the science.** Heterologous relapse requires a host's reservoir to contain
genotypes absent from the acute blood-stage infection. With one genotype per bite that can only
arise from separate bites. We therefore structurally suppress the single most characteristic
genetic feature of vivax relapse — in a model built to study relapse identifiability.

---

## 2. Evidence

### 2.1 Fate determination — the crux parameter

[Maher et al. 2022, liver-stage fate determination from patient isolates](https://pmc.ncbi.nlm.nih.gov/articles/PMC9539820/)
measured the **proportion of liver forms that become hypnozoites**: 57.9 ± 0.32%, 61.8 ± 0.32%,
62.89 ± 0.31% across three isolates; strain-dependent and varying with geographic origin, with
some of 51 isolates "remarkably high". Total parasites scale with inoculum (211.94 ± 16.66 per
well at 30×10³ sporozoites ≈ 0.7% establishment in vitro).

→ **p_hyp ≈ 0.6**, strain-specific. This is also the natural home for the tropical/temperate
relapse phenotype, which we currently encode via `ALPHA_ACTIVATE` / `K_LATENCY` only.

### 2.2 Liver-stage infections per bite — three independent routes converge

| route | calculation | N_liver |
|---|---|---|
| White's ν ÷ measured fating | 5 hypnozoites ÷ 0.6 | **8.3** |
| Inoculum × establishment efficiency | 862 spz ([Andolina 2024](https://www.nature.com/articles/s41467-024-44962-4)) × ~0.7% | **~6** |
| Saimiri (needle, 13,000 spz) | "<10 hepatic schizonts", 2–4 genotypes ([Nat Commun 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC10896380/)) | lower efficiency, but ≥10 liver forms |

Routes 1 and 2 are independent and land on **6–8**. Route 3 is a needle inoculum into a monkey
and gives lower per-sporozoite efficiency, as expected.

**The consistency check that makes this worth doing:** if N_liver ≈ 8 and p_hyp ≈ 0.6, then
ν emerges as ≈ 4.8 — essentially White's 5, recovered from liver biology rather than imported
from a genotype-free model. And ~3 of the 8 become schizonts, i.e. **2–3 blood-stage genotype
draws per bite**, which is where co-transmission comes from.

Caveat to carry: sporozoite inoculum is contested — median 15 (range 0–978) in older work vs
geometric mean 862 (IQR 354–2,666) in Andolina. In-vitro hepatocyte efficiency is not in-vivo
efficiency. N_liver should be a *parameter with an evidence-based prior*, not a fixed constant.

### 2.3 Oocysts per infected mosquito

Ethiopian membrane feeding with wild *An. arabiensis*
([Tadesse et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10457422/)): median **7.2
oocysts/midgut** in wild, 17.8 in colony mosquitoes; 63–74% of gametocyte-positive samples
infective. Our effective default is **1**. Symptomatic-patient DMFA overestimates relative to
natural infections, so treat 7.2 as an upper anchor — but 1 is clearly too low.

### 2.4 Validation targets (not parameters — these are what we should reproduce)

**[Popovici et al. 2025, Cambodia](https://pmc.ncbi.nlm.nih.gov/articles/PMC13060509/)** —
59 patients relocated to a malaria-free area, so recurrences are relapses by construction:

- 81% (48/59) had ≥1 recurrence; mean 2.25 per patient; median 20 d to first (IQR 14–34)
- 45% (21/47) polyclonal at enrollment; 42% (22/52) polyclonal relapses
- **Polyclonal enrollment → median 3 relapses; monoclonal → median 1 (p < 0.0001)**
- Monoclonal proportion falls 100% (no relapses) → 0% (4 relapses)
- 3 relapses sibling/parent-offspring related, **explicitly attributed to the same mosquito bite**
- Relapses are infectious: 47–67% transmission vs 80% at enrollment

**[Popovici et al. 2018, mBio](https://pmc.ncbi.nlm.nih.gov/articles/PMC5784252/)** — 20 relocated
patients, 60% relapsed in 60 d; 10/20 polyclonal at enrollment; ≥4/8 sequenced relapses polyclonal;
relapse clones identical, meiotic-sibling, or unrelated to enrollment.

**[Imwong et al. 2007](https://pubmed.ncbi.nlm.nih.gov/17330790/)** — novel genotypes in 72% of
recurrences (78% Thailand / 75% Myanmar / 63% India); in **61% of Thai and Burmese patients the
recurrence contained none of the acute genotypes**; mean MOI 2.67 / 2.63 / 1.84; "a high rate of
mixed-genotype infections even in settings where transmission is low".

> **Caveat — do not use the 61% as a clean relapse target.** Imwong's patients were followed in
> endemic settings, so recurrences conflate relapse with **reinfection**, and a reinfection is
> heterologous almost by definition. The figure is an *upper* bound on heterologous relapse.
> **Popovici's relocated cohorts are the clean target** — patients moved to a malaria-free area,
> so every recurrence is a relapse by construction. Validate on Popovici; treat Imwong as
> supporting context for mixed-genotype prevalence, not as the number to hit.

**[Nat Commun 2024, Saimiri](https://pmc.ncbi.nlm.nih.gov/articles/PMC10896380/)** — sequential
inoculation ~14 d apart fails to establish polyclonality; simultaneous succeeds. One
salivary-gland genotype never reached blood stage (a pre-erythrocytic bottleneck).

Our current output at moderate transmission: COI 1.15–1.27, 0–17% multi-genotype bites.

### 2.5 Where this sits in the modelling literature

Two literatures have developed in parallel and have not been joined. **The fate probability is
the hinge between them, and it is the one piece neither has.**

#### Lineage A — Pv hypnozoite models: rich reservoir mathematics, no genotypes

| model | inoculation step | genotypes? |
|---|---|---|
| [White 2014](https://elifesciences.org/articles/04692) | primary infection **+** Geom(ν) batch per bite | no |
| White 2016 / 2018 (agent-based) | same; adds exposure heterogeneity | no |
| [Mehra 2022 superinfection](https://pmc.ncbi.nlm.nih.gov/articles/PMC10692056/) | same; adds MOB / superinfection | no |
| [Mehra 2024 hybrid](https://pmc.ncbi.nlm.nih.gov/articles/PMC11108905/) | same; adds immunity | no |
| Mehra & Taylor 2025 branching process | same; i.i.d. batches across bites | no |
| Anwar 2022 multiscale | White-derived | no |
| [Champagne 2022 (Swiss TPH, `VivaxModelR`)](https://pmc.ncbi.nlm.nih.gov/articles/PMC8786669/) | extends White 2016; case management (liver/blood stage) + importation; calibrates transmission from observed incidence, derives Rc (Panama) | no |

**The Mehra/Flegg papers state our mechanism explicitly and then do not implement it.** Mehra 2022:
"a *P. vivax* parasite (sporozoite) injected into a human host has **two possible fates**: it either
gives rise to a primary (blood-stage) infection within approximately 9 days of the bite itself
(Mikolajczak 2015), or develops into a hypnozoite (Mueller et al. 2009)." Mehra 2024 reasons the
same way for immunity — "due to the potential for **each sporozoite** to develop into a hypnozoite,
pre-erythrocytic immune protection has been hypothesised to substantially mitigate the relapse
burden". But the formal model in both is "(i) an additional primary blood-stage infection; and
(ii) a geometrically-distributed hypnozoite inoculum of mean size ν" — primary **and** batch, in
parallel, never a competing fate.

That is not an oversight on their part: without genotypes there is nothing a fate split would buy
you, since 1 primary + ν hypnozoites and N_liver × fating are indistinguishable at the level of
counts. It only becomes observable once the forms carry identity. Champagne's work sits further
downstream again — it consumes White's structure to calibrate transmission from surveillance data,
so the inoculation step is inherited rather than revisited.

The [Anwar 2024 scoping review](https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1011931)
confirms the pattern across all 48 Pv models (38 compartmental, 9 agent-based): **none represent
genotypes**, and White 2018 and descendants "do not account for the variability in hypnozoites
across mosquito bites".

#### Lineage B — genotype-resolved vector genetics: no hypnozoites

| model | inoculation step | hypnozoites? |
|---|---|---|
| [Wong 2018](https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1005923) | oocyst count **and infected hepatocyte count** drawn; sporozoites sampled abundance-weighted into hepatocytes; hepatocyte count = 1 excluded as degenerate | no (Pf) |
| [IDM within-vector 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5440164/) | 3000 spz/oocyst, 20% reach glands; **10 sporozoites per bite drawn at random from the salivary pool, duplicate genotypes RETAINED, not collapsed** | no (Pf) |
| EMOD (genetics module) | genotype-resolved Pf: explicit chromosome structure, 25 loci, trackable recombination events | **vivax module not verified** |
| SIMPLEGEN / PlasmoSim | not verified — SIMPLEGEN's genetic-model docs are still TODO/404 | — |
| recombuddy | no mosquito stage at all; downstream painting, not a comparator | — |

The IDM within-vector model is the closest existing precedent for the *sampling* half of the
proposal: draw a fixed number of sporozoites from the salivary pool at random and **keep
duplicates**. That is exactly the step we currently collapse.

#### What the proposal actually is

Not an invention — a join. Take the fate dichotomy that Mehra/Flegg state but cannot use, give it
the measured value from §2.1, and apply it to Wong's hepatocyte primitive with IDM's
don't-collapse sampling. Hypnozoites enter a genotype-resolved vector-genetics model for the first
time, and in exchange ν stops being an input and becomes a prediction we can check against White.

Per-bite (rather than per-genotype) batches are standard across lineage A — White, malariasimulation,
OpenMalaria "brood", Mehra/Anwar — but every one of those lacks genotypes, so that convention is an
absence of capability, not evidence for the structure.

---

## 3. Proposed design — two changes, A and B

**ν and `NU_HYPNOZOITE_BATCH` are kept exactly as they are.** The per-bite Geom(ν) batch stays.
This is deliberate: ν is the one hypnozoite parameter with a literature value, a calibration
history (White, Mehra) and our own fit. Making it emergent would trade one well-identified
parameter for two poorly-identified ones and break comparability with the models we validate
against. See §3.3 for what we are *not* doing and why.

### A. A bite delivers several genomes — raise the default, keep the collapsing

Currently `sampleTransmissibleMosqBroods` draws `n ~ ZT(LAMBDA_SPOROZOITES)` products weighted by
`spz_weight`, then discards duplicates (`draws[first]`), and `LAMBDA_SPOROZOITES = 0` pins `n` to 1.

Change: **give `LAMBDA_SPOROZOITES` a non-zero default. That is all.**

**Correction to an earlier draft of this doc, which proposed also removing the duplicate collapse.
That was wrong and is withdrawn.** Collapsing does not reduce *diversity* — the number of distinct
lineages among `n` draws is identical whether or not duplicates are then dropped. Collapsing only
discards *multiplicity*. So it is not a cause of `k = 1`; `n = 1` is the sole cause, and raising the
default fixes it completely.

Removing the collapse would also be actively harmful, though **not** via Simpson's index —
`computeOutcrossingProb` is now **logged only** ("the REALISED outcome now emerges from the
independent gamete draws below", [`:2007`](inst/slim/malariafwd_vivax.slim#L2007)). The binding
mechanism is the gamete draw itself:

```
gam = sample(active, 2, replace=T);   // UNIFORM over Individuals in the gametocyte pool
```

Uniform over **Individuals**, so a lineage present as two Individuals has twice the chance of
being drawn as a gamete. Duplicates would therefore raise that lineage's selfing probability,
cutting recombination and new-lineage generation — worsening the monoclonality trap this change
exists to fix. The comment at [`:1865`](inst/slim/malariafwd_vivax.slim#L1865) was right:
"duplicates would inflate that lineage's weight in the recipient's gametocyte pool".

(Noted in passing: that uniform draw uses infection counts as a proxy for gametocyte density.
Real outcrossing depends on each clone's gametocyte density, which replication and immunity
govern, not on how many inoculation events delivered it. Out of scope, worth recording.)

### B. Each hypnozoite is genotyped independently from the same pool

Currently the Geom(ν) batch is split uniformly across the `k` inoculated genotypes, and each of
those `k` also starts a primary with probability `P_PRIMARY` (**default 1.0**). So at the default
every genotype carrying hypnozoites also appears in blood: **reservoir genotypes are a subset of
blood-stage genotypes** — and at `k = 1`, a subset of one.

Precisely: the nesting comes from `P_PRIMARY = 1.0` combined with the batch being divided among the
inoculated set. `P_PRIMARY < 1` already allows a dormant-only genotype, but still only one drawn
from the same `k`.

Change: draw each hypnozoite's genotype **independently from the same weighted salivary pool the
primaries came from**, not from the `k` that happened to be sampled for primaries.

```
primaries   = sampleTransmissibleMosqBroods(m)                                // A: default n > 1
bite_batch  = drawHypnozoiteBatchSize(effectiveNuBatch())                     // unchanged, Geom(nu)
hypnozoites = sample(ready_pool, bite_batch, replace=T, weights=spz_weight)   // B: own draws
```

Equivalently: `n + ν` liver infections are drawn from one pool, the first `n` develop and the rest
go dormant. That is Wong's structure with the fating step left deterministic — which is why B gets
most of the benefit of §3.3 without touching ν.

### 3.1 What this buys

- **Co-transmission emerges** at any transmission level from the mosquito's composition — no switch.
- **Heterologous relapse emerges from a single bite.** This is B's work, and it is the one that
  matters most: Imwong's 61% of Thai/Burmese recurrences containing *none* of the acute genotypes
  is hard to produce when the reservoir is a subset of the blood-stage infection.
- **The uniform-split question dissolves** — abundance weighting now applies consistently at both
  steps, so there is no batch to divide.
- ν, the calibrated quantity, is untouched.

### 3.2 Defaults — and why the §2.2 bracket does NOT transfer

**Correction.** An earlier draft pointed `LAMBDA_SPOROZOITES` at the 6–8 "liver forms per bite"
bracket from §2.2. **That was an error introduced when the fate split was deferred.** The 6–8
figure is *total* liver-stage infections (schizonts **plus** hypnozoites). Under A-only those are
two separate draws: liver forms = `k` primaries **+** `ν` hypnozoites. With ν ≈ 3–5 already
accounted for, a total of 6–8 leaves only **k ≈ 1–3 primary-forming genotypes**, not 6–8.

Direct anchor instead: the Saimiri animals established **2 and 4 distinct genotypes** in blood
from "<10 hepatic schizonts". That brackets `k` at roughly **2–4** — and, importantly, from a
needle inoculum of 13,000 sporozoites dissected from **20 mosquitoes**, i.e. ~20 bites' worth. For
a single natural bite the honest read is **lower**, so k ≈ 1–3.

So: `LAMBDA_SPOROZOITES` chosen for a realised mean of roughly **2–3 draws** (λ ≈ 2–2.5; remember
λ is **pre-truncation**, so λ = 2 realises 2.31). Use `uptake_moments()` to set it rather than
reading λ as the mean.

**`LAMBDA_OOCYSTS`: no value proposed, and it may not need one.** The Ethiopian median of 7.2
oocysts/midgut is from **symptomatic-patient membrane feeding**, which overestimates natural
infections, and it is a realised count whereas λ is pre-truncation. Our own sweep also shows it
does little for co-transmission on its own (λ_spz 2 alone: 6.7% multi-genotype bites; adding
λ_ooc 2: 8.0%) — because a *single* outcrossed oocyst already yields four distinct sibling
lineages. Raising it mainly adds half-sib and unrelated combinations plus computation. **Treat as
a separate decision, not part of A.**

### 3.3 Considered and deferred: fate determination

Replacing both draws with one primitive — `N_liver` liver forms per bite, each independently
schizont or hypnozoite with probability `P_HYPNOZOITE` ≈ 0.6 (§2.1) — is more elegant and would
make ν emergent (≈ N_liver × p_hyp ≈ 4.8, i.e. White's 5 recovered from liver biology).

**Deferred, for a calibration reason that is decisive.** White and Mehra model primary **and**
batch in parallel rather than as competing fates, and without genotypes the two are
*observationally identical* — same counts, nothing to gain. For us they would be distinguishable,
but ν is well identified and `N_liver`/`P_HYPNOZOITE` are not, so the swap loses more than it
gains. A and B fix both scientific problems without it.

Also deferred: making the inoculum scale with the biting mosquito's own oocyst/sporozoite burden
(Andolina). Biologically right, needs the field data properly, should not ride along.

### 3.4 B is conditional on A failing — decide it empirically

B is **not** in the same category as the fate split. The fate split changes hypnozoite *counts*,
which we can never observe. B changes the reservoir's *genotype composition*, which we observe
indirectly but genuinely — **through relapse genotypes**. Imwong's 61% of Thai/Burmese recurrences
containing none of the acute genotypes, and Popovici's monoclonal initial infections relapsing as
a different single genotype, are exactly the signal B would move.

But A may be enough on its own. With `k > 1` the reservoir already holds several genotypes, and
because blood-stage infections clear in ~60 days while hypnozoites persist for months, a relapse
from an earlier bite can easily be heterologous to whatever acute infection is present at
enrollment. **A alone therefore produces some heterologous relapse; the question is how much.**

So: **implement A, measure the heterologous fraction against Imwong's 61% and the Popovici
patterns, and only add B if A cannot reach them.** Parsimony first, and the decision is empirical
rather than aesthetic. If A suffices, B never ships and the model stays closer to White/Mehra.

### 3.5 Noted for a later phase: per-oocyst sporozoite productivity

`spz_weight` deliberately stores one Individual per distinct lineage and records how many of the
oocyst's four products it stands for — a storage optimisation that keeps the object count from
scaling, corrected for at sampling. **Keep it.**

It does carry one assumption: every oocyst contributes equally (total weight 4). The
[IDM within-vector model](https://pmc.ncbi.nlm.nih.gov/articles/PMC5440164/) instead gives each
oocyst **~3,000 sporozoites drawn from a Poisson**, of which **20% reach the salivary glands**
(binomial), then draws 10 sporozoites per bite from the resulting pool. That makes salivary-gland
composition depend on per-oocyst productivity and timing, so an early or productive oocyst
dominates the inoculum.

Not needed for A/B, and it could be added **without changing the storage scheme** — productivity
would just become a per-oocyst multiplier on `spz_weight` instead of a fixed 4. Worth a
sensitivity test in the next phase: does unequal oocyst productivity change realised
co-transmission and within-host relatedness?

---

## 4. Migration and risk

| item | impact |
|---|---|
| `NU_HYPNOZOITE_BATCH` | **unchanged.** Still set, still calibrated, still comparable to White/Mehra. |
| `LAMBDA_SPOROZOITES` | meaning unchanged, collapsing unchanged; **only the default moves**. |
| Gametocyte pool weighting | **no longer at risk** — A keeps the collapse, so `computeOutcrossingProb` is untouched. This was the reason for withdrawing the no-collapse variant (§3 A). |
| `P_PRIMARY` | unchanged (no fate split in A+B). |
| `MAX_MOI` | gate currently skips a whole inoculation; semantics need rechecking when a bite delivers several genotypes. |
| Calibration | **COI, prevalence and relapse rate all move.** Must land *before* the full-scale biting recalibration, not after. |
| Performance | more Individuals per bite; tree size and IBD cost grow. Needs a timing check at full scale. |
| Version | behaviour-changing → **0.2.0**, NEWS entry, regression tests for the new defaults. |

**Identification note.** Choose `LAMBDA_LIVER` and `P_HYPNOZOITE` from liver biology (§2.1–2.2),
then calibrate `BITING_RATE` to prevalence/EIR. Do **not** fit the liver parameters to COI: COI is
jointly determined by them and by biting rate, so fitting both to COI leaves them confounded — the
same coupling already recorded in the relapse-identifiability work.

**Effect on existing results.** None of the completed analyses are invalidated by the *current*
per-bite batch rule (they ran at ~1 genotype/bite, where every version coincides). They would,
however, all need re-running under this restructure, which is the real cost of the change.

---

## 4b. Decisions NOT grounded in data or discussion

Flagged explicitly so they are reviewed rather than inherited.

| decision | grounding | status |
|---|---|---|
| `LAMBDA_SPOROZOITES` realised mean **2–3** | Saimiri 2–4 genotypes, but from ~20 bites' worth of sporozoites, so a loose upper bound. No measurement of genotypes per *natural* bite exists. | **Weakly grounded — sweep it, don't trust the point value.** |
| `LAMBDA_OOCYSTS` value | Ethiopian median 7.2 is symptomatic-patient DMFA (overestimates) and is a *realised* count vs our pre-truncation λ. | **No value proposed.** Separate decision. |
| B is conditional on A | My call. Parsimony + the §3.4 argument that A may suffice. Not from literature. | **Judgement, not evidence.** Reasonable but not forced. |
| "~60 days of blood-stage infection makes A sufficient" (§3.4) | Uses our own `MEAN_BLOOD_STAGE_DURATION = 60` default, which is itself a modelling choice, not a measurement. | **Circular if used as proof.** It motivates the test; it does not settle it. |
| p_hyp ≈ 0.6 (§2.1) | Read from a fetched summary of Maher 2022, not the primary PDF. Three isolates quoted; the 51-isolate spread was never quantified. | **Deferred anyway**, but re-read the paper before it is ever used. |
| ν unchanged | Our own fit returned NU = 3 inside all published bands; White/Mehra comparability. | Grounded. |
| Keep the duplicate collapse | Verified from code: collapsing loses multiplicity not diversity, and removing it raises selfing via the uniform gamete draw. | Grounded. |
| Land before the full-scale biting recalibration | COI/biting confounding, already recorded in the relapse-identifiability work. | Grounded. |
| 0.2.0 + NEWS + regression tests | SemVer 0.x convention. | Convention. |

## 5. Suggested order

1. Agree the `LAMBDA_SPOROZOITES` / `LAMBDA_OOCYSTS` defaults from §2 — **this doc, now**.
2. **Implement A only** (raise the default; collapsing and ν unchanged).
3. Regression: ν per bite must be **unchanged** (≈ `NU_HYPNOZOITE_BATCH`) — A redistributes
   hypnozoite *genotypes*, it must not change hypnozoite *counts*.
4. Validate against §2.4 — **primary target: Popovici's relocated cohorts** (polyclonal-vs-
   monoclonal relapse counts, proportion of relapses unrelated to enrollment, relapse
   infectivity). Imwong's 61% is an upper bound only; see the caveat in §2.4.
5. **Decision point (§3.4):** if A reproduces the Popovici patterns, stop — B never ships.
   If it falls short, add B and re-validate.
6. Recalibrate `BITING_RATE` to prevalence/EIR at full scale.
7. Re-run MMCA and China on the updated model.

Deferred to a later phase: per-oocyst sporozoite productivity (§3.5), burden-scaled inoculum
(§3.3), fate determination (§3.3).
