# malariafwd v0.1.0 — validation

How the model is checked: the validation scripts and what each found,
and the numeric audit of the recombination machinery. The [model
description](https://shaziaruybal.github.io/malariafwd/MODEL_DESCRIPTION.html)
says what the model does, and its section 13 gives which parameters are
calibrated to data and which are not.

## Validation suite

`inst/validation/` ships five scripts: assertion checks on a run’s event
log and genealogy, the calibration of `NU_HYPNOZOITE_BATCH`, and two
numeric audits of the recombination machinery.
`inst/validation/RESULTS.md` records the last full run — its date,
commit and versions, and how many checks each script passed. The counts
are there rather than here, so that they come from a run rather than
from editing by hand.

| script | what it checks |
|----|----|
| `audit_events_and_ibd.R` | event-log integrity, transmission links, the ARG and IBD, and genotype painting, on each of three configurations |
| `check_parameters.R` | that the structural invariants hold across 34 configurations, and that each parameter moves the model the way it should |
| `calibrate_nu.R` | that every run completes, that the shipped ν does not die out, and that it sits inside all four targets; it makes the table below |
| `audit_arg_recombination.R` | a numeric audit; see below |
| `audit_arg_junctions.R` | a numeric audit; see below |

### Event-log integrity

`audit_events_and_ibd.R` runs three configurations of 80 hosts and 160
mosquitoes for 400 days, with 8 founder lineages: biting at 0.12 with
one oocyst an uptake and one genotype a bite, every other setting at its
default; the same with several oocysts an uptake, negative-binomial
counts, several products a bite and imports; and biting at 0.16 with
several oocysts, blood-stage treatment whose failures recrudesce and
radical cure. The observed column gives the three in that order.

| property | expectation | observed |
|----|----|----|
| every emitted `event_type` documented in [`event_types()`](https://shaziaruybal.github.io/malariafwd/reference/event_types.md) | 0 undocumented | 0 |
| `event_id` unique | no duplicates | 3,896, 14,707 and 17,014 ids, all distinct |
| creation events mint unique `pedigree_id`s | no reuse | 577, 977 and 1,085 created, all distinct |
| each `inoculation_id` is one bite | 1 tick, 1 host | 197, 269 and 403 bites |
| **no batch draws more hypnozoites than it held** | **0 over-draws** | **150, 242 and 321 batches, 0 over-drawing** |
| every activation references an established batch | 100% | 218, 362 and 265 activations |
| unresolved hypnozoites == final reservoir in `state.tsv` | exact match | 230 vs 230; 313 vs 313; 282 vs 282 |
| each uptake has its own `oocyst_id` | unique | 1,328, 4,896 and 5,945 oocysts |
| oocysts per blood meal | ≥ 1; 1 at λ = 0 | mean 1.00 (max 1), 2.94 (max 22) and 3.05 (max 28) |
| recombinant uptakes record two distinct parents | always | 9.4%, 23.9% and 13.0% recombinant |

The over-draw check guards the retirement of exhausted batches, which
the Eidos scoping trap in the [implementation
notes](https://shaziaruybal.github.io/malariafwd/IMPLEMENTATION_NOTES.html)
would break silently: an unretired batch keeps activating after it has
run out.

### ARG and IBD

Checks that the recorded genealogy keeps every case and gives the
relatedness that the transmission history implies.

| property | expectation | observed |
|----|----|----|
| tree sequence loads | — | 4,406 trees / 4,560 nodes, 42,139 / 18,778 and 27,531 / 19,340, each 24.21 Mb |
| logged cases retained in the tree sequence | 100% | 100.0% of 415, 723 and 1,079 cases |
| ARG recombines ([`arg_summary()`](https://shaziaruybal.github.io/malariafwd/reference/arg_summary.md), over the cases and founders) | \> 1 local tree | 271, 1,581 and 1,022 trees |
| **genealogy is a forest** | founders never coalesce | 8, 51 and 8 roots a tree, **0% coalesced** |
| **relapse-to-origin IBD** | **exactly 1** | **at least 1.000 over 218, 362 and 265 relapses** |
| distinct founder lineages unrelated | IBD ≈ 0 | median 0.000, max 0.000 over 28 pairs in each |
| clonal founder copies | IBD = 1 | 1.000 |

### Parameter-response behaviour

Monotonicity and switch checks — that each parameter moves the model in
the direction it should.

| property | observed |
|----|----|
| higher biting → higher prevalence | 0.92 vs 0.02 |
| higher biting → higher COI | 1.90 vs 1.00 |
| more founder lineages → more circulating lineages | 45 vs 2 |
| `NU_HYPNOZOITE_BATCH = 0` → no relapses at all | 0 relapses, reservoir 0.00 |
| faster activation → higher relapse fraction | 0.70 vs 0.41 |
| faster activation → **smaller** standing reservoir | 0.58 vs 1.03 (3-seed means) |
| faster hypnozoite death → smaller reservoir | 0.43 vs 1.22 (3-seed means) |
| activation boost → higher relapse fraction | 0.67 vs 0.52 |
| batch boost enlarges the reservoir | 3.41 vs 0.73 |
| radical cure → smaller reservoir | 0.02 vs 1.22 (3-seed means) |
| weaker radical cure leaves a bigger reservoir | 0.14 at `P_RAD_SURVIVAL` 0.8 vs 0.02 at 0.05 (3-seed means) |
| blood-stage treatment → lower prevalence | 0.04 vs 0.24, 63 cases treated |
| weaker blood-stage treatment leaves higher prevalence | 0.07 vs 0.04 |
| vector control → lower prevalence | 0.02 vs 0.24 |
| resurgence trigger → higher prevalence | 0.77 vs 0.24 |
| temperate latency delays each batch’s first relapse | median **98 d** vs **23 d** tropical |
| co-transmission (`LAMBDA_SPOROZOITES = 2`) → more outcrossing, 3-seed means at `BITING_RATE` 0.2 | 49.8% vs 42.1% |
| all three importation modes produce imports | 105 / 957 / 54 (a boosted rate, a fixed count, a pulse) |
| `MAX_MOI` gates superinfection | mean COI 1.87 at `MAX_MOI = 2`, 4.28 without |
| detection layer thins the observed list only | baseline 157/157 detected; 68/198 with lower detection |
| `GENOME_MODE = "core"` shortens the genome | 21.8 vs 24.2 Mb |
| `CHROMOSOMES` subset | 3 chromosomes simulated |
| [`run_grid()`](https://shaziaruybal.github.io/malariafwd/reference/run_grid.md) completes and returns one row per (set × seed) | 4 runs, 4 rows, 4 dirs |

## Recombination machinery audit

| property | expectation | observed |
|----|----|----|
| all 13 chromosome junctions used | yes | **yes** |
| junction share of crossovers | 26.6% (6.5 junction vs 17.9 within per product) | consistent |
| within-chromosome breakpoint density uniform across chromosomes | CV ≈ 0 | **CV = 0.008** (7,376–7,628 per Mb) |
| within-chromosome breakpoint position uniform | KS p \> 0.05 | **D = 0.002, p = 0.799** |
| a product’s parent across a junction (independent assortment) | φ ≈ 0 | **φ = 0.01 ± 0.03** at 50.4% the same |
| a product’s parent across 200 kb within a chromosome | φ = e^(−2r × 200 kb) = 0.744 | **φ = 0.712** at 85.6% the same |
| IBD segments two products of one oocyst share (two unrelated founder parents) | mean L/(14 + 2rL) = 0.49 Mb, 16.6 a pair | **0.48 ± 0.02 Mb**, 16.8 a pair (85 pairs from 45 oocysts) |
| IBD segments a product shares with its parent | mean L/(14 + rL) = 0.76 Mb, 16.0 a product | **0.78 Mb**, 15.9 a product (100 products; KS p = 0.07 against the map) |

Segments are counted as the [derived
metrics](https://shaziaruybal.github.io/malariafwd/METRICS.html) define
them, with no length floor. Where two products of one oocyst carry the
same parent they carry its two chromatids, so every chiasma (rate 2r)
ends a shared stretch, and where they differ one begins at rate r; a
product against its parent changes only at a crossover (rate r).

The junction-vs-control contrast, measured on each product against its
parent (φ about 0 against 0.71), is the direct evidence that independent
assortment is working: a product’s parent at the end of one chromosome
says nothing about the start of the next, while within a chromosome it
stays linked over 200 kb. Pooled over pairs of products of one oocyst
the IBD state is correlated across junctions (φ about 0.09 expected),
because meiosis I pairs the products alike on every chromosome.
