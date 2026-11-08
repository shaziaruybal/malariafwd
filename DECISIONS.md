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

## D76 — The 0.1.0 API: one function per question, every reader takes a run, no plumbing · 2026-11-08 · Accepted
**Context.** The pre-release review tabulated the 55 exports (Q24): some overlap (two functions
answer one question), some are plumbing the user never needs, and some argue against themselves in
their own help. The author accepted the table's proposals; each export change below is one commit.
**Decision.** The changes, one row per commit:

| was | is | why |
|---|---|---|
| `arg_diversity()`, `arg_divergence()` | removed | branch statistics are no relatedness measure on the forest the model records, as their own help said; `arg_tmrca()` is (the author's answer to question 12) |
| `classify_who()`, labels `Imported` / `Indigenous` | `classify_travel()`, labels `imported (recent travel)` / `local (no recent travel)`; takes a run, and `window` is checked | it applies a travel rule, and "WHO" and "Indigenous" read as the truth, which `classify_chain()` gives (review m16); lower case like the package's other classes |
| `model_version()` | removed: a run's versions are `load_run()$malariafwd_version` and `$slim_version` | it told two model generations apart by a column every run of this release has |
| `moi_vs_relapse()`, `plot_moi_vs_relapse()` | removed | one study's figure; a host's peak `moi_active` is in `load_state()`, its relapses in `classify_cases()` |
| `pick_host()` | removed: `ev[event_type == "hypnozoite_activation", .N, by = human_id][order(-N), human_id[1]]`, the host with the most relapses | one line that picks a figure's host, and nothing in the package calls it; the vignette and `quickstart.R` write it out |
| `as_slim_args()`, `slim_path()` | internal | plumbing: `run_sim()` builds SLiM's arguments, and `check_setup()` shows the SLiM path |
| `theme_malariafwd()` | internal | every plot applies it; a user restyles a plot with ggplot2's own themes |
| `relapses_per_person_year()` | internal: `summarise_run()$relapses_per_py_at_risk`, over its window | one summary of a run; `?summarise_run` defines the person-years at risk |
| `incidence_by_origin()` | `force_of_infection(detected_only =)`: its `n` is the incidence, in every window and origin | two functions counted the same cases by window and origin; the merged one has D62's rows and `rate` (review n57) |
| `reservoir_trajectory()` | `host_state_prevalence()$hypnozoites_per_host`, the day's mean on each of its rows | one table of the hosts' states by day, the reservoir among them, named as `summarise_run()` names it |
| `ibd_relapse_origin()` | removed: the `relapse vs its origin` rows of `ibd_truth_pairs()`, or `ibd_lookup()` on each relapse and its batch, `parent1_pedigree_id` | it checked a property that holds by construction, against an earlier primary of the lineage rather than the batch the relapse woke from |
| `new_infection_diversity()$clonal_frac`, a percentage | a fraction, from 0 to 1 | every other `_frac` column is a fraction (D69) |
| `ibd_trajectory()` | `new_infection_diversity(bin =)`, with `from`, `to` and `min_pairs`: a row for every window, its `days` and `n_pairs` | one question, over the run or by window; it counts new infections only, where the trajectory counted relapses too, and every window keeps its row |
| `plot_host_trajectory(state, events, host_id)` | `plot_host_trajectory(run, host_id)`, its x axis to `SIM_DAYS` | it read both of a run's tables, so it takes the run, as the other plots do; its summary layout takes the package theme (review n11) |

**Rationale.** A user looking for how to answer a question should find one function, and a function
whose help tells the user not to use it should not be exported.
**Alternatives rejected.** *Keeping removed functions internal behind `:::`* -- no caller exists, and
internal copies would still need maintaining.
**Consequences.** Runs are unchanged. Code that called a removed name stops with R's "could not find
function"; the table above says what to use instead.
**How to verify.** `tests/testthat/test-exports.R` lists the exports exactly.

---

## D75 — Genotype data are painted on the root genomes of the tree sequence · 2026-11-08 · Accepted
**Context.** The model is a mutation-free scaffold, so its tree sequence holds no sites, while the
methods a run is meant to test read marker genotypes. Every genome descends, at each position, from
exactly one root genome, a founder template or an import (D17, D18), and root genomes are unrelated
(review M23, Q23).
**Decision.** `simulate_genotypes()` takes biallelic loci as PvP01 chromosome, 1-based position and
PLAF (the ALT frequency). Each root genome draws its allele at each locus independently, and each
genome carries its root genome's allele; the truth is the root genome, by pedigree id, behind every
genome's allele at every locus. Samples are clones (0/1) or host samples grouped by `sample_id`
(0/1, 0.5 for both present). An error replaces a call (a clone call flips; a host call becomes one of
the other two at random), then a missing call is `NA`. The root genomes' alleles are drawn first,
with R's generator, for every root genome of the run, and `seed` seeds one call without touching the
session's stream. `write_vcf()` writes VCF 4.2 with one contig per PvP01 chromosome, haploid clone
calls, unphased diploid host calls and the painting frequency as `PLAF` in `INFO`. The author chose
painting (Q23); the host-sample error model and the `sample_id` grouping are the design's defaults.
**Rationale.** In a neutral, mutation-free model, painting the root genomes is exact; PLAF is what
users have from a panel; pedigree ids survive simplification; drawing for every root genome makes a
genome's alleles independent of the company it is asked for in.
**Alternatives rejected.** *msprime.sim_mutations()* -- it never puts a mutation on a root, so no
site separates two root genomes (0 of 1,657 and 0 of 16,806). *Recapitation, then mutations* -- a
coalescent prior, not the panel's PLAF; pyslim refuses roots minted mid-run, and msprime needs a
demography for 1,298 populations. *Painting in tskit and reading `genotype_matrix()`* -- by default
it reads a genome alone under its root as missing wherever the root is REF (2% of case calls), and
`ts.write_vcf()` writes one contig in concatenated coordinates. *Drawing only for the root genomes
asked for* -- a genome's alleles would depend on the others asked for.
**Consequences.** Runs are unchanged; two new exports. Painting has no background relatedness or
linkage disequilibrium among root genomes, gives imports the local frequencies, and records which
alleles a host holds, not their proportions. A run with `RETAIN_FULL_HISTORY = FALSE` limits the
genomes to those alive at its end.
**How to verify.** `tests/testthat/test-genotypes.R`: every genome carries its root genome's
allele, by the bookkeeping, by painting the root nodes in tskit, and by walking the unsimplified
genealogy; root genomes alone are their own roots; the same root genome exactly where
`ibd_segments()` finds a common ancestor; the root genomes and their days; the seed's invariance;
host samples; the VCF; refusals; the schema. `inst/validation/audit_events_and_ibd.R`: four checks on
each audit run.

---

## D74 — An IBD segment is a maximal stretch from one common ancestor, cut at the chromosome junctions; the pair fraction counts every shared stretch · 2026-11-08 · Accepted
**Context.** `ibd_pairs(segments = TRUE)` returned tskit's segment rows as they came, which are
pieces of the paths through the trees: how a pair's shared ancestry was split into rows depended on
the genomes the tree was simplified to (4,154, 2,493 or 2,394 rows for one set of 400 case pairs),
and 3% of them ran across a chromosome junction, since in the tree sequence a stretch continues
across one when both genomes carry the same parent on both sides. `ibd_segments()` cut at the
junctions but kept the path pieces, so the two functions gave different rows for one pair. Both
applied `min_span` to tskit's rows, before any cut, and the pair fraction inherited that floor, so
it too depended slightly on the other genomes passed. `ibd_segments()` numbered its pairs from 0,
needed pandas, and gave a genome with itself one row over the whole concatenated genome (review
M14, m32, m38, m39; Q17).
**Decision.** A segment is a maximal stretch over which two genomes descend from one common
ancestor genome, cut at each chromosome junction: abutting rows with one ancestor (the node, not
its time) merge, and a change of ancestor starts a new segment even with no gap. `min_span` drops
the pieces shorter than it after the cut, inclusively. A pair gets the same segments whatever other
genomes are passed, and both functions give the same rows, with the chromosome, both coordinate
systems, the length, the ancestor and its age. The pair fraction counts every shared stretch, of
any age and any length, and `min_span` does not change it (the author's answer). `ibd_segments()`
reads only its listed pairs, numbers them from 1, and gives a genome with itself one row per
chromosome. Nothing imports pandas, so the Python requirement is tskit (>= 1.0) and numpy (the
author's answer).
**Rationale.** A segment that depends on which other genomes were asked for, or that runs across
chromosomes, cannot be compared, drawn or counted.
**Alternatives rejected.** *tskit's rows as they are* -- path-dependent. *Merging on the ancestor's
time* -- two ancestors born on one tick would merge. *A floor on the pair fraction* -- tskit's floor
depends on the other genomes passed, and the segments' floor costs 4-5 times the time and memory.
**Consequences.** Runs are unchanged. On a 150-host 548-day run with co-transmission (seed 1),
1,242 of 135,981 pair fractions move, by at most 3.6e-5. Segment rows are merged by ancestor, cut
and floored after the cut; `chrom`, `start`, `end`, `gstart`, `gend` and `length` replace `left`
and `right`, and `ibd_segments()$pair` counts from 1. `ibd_segments()` for 300 pairs of that run
takes 0.06 s, against 0.66 s.
**How to verify.** `tests/testthat/test-ibd.R`: merging, cutting and flooring on a synthetic table;
on a run, the same segments from any set of genomes and from both functions; no segment across a
junction or below `min_span`; the segments adding up to each pair's IBD; one schema for every
pair; a genome with itself; pairs counted from 1.

---

## D73 — Import classes come from the true transmission chain · 2026-11-08 · Accepted
**Context.** The WHO terminology classes a case as imported, introduced or indigenous by where it
was acquired and its chain of transmission from an imported infection. The package's only import
classification, `classify_who()`, applies a travel rule, what surveillance sees, and has no truth
to score it against. The analysis repository's `chain_classify()` walked each case's parents by
hand and lost the chain at a recrudescence, so every recrudescence and every case downstream of
one came out of a local lineage (review m16).
**Decision.** `classify_chain()` gives every case (a row of `classify_cases()`) its generation, the
fewest transmission links from a clone acquired by importation, through both routes (D70), a
recrudescence counted as the clone it came from and an outcrossed oocyst's shorter chain counting;
the import row the chain starts from (the batch, for an imported relapse; the smaller id on a
tie); and its class: `imported` (generation 0), `introduced` (1), `indigenous (from an import)` (2
or more) or `indigenous (local lineage)` (no import in the chain). It searches breadth first from
every import at once.
**Rationale.** The true class is what a travel rule, or a genetic classifier, is scored against.
**Alternatives rejected.** *A walk up each case's parents* -- the same classes once fixed, slower,
and where the recrudescence bug lived.
**Consequences.** Runs are unchanged. A new export. On a 150-host 548-day run with imports and
recrudescence (seed 2), 61 of 698 cases are classed differently from the analysis repository's
walk: 26 recrudescences and 35 cases downstream of one.
**How to verify.** `tests/testthat/test-imports.R`: the fixture's classes, generations and import
rows, an imported relapse from its batch and an import's recrudescence; an empty log; on a run,
the imported class exactly the imported cases, and every generation one link from the one before.

---

## D72 — Transmission chains follow both routes, keep a recrudescence with its clone, and are columns of the links · 2026-11-08 · Accepted
**Context.** `transmission_chains()` took the components of the direct links only, so a relapse in
another host never joined its infector's chain; it started a new chain at every transmitting
recrudescence, which is the same infection as the clone it came from (D25); and its union-find,
without path compression, was quadratic in a chain's length, while liver links make chains of
thousands of genomes. It returned sizes alone, which the link table could not be joined to (review
M9; Q15, and the export table's merge of chains into the links).
**Decision.** A chain is a connected component of the links, each taken in either direction, with
each source taken as the clone it came from through repeated recrudescences. `transmission_links()`
gives each link its `chain_id`, numbered in the order of the chains' first links, and the chain's
`chain_size`, its genomes, founders included. Chains follow `kind`. A genome with no link is in no
chain. The components come from a union-find over `match()` indices with path halving.
`transmission_chains()` is removed: `unique(links[, .(chain_id, chain_size)])` gives its table.
**Rationale.** Liver links are two-thirds of transmission, and a recrudescence's onward
transmission belongs to its clone's chain. On the links, a chain can be followed and joined.
**Alternatives rejected.** *Label propagation with pointer jumping* -- it hooks nodes, not roots,
so it needs about as many passes as a chain is long: 292 s for a shuffled path of 100,000 genomes,
against 0.04 s. *Singleton chains for unlinked genomes* -- an edge table cannot hold them;
`classify_chain()` labels every case with its introduction.
**Consequences.** Runs are unchanged. In a 5-year run of 200 hosts (seed 1), 466 chains (the
largest 40 genomes) become 3 (the largest 3,886); the direct links alone give the earlier 466.
The links of that run take 0.02 s.
**How to verify.** `tests/testthat/test-transmission.R`: the fixture's chains for each kind, the
recrudescence kept with its clone and the outcrossed oocyst joining two chains; `.link_components()`
on a long path, disjoint and star graphs and 200 random graphs against a reference union-find; an
empty log; on a run, chains partitioning the linked genomes.

---

## D71 — Rc is the case reproduction number by true origin, each secondary case counted once · 2026-11-08 · Accepted
**Context.** `rc_by_origin()` counted only the direct links, so the relapses a case's bites caused
in other hosts, two-thirds of all transmission in a 5-year run, were not its secondary cases. It
credited both parents of an outcrossed oocyst in full, so a target counted twice, gave rounded
values under an undefined percentage, could not leave out the cases near the end of the run that
had not had time to transmit, and compared itself with an Rc that is not the same quantity (review
M9, m41, m82, m88; Q15).
**Decision.** Rc is the case reproduction number (Fraser 2007): the mean number of secondary cases
each case caused, by the case's true origin. A case is a clone onset (a row of `classify_cases()`);
its secondary cases are the targets of its links through both routes (D70), each counted once, 1/2
to each parent of an outcrossed oocyst (the author's answer). `kind` chooses the routes. The
function takes the run: `from` and `followup` choose the onsets counted, from day `from` to
`SIM_DAYS - followup`, both 0 by default. Every origin has a row; the columns are `n_cases`,
`n_secondary`, `Rc` and `onward_frac`, unrounded, `NA` with no case. `unit = "acquisition"` walks a
recrudescence back to its clone through repeated recrudescences. The help cites Fraser (2007) and
no longer compares the acquisition reading with Champagne et al. (2022), whose Rc is a model's
next-generation number under case management (the author's answer); this supersedes D25 item 3's
Champagne sentence.
**Rationale.** A case reproduction number counts what each case actually caused. Half credit keeps
each secondary case once, so that Rc is near 1 at a steady state.
**Alternatives rejected.** *Full credit to both parents* -- it counts an outcrossed secondary case
twice. *A default follow-up* -- the right one depends on the run's phase, which the help explains.
*An events table and a separate run length* -- the run gives both.
**Consequences.** Runs are unchanged. Rc roughly triples: in a 5-year run of 200 hosts (seed 1)
local reinfection goes from 0.355 to 0.966 and relapse from 0.354 to 0.994; over a steady window
(`from = 730, followup = 365`) the pooled Rc is 1.056, and 0.344 through direct links alone.
`followup` alone keeps the growing phase: 0.98, 1.07, 1.14 and 1.22 at 0, 180, 365 and 730 days.
The signature is `rc_by_origin(run, ...)`; `n_sources` is `n_cases`, and `pct_onward` (0-100) is
`onward_frac` (0-1).
**How to verify.** `tests/testthat/test-transmission.R`: the fixture's values for both units, half
credit and the recrudescence walked back; the routes add up; a run with no case; `from` and
`followup`, and their refusals; on a run, the weights of the links from cases and every case
counted.

---

## D70 — Transmission links include the relapses a bite causes in another host · 2026-11-08 · Accepted
**Context.** `transmission_links()` linked an oocyst's parents only to the blood-stage infections
its products started, so a relapse in the bitten host, from the batch the same bite left, had no
infector. Such links are most of them: 2,658 of 4,053 in a 5-year run of 200 hosts, with a median
lag from the infector's onset of 156 days against 59 for a direct link. The function returned bare
`source` and `target` pairs, without hosts, days or route, and counted both parents of an
outcrossed oocyst in full (review M9, n29; Q15).
**Decision.** A link joins each parent genome of an oocyst to each genome the oocyst's products
infected: `direct`, the bite's blood stage, or `liver`, each relapse of the batch the bite left,
which carries the batch's oocyst. Only local infections are targets, and a host's own relapse is
never linked to its earlier infection (the author's answer to Q15). Each row names both genomes and
hosts, the mosquito, the oocyst and four days (the source's onset, the feed, the bite and the
target), and a weight: 1/2 for each parent of an outcrossed oocyst (the author's answer), so a
target's weights sum to 1. `kind` keeps either route. A log cut short of the rows a link needs
gives one warning. In this commit `rc_by_origin()` and `transmission_chains()` read the direct
links with full credit, as before.
**Rationale.** A relapse in another host is a transmission delayed by the liver; leaving it out
undercounts onward transmission most for the infections whose batches relapse.
**Alternatives rejected.** *One row per target* -- it loses the second parent of an outcrossed
oocyst. *Carrying `detected` or `into_patent_lineage`* -- `classify_cases()` holds the case's
attributes, and the help gives the join. *A link from a primary to its own relapse* -- no mosquito
is involved.
**Consequences.** Runs are unchanged. `kind = "direct"` gives the earlier pairs under new names.
Every local case is the target of exactly one transmission.
**How to verify.** `tests/testthat/test-transmission.R`: the fixture's 13 links; the targets and
non-targets; weights; a self link; `kind`; empty logs; a cut log; on a run with co-transmission,
imports, recrudescence and batch-only bites, every local case a target and every link consistent.
`inst/validation/audit_events_and_ibd.R`: four checks on each audit run.

---

## D69 — summarise_run() is the one summary of a run, over one window of days · 2026-11-08 · Accepted
**Context.** Four of `summarise_run()`'s eight columns (the relapse fraction, relapses per
person-year, treatment courses and the imported share) covered the whole run, while prevalence,
COI, polyclonality and the reservoir covered its second half, so on a treated 400-day run the
relapse fraction read 0.608 at `burn_frac` 0 and 0.5 alike. The units were mixed and undocumented,
values were rounded, there was no EIR, FOI or molFOB, and the calibration was scored with the
analysis repository's own summary, under other definitions (review M11, m23, m29, m41; Q12).
**Decision.** One window, from day `round(SIM_DAYS * burn_frac) + 1` (or `from`) to `to` (by
default `SIM_DAYS`), for every column: counts of its days, means over its host-days, rates over its
person-years, the days after an extinction counting as days with no parasite (D61). One row of 28
typed columns named by their units (`_frac`, `_per_py`, `_tick`, `n_`), unrounded, `NA` where
there is nothing to divide. The EIR is the count of infectious bites (D47). The relapse rate is per
person-year at risk (D66); the EIR, FOI and molFOB are per person-year of every host; FOI counts
the new clones from bites and imports. molFOB has both D25 readings, and both relapse fractions are
given (D65). `exclude_masked` leaves D40's masked relapses out of every relapse column, and
`n_relapse_masked` counts them either way. The row records its window, its reading and how the run
ended. A run that stopped early without a record is refused unless `to` is at most its last day.
`plot_diagnostics()` takes its subtitle's numbers and its burn-in line from the same row, and its
force panel stacks the new clones only (review m17, m30, n10, n25). `new_infection_diversity()`
starts on the same first day for a given `burn_frac`.
**Rationale.** A summary pooled over a grid must mean the same thing in every column and every
cell. The bite count is exact, where the Ross-Macdonald formula needs a correction that depends on
mosquito mortality and the season.
**Alternatives rejected.** *Summing `force_of_infection()`'s rows* -- the window need not line up
with its bins. *Polyclonality in percent.* *One column per masked reading.* *A class with a print
method* -- `rbindlist()` drops it.
**Consequences.** Runs are unchanged. 28 columns replace 8, and every run's values change: the
window starts one day later, the event columns are windowed, nothing is rounded, and polyclonality
is a fraction. On the treated run (seed 1) the relapse fraction is 0.608 over the whole run and
0.640 over its second half, relapses per person-year at risk 5.77 and 5.93, courses 125 and 81, the
imported share 0.174 and 0.119. The glossary's FOI row, which said mosquito-borne only, gives the
FOI the functions compute, imports included. The NU = 0 validation check counts relapses, since its
window holds no onset. Section 13's figures are regenerated with the release (M33).
**How to verify.** `tests/testthat/test-summary.R`: names and classes on every kind of run; the
fixture's values over the whole run, the default window and two others; `exclude_masked`; logs
without the flag or the bite count; extinct, never-started and unfinished runs; the window's
arguments; on a SLiM run, every column against its definition over the window, and FOI and molFOB
against `force_of_infection()`.

---

## D68 — host_state_prevalence() counts subpatent carriers, and has every state every day · 2026-11-08 · Accepted
**Context.** `host_state_prevalence()` gave a row only for the (day, state) pairs that held a host,
so its stacked-area plot joined the rows on either side of an empty state, drawing hosts on days
that had none (the bug class of D62). It ignored `moi_subpatent`, so a host whose only clone was
suppressed below patency was counted as uninfected, as section 16.2 conceded (review m20).
**Decision.** One row per day of the state table and state, 0 where no host is in it, with `cls` a
factor of five levels: `uninfected`, `active only`, `active + reservoir`, `subpatent` (a subpatent
clone and no patent one, with or without hypnozoites) and `reservoir only (invisible)`. The count
column is `n`. A state table without `moi_subpatent` counts no subpatent host. The plot keys a
named five-colour palette by level.
**Rationale.** A subpatent carrier is infected and can infect a mosquito; counting it as uninfected
understates carriage exactly where treatment fails. A missing row is not a zero.
**Alternatives rejected.** *Subpatent carriers with hypnozoites in the reservoir-only state* -- they
carry a blood-stage clone, which a sensitive assay would find.
**Consequences.** Runs are unchanged. Every day has five rows; at the default
`P_RECRUD_GIVEN_FAIL = 0` the subpatent state is always 0. `N` is `n`.
**How to verify.** `tests/testthat/test-host-states.R`: every state on every day, the subpatent
carriers, a table without the column, an empty table, the plot. `tests/testthat/test-treatment.R`:
on a run with recrudescence, the subpatent count is every day's host rows with no patent clone and
a subpatent one.

---

## D67 — Time to first relapse is measured per hypnozoite batch · 2026-11-08 · Accepted
**Context.** `time_to_relapse()` timed each host from its first blood-stage infection, founder
seeding included, to the first relapse of any of its batches, possibly one a later bite laid down,
and ended every host's time on the last event's day. Founders carry no batch, so a founder's first
relapse always came from a later bite while its clock started on day 1: the curve mixed waiting
for reinfection into the time to relapse, and its shape depended on `N_FOUNDER_HOSTS` (review M13,
Q16).
**Decision.** One row per hypnozoite batch, timed from the day it was established to its first
relapse (`event = 1`), traced by the activation's `parent1_pedigree_id`. Without one, a batch is
censored on the day its last hypnozoite left it, by death or radical cure, or at `SIM_DAYS`
(`emptied = "censor"`, the author's answer to Q16); `emptied = "follow"` follows an emptied batch
to `SIM_DAYS`, the curve a trial reports. `exclude_masked = TRUE` passes over masked relapses,
whose hypnozoites still leave the batch (the author's answer to OQ2). The function takes the run;
its columns name the batch (`batch_pedigree_id`, `batch_tick`) and say why its clock stopped
(`end_reason`). `plot_relapse_km()` draws the fraction of batches not yet relapsed from (0, 1), and
stops with a message on a run with no batch. `.hypnozoite_exit_events()` names the rows other than
an activation that take hypnozoites out of a batch.
**Rationale.** A relapse comes from a batch, and the time from the infection that laid a batch down
to its first relapse is the quantity cohort and trial data report.
**Alternatives rejected.** *A host-level time from the first batch-establishing infection* -- a
host's first relapse can still come from a later batch. *Following emptied batches by default* --
the answer was to censor them; following stays available.
**Consequences.** Runs are unchanged. On a 100-host 548-day run (seed 1) the median time to first
relapse is 25.5 days per batch (200 of 271 batches relapsed), against 60 days per host (81 of 98
hosts); with masked relapses passed over it is 71 days. Section 13's latency row is re-measured
per batch: median 98 days (temperate) against 23 (tropical), where the host-level times gave 142
against 46. The function takes a run, not an events table.
**How to verify.** `tests/testthat/test-relapse-metrics.R`: the fixture's five batches, followed and
with the masked relapse passed over; the curve's steps; a run with no batch; on an import cohort,
every batch timed from its establishment and the batches radical cure reached first emptied; an
emptied batch in a run that went extinct; every liver-stage row that names a genome is an exit.

---

## D66 — Relapses per person-year count every relapse over the person-time of hosts ever infected · 2026-11-08 · Accepted
**Context.** `relapses_per_person_year()` took the events alone, so it ended every host's time on
the last event's day, not the run's. It started a host's time at its first blood-stage primary,
founder seeding included, and left out a host with no blood-stage primary, and its relapses (56
of 438 at `P_PRIMARY = 0.5` in a 100-host 548-day run). Section 16.3 defined the metric over every
host instead (review M12, Q14, m22).
**Decision.** Relapses per person-year are every relapse over the person-time of the hosts ever
infected. A host's time starts the day after its first establishment, a primary or a hypnozoite
batch from a bite or an import, and runs to `SIM_DAYS`, also after an extinction; founder seeding
does not start it (the author's answer to OQ4). `relapses_per_person_year(run, from, to)` takes
the run and an optional window, which restricts both the relapses and the person-time, and gives
`n_relapse`, `n_hosts`, `person_years` and `relapses_per_py`, unrounded, `NA` with no person-time.
An internal core takes the events and the run's length and can leave the masked relapses out of
the numerator, for `summarise_run()`.
**Rationale.** A host that has never been infected cannot relapse, and a host whose first
infection left only a batch can. The run's length, not its last event, is the period observed.
**Alternatives rejected.** *Activations over `N_HOSTS` times the run's years* (section 16.3) -- it
counts hosts that were never at risk, so the rate falls with the share never bitten. *Every
ever-infected host over the window's whole length* -- it counts a host's days before its first
infection.
**Consequences.** Runs are unchanged. The value changes on every run: 6.12 to 6.80 per person-year
on a 100-host 548-day run (seed 1), 5.64 to 5.51 at `P_PRIMARY = 0.5`, where 438 relapses are
counted, not 382. The function takes a run, not an events table. Section 13's band table states its
own definition (relapses after burn-in per year of the hosts ever infected) and is re-measured at
the next band re-check.
**How to verify.** `tests/testthat/test-relapse-metrics.R`: the fixture's relapses and person-time,
with and without a window; refusals; a run with no infection; on a run every activation counted,
and person-time to `SIM_DAYS` after an extinction.

---

## D65 — relapse_fraction() gives relapse_frac_episodes and relapse_frac_acquisitions · 2026-11-08 · Accepted
**Context.** `relapse_fraction()` returned `relapse_fraction` over acquisitions and
`relapse_frac_episodes` over every case: one reading named for the function and one for its unit,
where D25 had named both by unit. Its second count, `n_reinfection`, held imports too, and a run
with no case gave 0 for both fractions, a value it had not measured (review M12, Q13).
**Decision.** The fractions are `relapse_frac_episodes` (relapses over every case, recrudescence
included) and `relapse_frac_acquisitions` (relapses over imported, local reinfection and relapse
cases), both counted per clone onset, each `NA` when its denominator is 0; the second count is
`n_new_infection`. The acquisitions reading is the unit of the published 0.55-0.85 band. Whether
masked relapses are counted is the caller's choice (`exclude_masked`): the help says ν was fitted
with every relapse counted (0.66 at ν = 3) and that leaving them out gives 0.39, below the band,
and names neither as the band's own reading (the author's answer to OQ1). `summarise_run()`'s
`relapse_frac` is the acquisitions reading; the `plot_diagnostics()` subtitle shows the episodes
reading. This supersedes D23's clause that made `relapse_fraction()` one reading over every case.
**Rationale.** Both readings are needed once recrudescence is on, and a name that says its unit
keeps them from being swapped. A fraction with nothing to divide is not 0.
**Alternatives rejected.** *Rename the per-onset readings `_cases`* -- the episode names stay, defined
per clone onset (the author's answer). *Name the masked-out reading the band's* -- below the band at
every ν tried, while the fit counted every relapse.
**Consequences.** Runs are unchanged. Two columns are renamed and one count is renamed; values change
only on a run with no case, which gives `NA` (and `summarise_run()$relapse_frac` `NA`). The two
readings are equal on any run with no recrudescence, as at the defaults.
**How to verify.** `tests/testthat/test-relapse-metrics.R`: names and values on the fixture, with and
without the masked relapse; a run with no case; a run and its events table alike.
`tests/testthat/test-treatment.R`: on a run with recrudescence, each reading equals the share of
relapses among its cases, and the acquisitions reading is the larger.

---

## D64 — A case is a clone; episodes and recurrences are built in the package · 2026-11-08 · Accepted
**Context.** `classify_cases()` gave one row per clone onset without saying which rows belong to one
episode, and every analysis paired recurrences with earlier episodes by hand, pairing rows: a bite
that delivered two genotypes gave two rows on one day, which a row-based pairing paired with each
other. Clone spans were rebuilt in three places (the trajectory plot, two analysis scripts), each
missing some way an episode can begin or end (review M10, Q18).
**Decision.** An episode is the clones that became patent in one host on one day, named by
`episode_id` (`"<human_id>:<tick>"`). `classify_cases()` gives it on every clone row. `episodes()`
gives every blood-stage clone's span from its onset to its first end, founder seeding included, from
the event dictionary's onset and end types. `recurrences()` pairs each clone of every episode after
a host's first with each clone of the earlier episodes (`prior = "all"`) or the latest
(`"previous"`), never within one episode, and reads IBD from the run's `vivax_ibd` object; founder
seeding is not an episode (the author's answer to OQ4). The host trajectory plot is drawn from
`episodes()`.
**Rationale.** The unit a recurrence study compares is the episode, and the pairing is where
hand-written analyses went wrong; one implementation, tested against `state.tsv`, removes the
copies.
**Alternatives rejected.** *One row per episode* -- the genotypes of one episode differ in origin
and in IBD, which a classifier needs. *Founder seeding as an earlier episode* -- every founder
host's first acquired episode would become a recurrence, paired with an unrelated founder genome.
**Consequences.** Runs are unchanged. `classify_cases()` gains `episode_id`. In a 100-host run 1.5%
of clone rows share a host and a day, 15% with `LAMBDA_SPOROZOITES = 2`, and a row-based pairing
loses 7 of 5,619 pairs, 67 of 3,286 with co-transmission. The trajectory plot draws a clone still
patent at the end of the run to the last day of the state table.
**How to verify.** `tests/testthat/test-episodes.R`: the fixture's episodes, spans and pairings; IBD
from the object; tables, runs and directories alike; typed empty tables; with co-transmission,
several clones to an episode, never paired with each other; and on a run, every clone's span covers
exactly the days it counts in `moi_active`.

---

## D63 — ibd_pairs() returns a pedigree-keyed vivax_ibd object; unknown ids read NA · 2026-11-08 · Accepted
**Context.** `ibd_pairs()` returned a bare list keyed by the node ids of a tree it had simplified
and did not return, so every lookup went through `ibd_index()`, a list holding a closure that
printed its own source. The genomes' node times, which D22's recency filter needs, were not
returned, and an id the tree sequence did not hold scored 0, the same as "unrelated" (review M14).
The Python ran as code pasted into `__main__`, with tables passed through CSV files (m57, m58).
**Decision.** `ibd_pairs()` returns a `vivax_ibd` object: `pairs` (every pair of the genomes found,
`a < b`, zeros included), `samples` (`pedigree_id`, `node_time`, `tick`), `missing`, `segments`
(with `mrca_id`, the common ancestor's pedigree id, and `depth`, its age measured from the pair),
and the run's `end_tick`. `ibd_lookup()` reads pairs (1 for a genome with itself, `NA` for one not
in the object), and `as.matrix()` and `as.data.frame()` give the table. Requested ids the tree
sequence does not hold warn and are listed; asking only for such ids is an error, which gives the
reason the run's record holds (it never started, or `RETAIN_FULL_HISTORY = FALSE`), and the readers
built on it, `new_infection_diversity()` and `ibd_truth_pairs()`, give their `NA` row or empty table
instead. `ibd_index()` is removed. The Python side is one private module, called with arguments and returning arrays.
**Rationale.** An object keyed by what callers hold, pedigree ids, needs no translation, and `NA`
for an unknown genome cannot be mistaken for an unrelated one.
**Alternatives rejected.** *A deprecated `ibd_index()` wrapper* -- nothing has been released. *A
sparse table with implicit zeros* -- every reader would have to fill them. *0 for absent ids* -- the
confusion this removes. *Returning the simplified tree* -- a live Python object, which cannot be
saved.
**Consequences.** Values are unchanged, except a few that the CSV round trip had moved by one unit
in the last place (at most 5.6e-17); the design's measurement put it about three times as fast (0.35
against 1.02 s for 800 genomes). Code that used `ibd_index()` or `getibd()` changes, six files of
the analysis repo among it. `ibd_pairs()` no longer uses pandas.
**How to verify.** `tests/testthat/test-ibd.R` (reading the object back; its shape on a run; lookups
of relapses against their origin; segment ancestors and ages; ids not in the tree sequence; a run
with no parasite) and the IBD checks of `inst/validation/audit_events_and_ibd.R`.

---

## D62 — Windowed series have a row for every window of the run · 2026-11-08 · Accepted
**Context.** `force_of_infection()` had rows only for the (window, origin) cells with a case, and
the dashboard and `plot_force_of_infection()` zero-filled only over windows that already had a row,
so a window with no case vanished and the lines and areas were drawn straight across it: the
artefact D21 set out to remove, at the scale of whole windows. A run with no case gave no rows, and
its plot failed. The first window, which holds days 1 to `bin - 1`, was divided as `bin` days
(review M15, n6).
**Decision.** `force_of_infection()` gives one row per window of the run, to the window holding
`SIM_DAYS` (D61), and per origin, with `n = 0` where there is none; `win` and `days` are integers,
and `days` is what the window covers within the run, `bin - 1` for the first. `bin` must be a whole
number of days, at least 1. `plot_force_of_infection()` draws both series in every window, and
`plot_diagnostics()`'s force panel divides each window by the same days. `incidence_by_origin()` is
left as it is, to be merged into `force_of_infection()`.
**Rationale.** A window without a case is a measurement of zero, and a series should say so rather
than leave it out.
**Alternatives rejected.** *Zero-filling only in the plots* -- every other reader of the table would
meet the gaps again.
**Consequences.** Runs are unchanged. The table gains the zero rows; the first window's rate rises
by `bin / (bin - 1)`; full windows are unchanged.
**How to verify.** `tests/testthat/test-windows.R` (every window and origin, the days of the first
and last windows, a run with no case, the plot of one, an extinct run's zeros).

---

## D61 — Run-length windows use SIM_DAYS, and days after an extinction count as days with no parasite · 2026-11-08 · Accepted
**Context.** D37 ended an extinct run on its extinction tick, and the R functions that summarise a
window of a run took its end from the last day of `state.tsv` or the last event (D21), so a run that
went extinct had a shorter window: six runs extinct by day 196 of 400 reported a prevalence of 0.005
to 0.045 for days 200 to 400, and the windows holding their extinction day were 10 to 25 days long
(review M7, Q8).
**Decision.** `summarise_run()`, `force_of_infection()`, `new_infection_diversity()` and
`plot_diagnostics()` take the run's end from `SIM_DAYS`. Prevalence and the reservoir are means over
every host-day of the window, so the days after an extinction count as days with no infected host;
COI and polyclonality stay means over infected host-days, `NA` when there are none. The dashboard's
panels run to `SIM_DAYS` on a grid of every window, at zero after the extinction, and its subtitle
names the extinction day. The per-day tables (`coi_distribution()`, `host_state_prevalence()`,
`reservoir_trajectory()`) report the days `state.tsv` holds, which is not padded.
**Rationale.** A run's window is the period the study asks about, not the period the parasite
lasted; a run that died out had no infection on its remaining days.
**Alternatives rejected.** *Padding `state.tsv` to `SIM_DAYS`* -- the file would hold rows for days
the model never ran.
**Consequences.** Runs are unchanged. A run that reaches `SIM_DAYS` gives the same values; one that
went extinct gives zero prevalence for its days after the extinction, and full-length windows.
`force_of_infection()` no longer reads `state.tsv`. The molFOB citation is corrected to Koepfli et
al. (2013) (review m89).
**How to verify.** `tests/testthat/test-windows.R`.

---

## D60 — Parameters are documented from one dictionary; an edited set is checked when it is used · 2026-11-08 · Accepted
**Context.** `?vivax_params` documented no parameter, and MODEL_DESCRIPTION section 10 missed ten of
them, with no link between the two (review M21). The facts about each parameter (its unit, the
values it may take, its source) were spread over prose, comments and hand tables that disagreed in
places; a parameter set had no table form (m55); the defaults' sync test skipped every default
written as an expression (m73). The review's export table asked for a set edited with `$<-` to be
checked at the edit, while M26 had settled on checking it when it is used.
**Decision.** One dictionary (`R/param-dictionary.R`) holds a record per parameter: the typed
default, unit, meaning, the single-value bounds, a joint rule in words and the source. The defaults
are taken from it, and it generates the Parameters section of `?vivax_params` and the tables of
MODEL_DESCRIPTION section 10 (`dev/update-doc-tables.R`, which also writes sections 16.1 and 16.2); `param_info()` gives it as a table, and
`as.data.frame()` gives a parameter set's. `validate_params()` stays the one check, and a test holds
every documented bound to it. An edited set is checked when it is used (`run_sim()`, `run_grid()`,
or `do.call(vivax_params, p)`), and `?vivax_params` says so; there is no `$<-` method. The sync test
evaluates each model default in order.
**Rationale.** A parameter's facts should be written once, where the code reads them, so the help,
the model description and the checks cannot disagree.
**Alternatives rejected.** *Generating `validate_params()` from the dictionary* -- its messages and
their order are pinned by the tests and the fuzz, and the joint rules are not single-value bounds.
*A `$<-` method that checks each edit* -- every caller that adds `SEED` or `OUTDIR` to a set, and
every test that builds a refused set on purpose, would have to strip the class first, for a check
that use already makes. *A TSV in `inst/extdata`* -- one more file to keep in step, for no reader
that R code is not. *Pointing section 10 at the help page* -- the model description is read on its
own. *A markdown table in the Rd* -- the console help does not wrap it.
**Consequences.** Runs are unchanged. `?vivax_params` documents every parameter, by group; section
10's tables are generated, so they change only with the dictionary. `param_info()` is a new export.
**How to verify.** `tests/testthat/test-param-dictionary.R` (a record and one group for each
parameter; the types are the model's; the bounds are the ones enforced, each refusal naming the
parameter; the help and section 10 hold the current tables; printing by value and the table of a
set) and the sync test in `tests/testthat/test-params.R`.

---

## D59 — A grid resumes from what is on disk, and grid.tsv records it · 2026-11-08 · Accepted
**Context.** `run_grid()` renamed every set `setNN` unless all were named, ran a set whose name
repeated into the directory of the first and reported both as ok, ran seed 2.5 as 2 and `c(1, 1)`
twice into one directory, found a refused set only cell by cell, ran every cell again on every call,
wrote no record of the grid, and advanced the caller's random-number stream (review M27, m46, n41).
**Decision.** Blank names are filled by position, and a name must differ from the others, also
ignoring case, and be usable as a directory name; seeds must be distinct whole numbers. Every set is
rebuilt, and every cell's directory checked against what `run_sim()` would do with it, before any
cell starts: a refused set, or a cell holding another run, one that did not finish, or files that
are not a run's, stops the grid with all of them listed. A cell that already holds its run is not
run again (`skip_existing = TRUE`, the default), so a second call resumes, and `overwrite = TRUE`
runs the cells of an edited set again. `grid.tsv`, written at the end with paths relative to the
grid's directory and merged with an earlier one, records each cell's outcome, versions and key, and
every parameter as SLiM was given it; `load_grid()` reads it. The result is a `data.table` with the
parameters that vary between the sets. The cells run with `furrr`'s seed off, so the caller's stream
is left as it was.
**Rationale.** A grid is launched, left, and launched again after a failure; it should cost only the
cells that are missing, and say what it holds.
**Alternatives rejected.** *A `vivax_grid` class* -- data.table keeps a subclass through subsetting,
so its methods would misreport subsets. *Keying cells on the SLiM version* -- a SLiM upgrade would
invalidate every finished grid; the versions are recorded instead. *Writing the record from the
workers* -- they race on one file. *Failing the cells whose conflicts are known before the start* --
the grid would run half its cells before saying what was wrong.
**Consequences.** A grid run again into the same directory skips its finished cells instead of
running them again; a changed set stops it until `overwrite = TRUE`. The grid's directory holds
`grid.tsv`. Cells run by other means, such as one `run_sim()` per job of a job array, are recorded
by calling `run_grid()` with the same sets and seeds afterwards. processx's supervisor, which stops
SLiM when R dies, opens two FIFOs the first time in each process; `run_grid()` muffles future's
report of them.
**How to verify.** `tests/testthat/test-grid.R` (names, seeds, the caller's stream, resuming, an
edited set, `grid.tsv` and `load_grid()` across a moved directory and a second call, a cell that
lost its run, a failed cell) and the grid tests in `tests/testthat/test-run.R`.

---

## D58 — A run is written beside its directory and moved into place when it has finished; DONE marks it · 2026-11-07 · Accepted
**Context.** `run_sim()` wrote into `outdir` as SLiM ran, so a run that failed, timed out or was
interrupted left old and new files mixed, which the loaders read without complaint; nothing marked a
finished run; a second run into the same directory replaced the first file by file; SLiM's error
output did not reach `run.log`; and SLiM kept running when R died (review M27, m60, m61). The author
chose that a directory already holding a run is refused unless `overwrite = TRUE`, or
`skip_existing = TRUE` for the same run (2026-11-07).
**Decision.** `run_sim()` writes the run into a new sibling directory, with SLiM started there
(`OUTDIR` is `'.'`) on the run's copy of the model, and renames it to `outdir` only after SLiM
exited 0, the four outputs exist and `DONE` has been written, last. A run that fails, times out or
is interrupted is moved to `<outdir>.failed`; a run being replaced is deleted only once the new one
is in place, and only run files are ever deleted. `DONE` records the time, the elapsed seconds, the
seed, a key of the parameters and seed (the md5 of the sorted arguments SLiM is given, without
`OUTDIR` and `VERBOSE`), the model's md5, the versions and every file's size. An `outdir` holding
the same run is returned with `skip_existing = TRUE`; one holding any run is replaced with
`overwrite = TRUE` and is otherwise refused (class `malariafwd_outdir_in_use`); one holding other
files is never used. `load_run()` gains `complete`, and every reader given a directory warns once
per directory per session (class `malariafwd_incomplete_run`) when `DONE` is missing or a file has
changed since. SLiM's error output goes to `run.log`, and SLiM is stopped when R dies.
**Rationale.** A run directory is ground truth: it should hold one whole run or none, say which, and
not be replaced by accident.
**Alternatives rejected.** *Deleting `outdir` first* -- a failed run then loses the one before it.
*`DONE` alone, with the run written in place* -- a failed run still leaves a mix. *Loaders that
refuse an incomplete directory* -- they would lock out direct `slim` runs and every run made before
`DONE`. *Hashing every output into `DONE`* -- sizes catch a truncated copy at no cost. *A marker
written by the model* -- the model cannot know that R has the whole run. *Replacing by default* --
re-running a script would replace ground truth without a word.
**Consequences.** Events, state and `params.tsv` are byte-identical; the tree sequence's provenance
records `OUTDIR='.'` and `model.slim`. A second run into the same directory stops unless asked, so
the tests write into `tempfile()` directories (n47). Directories made before `DONE` warn once each a
session. A killed R session leaves an `<outdir>.partial-*` directory to delete by hand; no SLiM
process remains (checked with `kill -9` of R during a run).
**How to verify.** `tests/testthat/test-run-dir.R`: a run stopped part-way leaves `outdir`
byte-identical; `outdir` appears only when the run has finished; SLiM's error output reaches
`run.log`; reuse and replacement only when asked; files that are not a run's are never touched;
`DONE`'s rows and key; an incomplete directory warns once; the directory repeats its own run; the
key in any locale.

---

## D57 — The tree-sequence functions reach Python through one check · 2026-11-07 · Accepted
**Context.** Three copies of the Python check (in `ibd_pairs()`, `ibd_segments()` and the ARG
functions) declared different modules (tskit with pandas; tskit with numpy), hid reticulate's
failures inside `try()`, suggested a `py_require()` that does nothing once Python has started, and
checked no version, so a tskit below 1.0 failed later with a Python traceback (review m45, m64,
M17). The author chose the hybrid: use the Python reticulate is bound to when it has the modules,
and let reticulate set up its own when none was chosen (Q21).
**Decision.** One check (`R/python.R`) declares tskit (>= 1.0), pandas and numpy with `py_require()`
before Python starts, which reticulate acts on only when no Python was chosen; it then starts Python
and reads each module's version, adding the modules to reticulate's own environment where they are
missing there. Every function that runs Python calls it, naming itself and the modules it uses, and
stops with class `malariafwd_no_python` and one paragraph: which Python, how it was chosen, what is
missing and what to do. A check that finds everything holds for the session. `tskit_available()`
answers without a condition, and `check_setup()` reports SLiM, the model and the Python side.
**Rationale.** Declaring changes nothing for a Python the user chose and is what reticulate needs to
build its own; one module set and one message replace three that disagreed.
**Alternatives rejected.** *Declaring in `.onLoad()`* -- it would add tskit to reticulate's
environment for a user who never calls an IBD function. *A bare `py_module_available()`* -- it is
FALSE before Python has started. *One Python module in `inst/python`* -- deferred to v0.2 with the
answer.
**Consequences.** Which Python is used does not change. A tskit below 1.0, or a missing numpy, is
refused with the message instead of failing later, and the error texts change (code matching the old
text should match the class). The validation scripts and `quickstart.R` use `tskit_available()`.
**How to verify.** `tests/testthat/test-python.R`: each branch mocked (no reticulate, an old one, a
Python that does not start, missing modules, an old tskit, reticulate's own environment, warnings
while starting, a function that needs only some modules, every Python-backed export stopping with
its name), `check_setup()` with a missing or old SLiM, and on a machine that has them, the real SLiM
and Python.

---

## D56 — The event vocabulary has no seed uptake, and mosq_death counts genomes · 2026-11-07 · Accepted
**Context.** `event_types()` listed `hv_uptake_seed`, a mosquito's founder genotype, but every
mosquito starts uninfected, so the type could never occur, and the model's mosquito-seeding branches
were dead code with comments saying otherwise. `mosq_death` was described as "mosquito dies", while
the model writes one row per genome an infected mosquito carried and none for an uninfected death,
so counting its rows as deaths gave vector mortality wrong by a factor of 4 to 40 (review m24).
**Decision.** `hv_uptake_seed` and the dead branches are removed; the dictionary describes
`mosq_death` as one row per genome lost when an infected mosquito dies, and says that uninfected
deaths are not logged. MODEL_DESCRIPTION 4 says that mosquitoes start uninfected.
**Rationale.** A dictionary is read as a promise about the log; a type that cannot occur, or a
description of the wrong unit, breaks it.
**Alternatives rejected.** *Seed infected mosquitoes, so the type can occur* -- no analysis asks for
it, and the founders are hosts by design. *Log every mosquito death* -- one row per death of a
mosquito without parasites would add most of a run's rows for a number the state table already
gives.
**Consequences.** Runs are unchanged: events, state, `params.tsv` and the tree tables are
byte-identical, and only the copied model and the provenance differ. `event_types()` has one row
fewer.
**How to verify.** `tests/testthat/test-load.R` (the vocabulary) and `tests/testthat/test-run.R`
(every type a rich run logs is in the dictionary).

---

## D55 — N_OOCYSTS_PER_BITE is retired · 2026-11-07 · Accepted
**Context.** `N_OOCYSTS_PER_BITE` was kept as a tripwire after `LAMBDA_SPOROZOITES` replaced it: any
value but 1 stopped the run, yet it was still a parameter, in every set and every `params.tsv`
(review n20). It capped the lineages a bite delivered, not oocysts. Retiring it is the review's fix,
taken as the Phase B plan's default; the author can reverse it at the B1 review stop.
**Decision.** It is retired as `N_FOUNDERS` is (D54): `vivax_params(N_OOCYSTS_PER_BITE = )`, a set
edited to hold it, and a direct run that defines it are refused, whatever the value, with a message
naming `LAMBDA_SPOROZOITES`; `params.tsv` has no row for it; a saved `params.tsv`'s row is left out
without a warning, as is the row of any retired parameter.
**Rationale.** A parameter that can take one value only sets nothing.
**Alternatives rejected.** *Keep the tripwire* -- one more name in every set and every record, for a
value that cannot change.
**Consequences.** Runs are unchanged (`params.tsv` loses one row). No analysis script sets it.
**How to verify.** `tests/testthat/test-params.R` and `tests/testthat/test-slim.R` (refused from R
and with `-d`, at the value 1 too) and `tests/testthat/test-load.R` (a saved row is left out
quietly).

---

## D54 — N_FOUNDERS is retired; the founders are N_FOUNDER_HOSTS and N_FOUNDER_LINEAGES · 2026-11-07 · Accepted
**Context.** `N_FOUNDERS` survived from before the founders were split into hosts infected at the
start (`N_FOUNDER_HOSTS`) and the distinct genomes among them (`N_FOUNDER_LINEAGES`): it only set
the default of `N_FOUNDER_HOSTS`, so an edit to it in a built set did nothing, and the documentation
called it legacy (review m48, n32). The author chose to retire it (2026-11-06).
**Decision.** `vivax_params(N_FOUNDERS = )`, a set edited to hold it, and a direct run with
`-d N_FOUNDERS=` are refused with a message that names `N_FOUNDER_HOSTS` and `N_FOUNDER_LINEAGES`.
`N_FOUNDER_HOSTS` defaults to 12 and `N_FOUNDER_LINEAGES` to `N_FOUNDER_HOSTS`; the message says
that editing `N_FOUNDER_HOSTS` in a set already built leaves `N_FOUNDER_LINEAGES` as it was, so both
should be set there. A set built by an earlier version holds `N_FOUNDERS`, so `run_sim()` tells such
a set how to drop it (`p$N_FOUNDERS <- NULL`), and `load_params()` reads a run object saved by an
earlier version again from its directory. `params.tsv` has no `N_FOUNDERS` row. A saved `params.tsv`
that has one still loads: the row is ignored beside an `N_FOUNDER_HOSTS` row, and otherwise, in a
file made before the founders were split, gives `N_FOUNDER_HOSTS` its value, with
`N_FOUNDER_LINEAGES` following, without a warning. A retired name gets its own message, so a later
retirement is one more entry.
**Rationale.** A parameter whose edit does nothing is a trap, and two names for one count invite
sets that disagree.
**Alternatives rejected.** *Keep it as an alias, re-deriving `N_FOUNDER_HOSTS` from it on every
edit* -- two names for one value, and the rule for which one wins is one more thing to learn.
*Refuse saved files that have the row* -- every saved run has it, and the value is recorded beside
it.
**Consequences.** Runs are unchanged (events, state and trees; `params.tsv` loses one row). Code
that sets `N_FOUNDERS` stops with the message. Renaming it to `N_FOUNDER_HOSTS` gives the same run;
where the code set `N_FOUNDER_HOSTS` too, as every analysis script did, deleting it does.
**How to verify.** `tests/testthat/test-params.R` (the refusal, with another unknown name, from an
edited set and from a set of an earlier version, and `N_FOUNDER_LINEAGES` following
`N_FOUNDER_HOSTS`), `tests/testthat/test-slim.R` (`-d N_FOUNDERS=12` refused),
`tests/testthat/test-load.R` (an older file's rows, and a run object of an earlier version) and
`tests/testthat/test-run.R` (no `N_FOUNDERS` row in `params.tsv`).
**Context.** Most rows of `state.tsv` are mosquitoes (100 of every 140 at the defaults), whose
`moi_active` counts the broods they carry, so a prevalence taken over every row,
`mean(moi_active > 0)`, gave 0.30 against a host prevalence of 0.54 in the review's run (m28). No
function of the package needs the mosquito rows. The author chose the host rows as the default
(2026-11-07).
**Decision.** `load_state(run, entity = c("human", "mosquito", "all"))` keeps the rows of the entity
asked for, the hosts by default; `"all"` gives both, in the order of the file. The columns and their
classes do not depend on `entity`.
**Rationale.** The default should be the rows nearly every analysis wants, so that the obvious
summary is the right one.
**Alternatives rejected.** *Every row by default* -- it keeps the mistake easy. *Separate loaders
for hosts and mosquitoes* -- two functions for one file, with one schema.
**Consequences.** Runs are unchanged. Code that used the mosquito rows of `load_state(run)` gets
none until it passes `entity = "all"` or `"mosquito"`; the package's own readers filter hosts
already and read less.
**How to verify.** `tests/testthat/test-load.R` (each `entity` gives its rows, with the same
classes; another value is refused) and `tests/testthat/test-run.R` (a run's host and mosquito rows a
day).

---

## D52 — load_params() returns the parameter set a run was made from · 2026-11-07 · Accepted
**Context.** `load_params()` returned every value of `params.tsv` as a string, with the record rows
(`ACTUAL_SEED`, `SEED`, the versions, `GENOME_LENGTH`, `N_CHROMOSOMES`, `RECOMB_WITHIN`) among them,
so `run_sim(load_params(run))` failed on the names it did not know, and comparisons of its values
were comparisons of strings (review M16). `load_run()` accepted any directory, and gave a run with
seed `NA` for one with no outputs (m49). The package's `%||%` masked base R's, with other semantics
(n19), and `chrom_map()` parsed `CHROMOSOMES` its own way (m34).
**Decision.** `load_params()` returns a `vivax_params` set: every parameter, in the defaults' order
and types, read from the text the model wrote (T or F, numbers stored by D49's rule, strings as they
are), and not validated. The record rows are typed fields of the `load_run()` object: `seed`,
`slim_version`, `malariafwd_version`, `genome_length`, `n_chromosomes` and `recomb_within`, beside
`extinct_tick` and `never_started`. A file that lacks parameters, or holds rows that are neither
parameters nor record rows, loads with the defaults for the first and without the second, and warns
once per directory per session (class `malariafwd_params_mismatch`, under `malariafwd_run_warning`).
`load_run()` refuses a directory without `params.tsv` or `events.tsv`. `chrom_map()` parses
`CHROMOSOMES` with `validate_params()`'s parser and stops when the chromosomes do not add up to the
recorded `GENOME_LENGTH`. `%||%` is gone. Printing a run shows its size, outcome, non-default
parameters and versions (m55).
**Rationale.** A run's record should rebuild the run, and the set `run_sim()` takes is the form that
does; the record rows describe the run rather than set it.
**Alternatives rejected.** *The record rows in an attribute of the set* -- most list operations drop
an attribute, and the run object is where a run's facts are. *Validating on load* -- a reader must
not fail on a file's content, and `run_sim()` checks a set before it runs. *Filling missing
parameters silently* -- every saved run of an earlier version lacks some, and the values filled in
are this version's defaults, not settings of that run.
**Consequences.** Runs are unchanged. `run_sim(load_params(run), seed = run$seed)` repeats a run
started from R exactly when each of its non-zero floats lies between 1e-130 and 1e15 in magnitude.
Code that read a record row from `load_params()` (`p$ACTUAL_SEED`) reads the run instead
(`run$seed`), and `as.numeric()` on a parameter does nothing. Each saved run of an earlier version
warns once a session until it is made again.
**How to verify.** `tests/testthat/test-run.R` (the set read from a run is the set it was made from,
and a run from it is byte-identical; the record is typed; the test helper writes `params.tsv` as the
model does) and `tests/testthat/test-load.R` (types, unreadable files, the warning given once, the
record rows, a directory that is not a run).

---

## D51 — A host's clones of the day count toward MAX_MOI, and do not stop a bite · 2026-11-07 · Accepted
**Context.** A blood-stage clone a host gains in `reproduction()`, by a relapse or a bite, is not in
the host's individuals until that callback is over, while the host is tagged infected at once
(D47). A further infectious bite on the host the same day went to the branch for an infected host
and an infected mosquito, which skipped a host with no visible transmissible clone, so the bite was
dropped both ways: 0-1% of infectious bites in untreated runs and 2.8-5.5% with case management
(`BLOODCURE_COVERAGE` 0.8), where hosts whose only clone is that day's relapse are common. The
`MAX_MOI` check (D38) counted only visible clones, so same-day bites took a host past it: at
`MAX_MOI = 2` with no relapses, blood-stage MOI reached 4. The author chose to fix it before
Phase B (2026-11-07).
**Decision.** Each host keeps the lineages of the blood-stage clones it gained that day, cleared at
the start of each day's `reproduction()`, and a bite's `MAX_MOI` check counts them with the clones
the host carries. A bite on a host tagged infected is inoculated as any other: only the mosquito's
uptake needs a transmissible clone, and it is skipped, with no draw, when the host shows none yet.
**Rationale.** The order of the bites within a day is an artefact of the loop; it should decide
neither whether a bite counts nor whether the cap holds.
**Alternatives rejected.** *Guard each direction of the branch separately, and nothing else* -- the
cap cannot see same-day clones either, so it would let a second clone in on the same day even at
`MAX_MOI = 1`. *Make the day's clones visible to the day's bites* -- SLiM adds a callback's new
individuals only once `reproduction()` is over.
**Consequences.** A run changes from the first infectious bite on a host tagged infected that shows
no transmissible clone, or on one that the day's clones have filled to `MAX_MOI`; a run with neither
is unchanged. With every infectious bite infecting and no host full, a host's inoculations a day are
its infectious bites. A clone a host gains during the day is first taken up by a mosquito the next
day, even at `GAMETOCYTE_DELAY = 0`. A founder host's copy of another founder's genome arrives on
day 2 after that day's bites, so the check does not refuse it a day-2 bite, as it does not refuse a
day-1 import (D39); relapses still bypass the check, and a bite that delivers several genotypes is
admitted whole. Re-checked, ν = 3 is still inside all four bands (D48 and MODEL_DESCRIPTION 13 give
the new figures); MODEL_DESCRIPTION's other measured figures, the representative run among them, are
refreshed before release (review M33).
**How to verify.** `tests/testthat/test-model-invariants.R`: at `P_TRANS_VH = 1` with case
management and no host full, each host's infectious bites a day are its inoculations that day; at
`MAX_MOI = 2` with no relapses and one genotype per bite, no host ever carries more than 2
blood-stage lineages.

---

## D50 — The loaders give every run the same columns, in the same classes · 2026-11-07 · Accepted
**Context.** `load_events()` and `load_state()` let `fread()` guess each column's class from the
values, so the classes drifted with what a run held: a header-only `events.tsv` (a run with no
event) loaded as 21 logical columns, a run of only no-take rows had its lineage columns logical,
and `lineage_ids` was integer in an extinct run and character in others. Code that joined or
compared these columns failed on some runs only, `moi_vs_relapse()` among them (review M15).
**Decision.** Both loaders read their file with the columns and classes the model writes, held in
one schema per file (`R/outputs.R`): a column the file lacks is added as `NA` of its class, a
column the schema lacks follows in the file's class, and a value that does not fit stops with the
file, the column and the line named. The `lineage_ids` column of older state files, named
`brood_ids`, is read under its present name. A missing or empty file, or an argument that is not a
run, is named in the error.
**Rationale.** A table's columns are part of its contract; a reader should not have to test classes
before using them.
**Alternatives rejected.** *Convert the -1 sentinels to NA* -- analyses test `!= -1`, and the
sentinel means 'not applicable', not 'missing'. *Refuse a file that lacks a schema column* -- every
saved run before 2026-11 lacks one of them.
**Consequences.** Runs and outputs are unchanged. A run with no parasite loads, and
`moi_vs_relapse()` gives an empty table for it. Saved runs gain `NA` columns for what they did not
record (`into_patent_lineage`, `moi_subpatent`, `infectious_bites`).
**How to verify.** `tests/testthat/test-load.R`: an empty log, content-independent classes, missing
and extra columns, a value that does not fit, a missing or empty file, and normal, extinct and
never-started runs with identical classes.

---

## D49 — params.tsv records floats to 15 significant digits, and a parameter set holds what SLiM is given · 2026-11-07 · Accepted
**Context.** `run_sim()` passes each double to SLiM with 15 significant digits, but the model wrote
`params.tsv` with Eidos's default 6, so `load_params()` read `ALPHA_ACTIVATE = 1/170` back as
0.00588235 and a set rebuilt from a run's record differed from the one it ran. `vivax_params()`
kept doubles as given, so a probability computed as `1 + 2e-16` was refused although SLiM would
have been given 1, and kept whole numbers given as doubles as doubles, so
`vivax_params(N_HOSTS = 40)` printed 40 as a non-default (review n18). The author chose 15 digits
(2026-11-06).
**Decision.** The model writes every float parameter to `params.tsv` with 15 significant digits,
keeping a decimal point on a whole value (`10.0`) and printf's exponent form (`7.4e-07`).
`vivax_params()` stores each value as its default's type: a whole number as an integer (a double
beyond the integer range), a double as the double SLiM is given (the value of its 15-digit text,
where R reads that text back to itself). The defaults are stored the same way, and
`validate_params()` works in doubles so that sums of large integer counts cannot overflow.
**Rationale.** The record should be the run, and a parameter set should hold what the run is given.
**Alternatives rejected.** *17 digits* -- it changes the text SLiM is given and so every run.
*Hexadecimal floats* -- Eidos cannot read them. *A decimal point inside exponent forms* -- R reads
`1.0e-25` and `1e-25` differently far below 1e-22. *Writing the model's own defaults as 15-digit
literals* -- the expressions (`1/170.0`) are the documented values.
**Consequences.** Runs are byte-identical: the text SLiM is given does not change. `params.tsv`
rows of floats with more than 6 significant digits change text (`ALPHA_ACTIVATE` and
`MU_HYPNOZOITE` at the defaults), and so do round exponent forms (`1.0e-06` becomes `1e-06`); R
reads both the same. A value built in the last digits of a cap's 1e-12 allowance can be refused,
since R checks the 15-digit value.
**How to verify.** `tests/testthat/test-params.R` (a double is stored as SLiM is given it; each
value has its default's type; the defaults are a parameter set; validation in doubles; a rebuilt
set is the set) and `tests/testthat/test-slim.R` (`params.tsv` records floats to 15 significant
digits, from R and from a direct run).
---

## D48 — NU_HYPNOZOITE_BATCH defaults to 3, the fitted value · 2026-11-07 · Accepted
**Context.** The default hypnozoite batch mean was 5, called a placeholder, while every documented
figure and the fit to four published vivax bands (relapses per person-year 0.5-5, molFOB 0.7-15,
median MOI 1-3, relapse fraction 0.55-0.85; malariafwd-analyses `01_calibrate_nu.R`) used 3. A
default run gave 6.8 relapses per person-year, above its band (review M32, question Q19).
**Decision.** The default is 3. `?vivax_params` says the 40-host defaults are an illustrative
preview, not a calibrated setting. MODEL_DESCRIPTION 13 gives the fit, re-measured at this version.
The author chose this, with a re-check of the bands after the other model changes (Q19,
2026-11-05).
**Rationale.** A default should be the value the documentation and its figures use, and the one the
fit supports.
**Alternatives rejected.** *Keep 5 and call it a placeholder* -- it is out of band and used
nowhere.
**Consequences.** Every run that leaves `NU_HYPNOZOITE_BATCH` at its default changes. The re-check
(400 hosts, 5 years, `BITING_RATE` 0.08, 3 seeds, re-measured after D51) puts ν = 3 inside all four
bands when every relapse is counted: relapses per person-year 2.32, molFOB 3.44, median MOI 1,
relapse fraction 0.66. With the relapses into a patent clone of their own lineage left out of both, as agreed for
the relapse-fraction band (two-thirds of relapses at ν = 3), relapses per person-year stay inside
their band up to ν = 12, but the relapse fraction is 0.39 at ν = 3 and stays below its band at
every ν up to 12 (0.48). Not settled here: which reading the band is compared with.
**How to verify.** `vivax_params()$NU_HYPNOZOITE_BATCH` is 3, as is the model's own default (the
test that keeps the two in sync); MODEL_DESCRIPTION 13's table is the re-check.
---

## D47 — state.tsv counts each host's infectious bites a day, which gives the EIR · 2026-11-07 · Accepted
**Context.** MODEL_DESCRIPTION said the package computes the EIR from `state.tsv`, with the
sporozoite rate read from `n_transmissible`. No function did, and that formula under-counts: the
snapshot is taken after the day's mosquito deaths, so it gives about 1 − 1/`MOSQ_LIFESPAN` of the
realised rate at constant biting (review m15), and with seasonality its error can go either way.
The EIR is the usual target for calibrating transmission intensity, and the author asked for it to
be part of `summarise_run()` (Q12).
**Decision.** Each host row of `state.tsv` gains `infectious_bites`: the bites the host received
that day from a mosquito carrying a product past the EIP, whatever the bite then did (`NA` on
mosquito rows). The model counts a bite as infectious exactly when
`sampleTransmissibleMosqBroods()` finds a product, which it already draws for every bite, so the
count draws no random number. The count is kept on the host and reset after each day's snapshot,
since `reproduction()` does not run on a day with no parasite. `summarise_run()` is to report the
EIR from it (Q12).
**Rationale.** Counting at the bite is exact, and per host it also gives the heterogeneity of
exposure; a snapshot formula needs a correction that depends on mosquito mortality and the EIP.
**Alternatives rejected.** *A column on mosquito rows* -- the EIR is per person, and the host row
is where a reader looks for exposure. *A separate per-day file* -- one more output to load and
document, for one number per day.
**Consequences.** `state.tsv` has one more column, last; `events.tsv`, `params.tsv`, the tree
sequence and every other column of `state.tsv` are byte-identical. In test runs the snapshot
formula gives 0.79-0.97 of the counted EIR at mosquito lifespans of 5-20 days. Writing the test
found that an infectious bite is dropped on a host whose only clones arrived earlier the same day,
by a relapse or by an earlier bite: a clone made in `reproduction()` is not visible until the day's
bites are over, while the host is already tagged infected, so the bite goes to the infected-pair
branch, which skips a host with no visible transmissible clone. That is 0-1% of infectious bites in
untreated runs and 2.8-5.5% with case management (`BLOODCURE_COVERAGE` 0.8) in the test's settings.
The same invisibility lets same-day bites take a host past `MAX_MOI`. D51 fixes both.
**How to verify.** `tests/testthat/test-model-invariants.R`: at `P_TRANS_VH = 1` with no host ever
full, each host's infectious bites a day are its inoculations that day, with case management too
(D51).
---

## D46 — Local vector control does not scale imports · 2026-11-07 · Accepted
**Context.** The Poisson import stream (`IMPORT_RATE_PER_TICK`) was scaled by the seasonal factor
of biting, by the vector-control factor 1 − `INTERVENTION_EFFICACY` inside the intervention window,
and by `IMPORT_RATE_BOOST` from the trigger; the fixed count (`N_IMPORTS_PER_TICK`) and the pulse
were not scaled at all, and none of this was documented (review m4, question Q6). Imports are
infections acquired elsewhere, which local vector control does not reach.
**Decision.** The vector-control factor no longer applies to imports. The seasonal factor stays on
the Poisson stream, and `IMPORT_RATE_BOOST` from `RESURGENCE_TRIGGER_TICK` and the import window
are the border measures. The fixed count and the pulse stay unscaled, as values the user sets
directly. MODEL_DESCRIPTION 2.1 and 8 list each route with its modifiers. `vivax_params()`'s
warning for an import rate above `N_HOSTS` follows the same factors. The author chose this (Q6,
2026-11-05).
**Rationale.** Bed nets and spraying reduce local biting, not the infections travellers bring back.
Travel, or transmission where travellers were infected, can follow the same season as local biting,
so the seasonal factor is kept.
**Alternatives rejected.** *Keep the factor* -- it ties a border measure to local vector control,
so an arm with vector control and open borders could not be expressed. *Scale the fixed count and
the pulse as well* -- they are set directly, and a user who wants them to vary can vary them.
**Supersedes.** D19's statement that `INTERVENTION_*` scales the importation rate (its separation
of case management from vector control stands), and the same statement where D33, D37 and D38 make
it in passing: full vector control as a closed border, and the vector-control factor in the
import-rate warning.
**Consequences.** A run with an import-rate window that overlaps an intervention window with
`INTERVENTION_EFFICACY > 0` brings more imports while the two overlap, and changes from the first
such day; every other run is byte-identical. The edge-value test closes the border with
`IMPORT_RATE_BOOST = 0` and a seasonal amplitude above 1, not with full vector control.
**How to verify.** `tests/testthat/test-edge-values.R`: with full vector control from day 10, no
bite or uptake happens from day 10 while imports keep arriving. `tests/testthat/test-params.R`: an
import rate of 40 a day for 30 hosts warns even under vector control.
---

## D45 — A run started with slim -d is checked as R checks it, and takes either spelling of a number · 2026-11-07 · Accepted
**Context.** `validate_params()` refuses every value the model cannot run (D33), and the model's
`initialize()` repeated only some of those checks for runs started directly with `slim -d`, which
bypass R. The rest ran: `RADCURE_MODE = 'Person'` silently ran the per-hypnozoite mode, a
probability above 1 was always true, and 12 founder hosts in a 10-host run seeded 10. And `-d X=7`
makes an Eidos integer where `-d X=7.0` makes a float, which Eidos will not mix in `min()`,
`rbinom()`, `seqLen()`, `setSeed()`, subsetting or a tick range, so a legal spelling of ten
parameters stopped the run: `ALPHA_ACTIVATE`, `MU_HYPNOZOITE`, `P_RAD_SURVIVAL`, `N_HOSTS`,
`N_MOSQ`, `N_FOUNDER_LINEAGES`, `N_IMPORTS`, `N_IMPORTS_PER_TICK`, `SEED` and `SIM_DAYS` (review
m62).
**Decision.** `initialize()` makes the rest of `validate_params()`'s checks: each parameter's type
and that it is a single value, whole numbers within ±2^53 for counts and days, finite numbers, the
probability ranges, positivity and non-negativity, the founder bounds, `MOSQ_LIFESPAN`,
`MEAN_BLOOD_STAGE_DURATION`, `MU_HYPNOZOITE`, `DELTA_LATENCY`, `RECOMB_RATE`, `SEASONAL_PHASE`,
`RADCURE_MODE`, a coverage against a closed treatment window, the mean oocyst count per uptake, an
overflowing biting rate, and `CHROMOSOMES` unreadable, listed twice or as a high-low range. Its
messages give a float to 15 significant digits and quote a string. The model converts a parameter
to the type it needs at each place Eidos is strict, so either spelling of a number gives the same
run. Not repeated: a misspelt name, which SLiM defines as a constant the model never reads, so the
model cannot see it; R's limits on what text can carry (UTF-8, control characters, values below
1e-300); and `SIM_DAYS` outside 2 to 1e9, INF, NAN or a vector, which SLiM refuses itself while it
reads the tick ranges, before `initialize()` runs, in messages that do not name `SIM_DAYS`.
**Rationale.** A direct run is the one place where a typo cannot be caught in R, and every one of
these ran to exit 0 meaning something else.
**Alternatives rejected.** *Refuse the other spelling* -- `-d TREATMENT_DELAY=7` is already
accepted on purpose (D32), and a rule that a float parameter must be written `7.0` is one users
would trip over. *Convert the constants once in `initialize()`* -- SLiM 5.2 cannot remove a
constant set with `-d`, so it cannot be redefined in another type.
**Consequences.** Every run started from R is byte-identical: R writes each parameter in its own
type and passes every check. In a direct run, the values listed are refused with the parameter
named. Across 318 pairs of direct runs on three configurations, each parameter in each spelling
gives byte-identical `events.tsv` and `state.tsv`.
**How to verify.** `tests/testthat/test-slim.R` refuses each new case with its parameter named, and
runs a 5-day configuration with latency, an MDA round in the per-hypnozoite mode and a pulse of
imports in both spellings, with identical outputs.
---

## D44 — MDA reaches subpatent clones and clones that became patent that day · 2026-11-07 · Accepted
**Context.** MDA goes through `applyTreatment()`, the case-management course (D20, D28), so it
inherited two filters meant for case courses: a course skips subpatent clones (a holder is below
patency, and re-rolling it emitted a second episode end), and it skips clones that became patent on
its own day, which cannot yet have caused symptoms (D24). So every suppressed clone survived a
round and later recrudesced (at 3a18dc4, the round in this entry's seed-11 run met 11 holders and
cleared none), although clearing undetected carriage is what MDA is for. D25 had left open what MDA
does to a holder (review m1).
**Decision.** An MDA course bypasses both filters. Every clone in a reached host, holders and
clones that became patent earlier that day included, is cleared with probability 1 −
`P_BLOODCURE_SURVIVAL`; the day's imports arrive after the round. A cleared holder is logged as
`mda_subpatent_clear`, which ends no episode (its patent episode ended at `bloodcure_subpatent`),
and its pending recrudescence goes with it; a holder that survives keeps its schedule and logs
nothing, since it is not suppressed again. A clone cleared on its first patent day is logged as
`mda_bloodstage_clear`, like any MDA clearance. Case courses keep both filters. Coverage keeps its
per-host draw among hosts not empty at the round (D28): the draw does not depend on infection, so
the coverage is that of the whole population, and `mda_bloodcure_given` and `mda_radcure_given`
count courses to hosts not empty at the round. A round on day 1 is refused when a founder host
carries a copy of another founder's genome: the copy is seeded on day 2, and a day-1 round cleared
its template and stopped the run. The author chose this (Q9, 2026-11-05).
**Rationale.** MDA is given to everyone, symptoms or not, so neither filter's reason applies to it:
the D24 floor stands for the delay from patency to symptoms and care, and the holder filter for a
course given for some other case. Drugs taken on the day a clone reaches the blood act on it too.
**Alternatives rejected.** *Keep MDA blind to holders and document it* -- a round would leave the
very reservoir it exists to clear. *Draw a surviving holder's suppression again* -- it would reset
the holder's schedule and log a second episode end.
**Supersedes.** D24's floor for MDA courses (case courses keep it), D28's statement that MDA
inherits that floor, and D25's open question of what MDA does to a holder.
**Consequences.** A run with an MDA round changes from the round's day whenever a reached host
carries a holder or a clone that became patent that day, since those clones now take a draw and
shift the random stream. Episodes that end on their first patent day exist again, on the MDA day
only. Runs without MDA are byte-identical. In a 50-host run in which every failed course suppresses
(seed 11, `P_BLOODCURE_SURVIVAL = 0.5`), the round met 11 holders and cleared 4; with every course
curing, no blood-stage clone but the day's imports is left after a full round, where the 4 clones
that became patent that day used to survive it.
**How to verify.** `tests/testthat/test-treatment.R` checks that holders at the round are either
cleared (`mda_subpatent_clear`, and never recrudesce) or left on their schedule (not suppressed
again, with no event that day), that the round does not clear a clone it has just suppressed, and
that nothing clears twice; and that with `P_BLOODCURE_SURVIVAL = 0` a round on day 2 clears every
founder copy seeded that morning, leaving no blood-stage clone, and `plot_host_trajectory()` draws
such a one-day episode.
---

## D43 — A founder-free run whose imports never take ends cleanly, as never started · 2026-11-06 · Accepted
**Context.** A founder-free run starts only from its imports. When every import failed to establish
a parasite, by chance (a Poisson rate that draws none, imports that take in neither compartment) or
because the settings let none take (`P_PRIMARY = 0` with a hypnozoite batch mean of 0, an import
rate of 0 throughout its window), the model's liveness check stopped the run with "this
configuration cannot bootstrap": `run_sim()` raised an error and the `run_grid()` cell was lost.
The check dates from a bug in which imports into an empty population were never seeded, so such
runs had finished silently with an empty event log. It fired only once the last import day had
passed, so a run whose import path lasted to its last day (the default `IMPORT_END`), or that was
shorter than 10 days, still ran to the end without a word. D33 left this case open.
**Decision.** Such a run ends cleanly on the first tick after its last import day, but not before
tick 10, or on its last day if that comes first: it writes its tree sequence, its `params.tsv`
ends with a `NEVER_STARTED` row, `run_sim()` says so, `load_run()` gives `never_started = TRUE`,
and `run_grid()` adds a `never_started` column and counts these cells in a message of their own.
Unlike extinction (D37), it is decided on the last day whatever the import schedule, so every
founder-free run that never had a parasite is marked. Settings that let no import take are not
refused. A founder-free run with no import path inside it at all is still refused before SLiM
starts (R and the model's `initialize()`). The author chose this (2026-11-06).
**Rationale.** A failure of importation is a result, as fade-out is: how often imports fail to seed
an outbreak is itself worth counting, and a sweep should not lose the cell. A regression of the
old bug would now show as never-started cells, which the grid's message and column make visible
rather than silent.
**Alternatives rejected.** *Keep the error* -- it loses sweep cells to chance. *Record it as an
extinction* -- a run that never started and one that died out mean different things. *Refuse a
founder-free set with `P_PRIMARY = 0` and a batch mean of 0* -- like a closed border (D33), it is a
legitimate null arm of a grid, and never started is its honest outcome.
**Consequences.** Only founder-free runs whose imports all fail change. Those the liveness check
stopped end with exit status 0 instead of an error; those that ran to their last day gain the
`NEVER_STARTED` row, with their events and state unchanged. Every other run is byte-identical. A
run that never started holds no parasite (a tree sequence with no nodes, and at most no-take rows
in its events). Every reader takes such a run, as it does an extinct one: its tables have their
usual columns and no rows, its counts are 0, and a quantity that needs a parasite is `NA`
(`tests/testthat/test-no-parasite.R`). `?run_grid` says to keep the `ok & !never_started` cells for
event and genetic summaries, while the cells that never started still count as outcomes.
**How to verify.** `tests/testthat/test-run.R` runs a founder-free set whose import rate draws no
import: `run_sim()` gives the message, `never_started` is `TRUE`, `extinct_tick` is `NA`,
`load_params()` has no `NEVER_STARTED`, and `run_grid()` counts both cells. A second test checks
the tick such a run ends on (10 for a day-1 pulse that cannot take, 16 for one on day 15, and the
last day for a window to the end, a pulse on the last day and an 8-day run), and that a
founder-free run whose import took as hypnozoites only, and died, went extinct.

---

## D42 — MEAN_RECRUD_DELAY = 10.5 d is a literature assumption, pending parameterisation · 2026-11-06 · Accepted
**Context.** D26 chose 10.5 d as "the reciprocal of Taylor's posterior median", replacing an
uncited 21 d. Taylor et al. (2019) fix the recrudescence rate with a strongly informative
Gamma(100, 1000) prior (mean 0.1 a day) for identifiability, and their posterior median, 0.095 a
day, barely moves from it (checked in the paper and their Stan code). So 10.5 d is close to their
prior, not an estimate from data.
**Decision.** Keep 10.5 d as the default and call it what it is: a literature assumption, pending
parameterisation. MODEL_DESCRIPTION 9a and 10 and the model comment say so, and section 14 lists
it for sensitivity analysis (S11). The author agreed (2026-11-06).
**Supersedes.** D26's description of 10.5 d as a fitted value; its choice of a geometric delay,
with selectable shapes, stands.
**Consequences.** No run changes. Analyses that report recrudescence should vary the delay.
**How to verify.** MODEL_DESCRIPTION 9a, 10 and 14 (S11), and the comment on `MEAN_RECRUD_DELAY` in
the model's `initialize()`.

---

## D41 — Meiosis I keeps sister chromatids together · 2026-11-05 · Accepted
**Context.** D34 dealt each chromosome's four chromatids to the four products at random,
independently per chromosome, as the simulator of Wong et al. (2018) does. A review of that change
found that this drops meiosis I: the two chromatids attached to one centromere go to the same pole,
so the two products of a pole carry the same parent at every chromosome's centromere. Under the
random deal any two products of one oocyst share about a third of the genome (sd about 0.08,
unimodal); under meiosis I the two pairs of meiosis II sisters share about 0.60 and the other four
pairs about 0.20 (overall sd about 0.20; about 0.06 within each group). IBD and relatedness
methods are scored against this ground truth, and the narrow distribution would flatter them. The
author chose to model meiosis I (2026-11-05).
**Decision.** For each chromosome, meiosis I sends the two chromatids that carry the same parent at
its centromere to one pole, products 0 and 1 or products 2 and 3, drawn independently for each
chromosome; meiosis II gives one of them to each product of that pole, in random order. Chiasmata
and crossovers are as in D34. The centromeres are the centres of the features PVP01_CEN01 to
PVP01_CEN14 of the PvP01 annotation (Sanger/GeneDB, released with Auburn et al. 2016), held as
fractions of each chromosome (`PVP01_CEN_FRAC`) so that `GENOME_MODE = "core"` keeps them. A
tetrad's work is linear in its chiasmata: they are de-duplicated without keeping their order
(they are sorted next), and each chiasma's two chromatids are appended to two vectors. A tetrad
takes about 24 ms at `RECOMB_RATE = 1e-3`, against about 320-440 ms before; the products and the
random stream are unchanged.
**Rationale.** This is canonical meiosis, and we know of no evidence that Plasmodium departs from
it. The split barely depends on where the centromeres sit: midpoints give 0.614 and 0.193, the
annotated centromeres 0.603 and 0.198, and centromeres anywhere in 20-80% of each chromosome about
0.605 and 0.197. No obligate chiasma is imposed (D34); one would move sisters to about 0.58 and
the others to about 0.21, about twice the effect of the centromere positions.
**Supersedes.** D34's assignment of chromatids to products, and the figures that follow from it
(unimodal around 1/3, sd about 0.08, "this is Plasmodium meiosis"); and D34's expectation that a
sibling pair's IBD state is uncorrelated across junctions. Pooled over pairs it is about 0.09
under meiosis I, since a pair that are sisters on one chromosome are sisters on all; a product's
own parent stays uncorrelated across junctions.
**Alternatives rejected.** *Keep the random deal* -- simpler and shared with Wong et al. (2018),
but it makes every same-oocyst pair look alike. *Pair products without centromeres* (a deal fixed
at each chromosome's start) -- gives sisters 0.50 and the others 0.25, so the centromere positions
matter.
**Consequences.** Every run with an outcrossed oocyst changes its tree sequence. On the full
14-chromosome map at the default rate `events.tsv` and `state.tsv` are unchanged: the genome is
neutral, every product is a fresh lineage, and the new deal draws as many random numbers as the
old (Eidos's sample() draws n - 1 numbers to order n items, so 3 either way). With a `CHROMOSOMES`
subset or a small `RECOMB_RATE` (every one of six seeds differed at 1e-8), an unrecombined product
is a copy that keeps its parent's lineage, and the deal decides which product that is, so those
runs' events change too. Each product's
crossovers and assortment are unchanged (17.9 M; a junction switch with probability 1/2). Two
products of one oocyst share a third of the genome on average for unrelated parents, sisters about
0.60 and the others about 0.20 (5-95% about 0.12-0.67 over all pairs), and a pair's sharing is
correlated across chromosomes. The README figures (karyotype, representative run) and the observed
values in MODEL_DESCRIPTION section 13 are regenerated with the release (M33).
**How to verify.** `tests/testthat/test-model-invariants.R` runs the model's own `seedTetrad()` for
400 tetrads with a probe at every centromere: 2:2 segregation at every probe, two sister pairs per
oocyst carrying the same parent at all 14 centromeres, sister and other pairs' sharing against
1/3 + 2/3 e^(-3rt) and 1/3 - 1/3 e^(-3rt) at map distance rt from the centromere, each product's
crossovers against the map (mean and variance), and a junction switch probability of 1/2.

---

## D40 — A relapse into a patent clone of its own lineage is counted, and flagged · 2026-11-05 · Accepted
**Context.** At the defaults about 80% of relapses start while a patent clone of the same lineage
is already in the host (70-82% in review M3's larger runs). Over four 3-year default runs that
clone was the primary from the same bite for 45% of them, an earlier relapse of the same batch for
46%, another batch's relapse for 8% and another inoculation of the lineage for 1%. Each such
relapse is logged as a new episode and counted by `relapse_fraction()` and the other metrics,
although a real observer could not see it as a recurrence. The model's own observation layer
draws its detection like any other relapse's, so it can be reported and treated (at the default
`P_DETECT_RELAPSE = 1` every one is reported), while a study would see only the infection already
there (review of M3).
**Decision.** The model is kept: a new genome per activation, counted as other relapses are. Each
`hypnozoite_activation` row gains `into_patent_lineage`: 1 if a patent clone of the same lineage,
born on an earlier day, was still in the host, 0 if not (-1 on other rows). The activation runs in
`reproduction()`, which sees no clone born in it, so a relapse that starts on the same day as
another clone of its lineage in the host (another batch waking, a recrudescence or a bite) is not
flagged for that clone, and two batches of one lineage that wake together, with no earlier clone of
it in the host, are both 0; this is rare at the defaults, and is documented rather than resolved
here (review of M3). In a later change, `classify_cases()`, `relapse_fraction()` and
`force_of_infection()` gain `exclude_masked = FALSE`, and the 0.55-0.85 relapse band is compared
with the acquisitions reading with these relapses left out.
**Alternatives rejected.** Merging such an activation into the ongoing clone. It would change
every run, and the flag gives each analysis both readings.
**Consequences.** `events.tsv` gains a last column; its other columns, and every other output, are
byte-identical.
**How to verify.** `tests/testthat/test-run.R` rebuilds each host's patent clones from
`events.tsv` in a run with treatment failures, recrudescence and an MDA round, and the flag
matches on every activation. It also matched on all 1,331 activations of eight other runs
(defaults, treatment failures, MDA, temperate latency, imports with `MAX_MOI = 3`). At the
defaults, 81% of 12,219 activations over eight 3-year runs are flagged (76-95% by run). On a single
chromosome the share is higher (86%), because an unrecombined product there keeps its parent's
lineage.

---

## D39 — IMPORT_TICK = 1 is an ordinary pulse · 2026-11-05 · Accepted
**Context.** `IMPORT_TICK = 1` created every host empty, silently dropping `N_FOUNDER_HOSTS` and
`N_FOUNDER_LINEAGES` whether or not `N_IMPORTS > 0`: a pulse on the first day turned the run into a
naive population (review M5). With `N_IMPORTS = 0` the run then stopped at the liveness check with
advice about `N_FOUNDER_HOSTS`. Since D17 a naive population is expressible directly.
**Decision.** The special case is removed. A pulse on day 1 lands like a pulse on any other day,
alongside the founders, with one difference (review of M5): a founder host that carries a copy of
another founder's genome (`N_FOUNDER_LINEAGES < N_FOUNDER_HOSTS`) receives the copy only on day 2,
since `addCloned()` needs `reproduction()`, which first runs then. On day 1 such a host is still
empty, so the `MAX_MOI` check does not refuse it a day-1 import, and from day 2 it carries both. With
20 hosts, 12 of them founders sharing 2 genomes, and a day-1 pulse of 20 at `MAX_MOI = 1`, only the
2 hosts of the founding genomes refused the import, and all 10 hosts carrying copies held two
lineages on day 2 (seeds 1-3). On day 2 the copy is taken from the founding genome in its host that
is not an import, since a day-1 import can sit beside it. A naive population, seeded only by imports,
is `N_FOUNDER_HOSTS = 0` with an import path.
**Consequences.** Runs with `IMPORT_TICK = 1` keep their founders, and so change. No analysis in
either repository uses `IMPORT_TICK = 1`, and no saved `params.tsv` has it. Other runs are
byte-identical. Picking the template by origin rather than as the host's first genome changes no
run: an import is appended after the founder, so the first genome was the founder already.
**How to verify.** `tests/testthat/test-run.R`: with `IMPORT_TICK = 1` the founders are seeded and
the imports arrive on tick 1.

---

## D38 — Imports go to uniformly random hosts, with the MAX_MOI check a bite gets · 2026-11-05 · Accepted; partly superseded by D46 (vector control does not scale imports)
**Context.** The importation block filled its quota from hosts without a blood-stage infection
first, using infected hosts only when those ran out, and checked no capacity (review M6). An
import therefore almost never arrived in a host that already carried a clone, which made import
detection and classification easier than the documented random placement, and the docs said
random hosts.
**Decision.** Imports go to hosts sampled uniformly at random, at most one per host a day as
before, so a day's imports beyond `N_HOSTS` are dropped. A host already carrying `MAX_MOI`
blood-stage lineages takes none, as a local bite does, and the refusal is logged as
`vh_inoculation_import_no_take`. `vivax_params()` warns when `N_IMPORTS` or `N_IMPORTS_PER_TICK`
exceeds `N_HOSTS`.

The review of M6 found the two gates unequal. A bite was checked only on a host tagged infected,
which needs an infectious clone, so at `SUBPATENT_INFECTIVITY = 0` a host whose clones were all
subpatent took bites past the cap while it refused imports; and `MAX_MOI = 0`, which
`vivax_params()` accepted, refused every import, an empty host's included, while bites still
infected empty hosts. A bite on a host tagged susceptible is now checked the same way, before its
`P_TRANS_VH` draw, so the gate is the same for every new inoculation: a host already carrying
`MAX_MOI` blood-stage lineages, patent or subpatent, takes none. `MAX_MOI` must be at least 1, in
`validate_params()` and in the model's `initialize()`; 1 means no superinfection, which is what 0
gave a bite. The warning also covers the Poisson rate, whose draws are cut at `N_HOSTS` a day too:
`IMPORT_RATE_PER_TICK` at its peak over the window's days inside the run, worked out as for the
biting rate (seasonality floored at 0, the intervention factor inside its window,
`IMPORT_RATE_BOOST` from `RESURGENCE_TRIGGER_TICK`), warns when, with `N_IMPORTS_PER_TICK`, a day
draws more than `N_HOSTS` imports with probability above 0.01, and the warning names the day, the
rate and the probability. A rate of 50 a day into 30 hosts loses about 20 imports on almost every
day (P = 0.998), and the old comment's "only by chance" did not hold for it. When the pulse and the
window's count each exceed `N_HOSTS`, on different days, the warning names both.
**Alternatives rejected.** *Placing every import, with replacement*, as review M6 proposed. At
D33's cap of 1e6 imports a day each host would take about 25,000 imports a tick at the default
`N_HOSTS = 40`, and a run started directly with `-d` at that cap never finished. *Saying in the
docs that a bite is gated only on a host with an infectious clone* -- two gates for what is
documented as one rule about parasite burden.
**Consequences.** A run changes once imports arrive while some host is infected: imports now
land in infected hosts in proportion to prevalence, so mixed infections with an import, and
treatment of imported cases, are more common. The random placement leaves runs without imports,
and a founder-free pulse into a wholly susceptible population, byte-identical; at `MAX_MOI = 0`,
now refused, that pulse was refused outright. The bite gate changes only runs with `SUBPATENT_INFECTIVITY = 0` in which a
host whose clones are all subpatent reaches `MAX_MOI`: whenever `SUBPATENT_INFECTIVITY > 0` a host
tagged susceptible carries no blood-stage clone, so the check never fires. The rate warning is
R-side only.
**How to verify.** `tests/testthat/test-run.R` imports 100 infections into 200 hosts, half of
them founders: about half the imports land in founder hosts (none did before), and with
`MAX_MOI = 1` none lands in a founder host that still carries its clone when the pulse arrives (a
founder host whose clone cleared first rightly takes one), and the refusals are logged. A pulse of
50 into 30 hosts brings exactly 30, with the warning. With `SUBPATENT_INFECTIVITY = 0`,
`MAX_MOI = 1` and every course failing into subpatency, no bite reaches a host after its clone is
suppressed (59 to 87 did over seeds 1-4 before). `tests/testthat/test-params.R` refuses
`MAX_MOI = 0`, warns for a rate of 50 into 30 hosts and follows the boost, season and intervention
timeline, and names both counts; `tests/testthat/test-slim.R` refuses `MAX_MOI = 0` in a run
started directly.

---

## D37 — A run that goes extinct records the tick and says so · 2026-11-05 · Accepted; partly superseded by D46 (vector control does not scale imports)
**Context.** A run in which every parasite and hypnozoite was gone, with no importation to come,
ended on that tick without telling R (review M7). Only run.log carried the EXTINCT line;
`run_sim()` returned as for a full run, and `run_grid()` reported the cell as ok. Metrics that take
the run length from the last tick then used a truncated window.
**Decision.** The extinction branch appends `EXTINCT_TICK` to `params.tsv`. `load_run()` reads it
into the run's `extinct_tick` (`NA` for a run that did not go extinct), `print()` shows it,
`run_sim()` gives a message, and `run_grid()` returns an `extinct_tick` column with one message
for the grid. `state.tsv` is not padded to `SIM_DAYS`. The metrics' windows move to `SIM_DAYS`, with
ticks after extinction counted as uninfected, in the R-side half of M7, still to come.

The review of M7 found four gaps, closed as follows. A run that dies out on its last day records
the tick too, so that `NA` means only that the run did not go extinct, which is what a grid counting
extinctions needs; the messages say "went extinct", not "ended early". `EXTINCT_TICK` is an outcome
rather than a parameter, so `load_params()` leaves it out of the list it returns, and every run's
list has the same names. The documented condition names what the check counts: no parasite in any
host or mosquito. And an import set for a day after `SIM_DAYS`, or a count left at its default day,
counted as still to come, so a run whose parasites had died out ran on, empty, to `SIM_DAYS` with
`extinct_tick` `NA`, while `has_import_path()` and the model's bootstrap check counted a pulse on
day 0 or earlier, which never comes, as an import path. One constant in `initialize()`,
`LAST_IMPORT_DAY`, the last day of the run on which an import can arrive (0 if none), now serves the
bootstrap check and the extinction check alike, and `has_import_path()` and the import warning in
`vivax_params()` share one R definition of the same days: a pulse counts when `N_IMPORTS > 0` and
`IMPORT_TICK` is a day from 1 to `SIM_DAYS`, and a window when it has a count or a rate and a day
from `max(1, IMPORT_START)` to `min(IMPORT_END - 1, SIM_DAYS)`. A window inside the run still counts
to its last day while its rate is 0 (a border closed by `INTERVENTION_EFFICACY = 1` or
`IMPORT_RATE_BOOST = 0`), so such a run goes extinct only once that day has passed; this is
documented rather than worked out from the rate.
**Consequences.** `params.tsv` of an extinct run gains one row. A run with an import set after
`SIM_DAYS` (a count whose day is left at the default included), or a window that closes before it
opens, goes extinct when its parasites die out, once its import days inside the run, if any, have
passed, instead of running on, empty, to `SIM_DAYS`; it loses only empty `state.tsv` rows. A
founder-free set whose only import falls before day 1 (a pulse on day 0 or earlier, or a window
that ends by day 1) is refused before SLiM starts instead of stopping at the liveness check on
tick 10. Every other run is
byte-identical.
**How to verify.** `tests/testthat/test-run.R` runs a set that dies out: the message appears,
`extinct_tick` is below `SIM_DAYS` and equals the last tick of `state.tsv`, `load_params()` has the
same names as for a full run, and `run_grid()` reports it; a full run has `extinct_tick` `NA`. The
same set with a pulse after `SIM_DAYS`, a window after it, or a count left at the default day goes
extinct on the same tick, and with `SIM_DAYS = 100000` the message gives the day as 100000.
`tests/testthat/test-params.R` refuses a founder-free set whose pulse falls on day 0 or whose window
ends on day 1, and `tests/testthat/test-slim.R` checks that the model refuses it when run directly.

---

## D36 — A bite that seeds only hypnozoites is no case · 2026-11-05 · Accepted
**Context.** `inoculateMosquitoToHuman()` drew detection once per bite and scheduled a course
whenever the bite established anything. With `P_PRIMARY < 1` a bite can seed only a hypnozoite
batch; it made no blood stage, yet a detected one got a blood-stage course and, with
`RADCURE_COVERAGE`, radical cure that removed the batch it had just seeded (review M2). Imports
already skipped this.
**Decision.** A course is scheduled only when the bite created a blood-stage primary (D20's
"decided where a case arises"). The detection draw is unchanged, so the random stream upstream is
the same.
**Consequences.** Runs with `P_PRIMARY < 1` change wherever a detected bite seeds only hypnozoites
while the case-management window is open, as it is for the whole run by default, even at zero
coverage: `decideTreatment()` draws the coverage (unless a course is already pending), at 0 too, so
the call such a bite no longer makes shifts the random stream. With coverage, such a bite also gets
no course: no phantom `bloodcure_given` or `radcure_given` row, no clearance of the host's other
clones, and no radical cure of the batch it has just seeded. Measured at `P_PRIMARY = 0.5` and zero
coverage over three 120-day runs (review of M2): each run differs from the model before, with no
course in any of them, and each is byte-identical to it with the window moved past the run
(`TREATMENT_START = 200`, `TREATMENT_END = 300`) or with `P_DETECT_PRIMARY = 0`. With
`P_PRIMARY = 1` (the default) every bite that establishes anything makes a primary, so those runs
are byte-identical.
**How to verify.** `tests/testthat/test-treatment.R`: at `P_PRIMARY = 0`, full coverage and
undetected relapses, no `bloodcure_given` or `radcure_given` row appears (there were 34 before).

---

## D35 — Latent hypnozoites die at MU_HYPNOZOITE · 2026-11-05 · Accepted
**Context.** With `K_LATENCY > 0` the latency branch only drew the batch's move to the next stage,
so hypnozoites could not die while latent (review M4). White et al. (2014), whose batch model the
reservoir follows, apply death at rate mu in every latent compartment. At `K_LATENCY = 2` and
`DELTA_LATENCY = 1/100` a hypnozoite relapsed with probability 0.722 rather than 0.480, about 1.5
times the relapses per bite, all of the extra ones late.
**Decision.** Each day a latent batch first loses Binomial(`batch_count`, `MU_HYPNOZOITE`)
hypnozoites, logged as `hypnozoite_clearance` like any other death; then, if any are left, the batch
moves to the next stage with probability `DELTA_LATENCY`. The batch still moves through latency as
one unit, which matches White et al., who also track latency per infection.
**Rationale.** It is the model the documentation cites. P(relapse) becomes
[delta(1 - mu)/(1 - (1 - delta)(1 - mu))]^k x alpha/(alpha + mu), close to
(delta/(delta + mu))^k x alpha/(alpha + mu).
**Supersedes.** D31's "Until hypnozoites can also die during latency, a temperate batch
(`K_LATENCY > 0`) leaves latency only at rate `DELTA_LATENCY`", given there as the reason that
`validate_params()` requires `DELTA_LATENCY` in (0, 1] when `K_LATENCY > 0`: a latent batch now
also empties by death. The requirement stays, because at `DELTA_LATENCY = 0` a temperate batch
would never leave latency, and so could never relapse.
**Consequences.** Runs with `K_LATENCY > 0` change: fewer relapses and a smaller reservoir.
The latent deaths leave tropical runs (`K_LATENCY = 0`, the default) byte-identical. Since
survival favours the batches that move through quickly, a hypnozoite that comes through latency alive has been latent
k/(1 - (1 - delta)(1 - mu)) days on average, 163 at k = 2 with the default `DELTA_LATENCY` and
`MU_HYPNOZOITE`, not the k/delta = 200 of a batch without deaths (MODEL_DESCRIPTION §5.1; review of
M4). The `hypnozoite_clearance` rows, latent deaths and the others, and the `human_clearance` rows
carry their lineage's `is_imported`, as D20 intends and MODEL_DESCRIPTION §16.1 defines; all three
wrote 0 (review of M4). In runs with imports only that column of those rows changes; the random
stream is untouched.
**How to verify.** `tests/testthat/test-run.R` imports 300 infections at `K_LATENCY = 2` with no
onward transmission, with `SIM_DAYS = 1500`, by which every batch has resolved: at the test's seed
the run ends by extinction on tick 1142. The share of resolved hypnozoites that activated must be
within 0.05 of 0.480 (0.47 at the test's seed); it was 0.73 before. Two runs seeded only by imports,
one whose batches stay latent throughout and one without latency, carry `is_imported = 1` on every
`hypnozoite_clearance` and `human_clearance` row.

---

## D34 — Outcrossed oocysts undergo per-chromosome four-strand meiosis · 2026-11-05 · Accepted; partly superseded by D41 (how chromatids reach the products)
**Context.** `seedTetrad()` built the four products of an outcrossed oocyst as two reciprocal
recombinants from one breakpoint draw plus two whole-genome copies of the parents. That applies a
single chromosome's tetrad to the whole genome. Half the products saw no crossover and no
independent assortment, the realised map averaged half the documented 17.9 M per product, a third
of same-oocyst pairs shared nothing, and the parental copies carried their parents' lineage onward
(review M1).
**Decision.** Per-chromosome four-strand meiosis without chromatid interference, as checked against
Walliker et al. (2005) and Camponovo et al. (2023). The zygote's two genomes replicate into four
chromatids per chromosome. Chiasmata fall at the union of two gamete-rate breakpoint draws, a
Poisson set with twice the per-product map length, within each chromosome; each joins a random
chromatid of each parent at that locus. The four chromatids of each chromosome then go to the four
products at random, independently across chromosomes. No obligate chiasma is imposed: neither
source claims one, and its effect on the map is negligible (on sibling sharing under meiosis I,
see D41). A product with no switch at all is a copy of its
starting parent and keeps that parent's lineage; every other product is a fresh lineage. That
happens with probability 2^-(n-1) e^-d for n chromosomes of total map d: about 2e-12 on the full
map, but often on a few chromosomes (about half of the products of chromosome 1 alone). Selfed oocysts are unchanged. The chromosome starts on D2's concatenated map are
recorded as `CHROM_STARTS`.
**Rationale.** This is Plasmodium meiosis. Every product carries the full map and assorts its 14
chromosomes independently, two products carry each parent's allele at every locus, and any two
products of one oocyst share about a third of the genome (sd about 0.08, 5-95% about 0.21-0.47),
none sharing nothing. The ground truth that IBD and relatedness methods are scored against is then
the biology rather than an artefact of the tetrad.
**Supersedes.** D2's mechanism of assortment, the r = 0.5 junction breakpoints passed through the
earlier `seedTetrad()` (the junction draws are now discarded, and the tetrad assorts the
chromosomes itself), and D2's validation sentence that the same-oocyst distribution "tightens
toward the center under 14 chromosomes"; D2's concatenated map stands.
**Alternatives rejected.** *Two independent reciprocal pairs* -- this restores the map length and
removes the parental copies, but two of the six same-oocyst pairs still share nothing by
construction. *An obligate chiasma per bivalent* -- no source supports it for P. vivax, and it
would barely change the map.
**Consequences.** Every run with an outcrossed oocyst changes, through the products and the random
stream; runs in which every oocyst is selfed are unchanged. Same-oocyst IBD is unimodal around 1/3
instead of split between 0 and about 1/2, and inoculated genotypes are no longer half parental
copies, so counts of distinct lineages can rise. The README karyotype figure and the observed
values in MODEL_DESCRIPTION section 13 come from the earlier tetrad and are regenerated with the
release (M33).
**How to verify.** `tests/testthat/test-model-invariants.R` runs the model's own `seedTetrad()` on
the PvP01 map for 400 tetrads: two products carry each parent's allele at every probe locus, none
is a whole-genome copy of a parent, no same-oocyst pair shares nothing, and the mean pair share is
within 0.02 of 1/3. `inst/validation/audit_arg_junctions.R` shows each product's parent uncorrelated
across junctions and correlated within chromosomes.

---

## D33 — Legitimate edge values run to the end; values the model cannot run are refused first · 2026-11-05 · Accepted; partly superseded by D43 (founder-free runs whose imports never take) and D46 (vector control does not scale imports)
**Context.** A search of the model's random draws, empty vectors, parameter edges and arithmetic
found about seventeen ways for a run to stop part-way with an Eidos error, or never finish, on
values `validate_params()` accepted. Some of those values are legitimate and the model mishandled
them: `SUBPATENT_INFECTIVITY = 0`, the competing-risks clamp, an import rate that reaches 0
(border closure, full border control, a seasonal amplitude above 1), the `"poisson"` delay at its
minimum, and the small `K_OOCYSTS` that `uptake_calibrate()` can return. The rest are values the
model cannot run, or runs as something else: probabilities outside [0, 1], negative boosts and
counts, non-whole counts and scheduling days, values of the wrong kind, `Inf`, and sizes far
beyond anything realistic.
**Decision.** The model takes every legitimate edge value at its natural limit: no imports at a rate
of 0, a one-day delay at the minimum, a death probability capped at 1, no infectiousness for a host
without a blood stage, a zero count from a gamma draw of 0, and an exact draw, by inversion, of the
oocyst count of an uptake whose 10,000 gamma-Poisson tries find none. `validate_params()` refuses
everything else before SLiM starts, with a message that names the parameter. Each value is a
number, `TRUE`/`FALSE` or a string, as its default is; numbers are plain and finite, and 0 or at
least 1e-300 in size. Counts, `SIM_DAYS`, the scheduling days (`MDA_TICK`, `IMPORT_TICK`,
`RESURGENCE_TRIGGER_TICK`, the `*_START`/`*_END` bounds) and the fixed delays `EIP_DAYS` and
`GAMETOCYTE_DELAY` are whole; mean durations, rates and `SEASONAL_PHASE` need not be.
Probabilities and coverages, `ALPHA_ACTIVATE` and `INTERVENTION_EFFICACY` lie in [0, 1],
`MU_HYPNOZOITE` in (0, 1] and `DELTA_LATENCY` in (0, 1] when `K_LATENCY > 0`; `BITING_RATE` is
positive; boosts, import rates, `LAMBDA_*`, `SEASONAL_AMP` and counts are non-negative;
`MOSQ_LIFESPAN` is at least 1 (D29); `RECOMB_RATE` is at most 0.5. `CHROMOSOMES` names chromosomes
1-14 once each, and strings are valid UTF-8 without control characters (C0 or C1). A few sizes are
capped: a mean of 1000 oocysts per successful uptake, with `K_OOCYSTS` 0 or at least 0.001;
`LAMBDA_SPOROZOITES` at 1000; the hypnozoite batch mean at 10^6, and times its boost (if above 1)
once the trigger falls inside the run; the import rate at 10^6 a day; `TREATMENT_DELAY` and
`MEAN_RECRUD_DELAY` at 10^6 days; and `SIM_DAYS` at 10^9, SLiM's last tick. The oocyst mean,
`LAMBDA_SPOROZOITES`, hypnozoite batch and import caps and the `K_OOCYSTS` floor allow a relative
rounding of 10^-12 (10^-11 in the model, whose values arrive rounded to 15 digits); the delay and
`SIM_DAYS` caps are exact. A run of 99999 days or more must move the days that default to 99999
past its end, and cannot schedule anything on day 99999 itself. Spellings that mean the same value
are stored in the form SLiM is given: a named number or a 1 x 1 matrix, a value wrapped in `I()`,
a whole number with floating-point noise (`1.1 * 100`, or `0.1 * 3 - 0.3` for 0), a number or
factor for a string parameter (`expand.grid()` makes factors), text in another encoding (latin1,
read as Windows-1252; UTF-8 marked as bytes; or a non-UTF-8 native encoding R can convert) as
UTF-8, and spaces around the commas and dashes of `CHROMOSOMES`. A number with any other class (a
lubridate duration, a `bit64::integer64`, a `ts`) is refused, since it need not mean its bare
value. The model's `initialize()` repeats, for runs started directly with `-d`, the checks that
guard its own draws: the delays; `K_OOCYSTS`, 0 or at least 0.001 and finite, and `LAMBDA_*` at
most 1000; `INTERVENTION_EFFICACY` in [0, 1]; the biting, import and activation rate factors and
`SEASONAL_AMP`, non-negative and finite; the peak import rate and hypnozoite batch mean; and the
"never" days. The oocyst mean cap, which bounds the work of an uptake rather than guarding a draw,
is checked in R only. `run_sim()` keeps SLiM's output whole in `run.log` and, on Unix, passes
strings to SLiM as UTF-8 bytes, so a label survives a C locale.
**Rationale.** A sweep must not lose cells to engine errors, and a refusal that names the
parameter is better than a stack trace from tick 200. The caps sit far above any realistic
setting and well below where SLiM's draws fail (binomial sizes are 32-bit, Poisson draws wrap
above 2^32, the geometric delay saturates near 4e9 days), and they bound the work of a single
draw. The oocyst cap is on the mean count per uptake rather than on `LAMBDA_OOCYSTS`, because a
small `K_OOCYSTS` gives a tail whose mean is many times lambda. `K_OOCYSTS = 0.001` already comes
within 0.5% of the largest variance the zero-truncated negative binomial reaches at its mean (the
K -> 0 limit), for any mean up to 1000, while at a given mean the gamma-Poisson tries an uptake
takes grow about as 1/K below it (about 530 at a mean of 3 at K = 0.001), so a smaller K costs run
time for almost no dispersion; in a run started with `-d`, where the mean cap is not checked, the
work would grow without bound. An uptake whose 10,000 tries all give 0 used to settle for one
oocyst, which made the count less dispersed than `uptake_moments()` says; the exact draw that
replaces it is reached by 37% of uptakes where a positive count is 1 in 10,000 (`LAMBDA_OOCYSTS`
about 1e-4, at any K), by 0.005-0.1% at 0.001, and by fewer than one in a million above about
0.003. The rounding allowance lets a fitted or computed value that lands a few ulps above a cap
(`uptake_calibrate()` at a mean of 1000) run, and a refused value then prints visibly above the
cap. 99999 is the day the defaults use for "never", so a run that reaches it would trigger every
default intervention on that day; since `run_sim()` passes every parameter explicitly, a day set
to 99999 on purpose cannot be told from the default.
**Alternatives rejected.** *Clamp silently in the model* -- an out-of-range value would run as
something else. *Refuse the legitimate edge values too* -- border closure, a null arm without
subpatent infectivity and the minimum delay are real scenarios. *Cap the import rate at
`N_HOSTS` a day* -- the model runs such rates, and the cap only refused small populations.
*Refuse alternative spellings* -- `1.1 * 100` and the factors of an `expand.grid()` sweep ran
correctly, and refusing them would cost sweep cells for nothing. *Draw the oocyst count exactly
wherever a positive count is rare* -- every accepted `K_OOCYSTS = 0.001` set needs at least 110
tries per positive count, so the switch would change runs that never exhaust the tries.
**Consequences.** Every run that completed before and is still accepted gives byte-identical
events and state, except a run in which some uptake found no positive oocyst count in 10,000
tries (in practice only with `K_OOCYSTS > 0` and `LAMBDA_OOCYSTS` below a few thousandths; no saved
analysis set has `K_OOCYSTS > 0`): that count is now drawn exactly, and the run differs from there.
Many of the refused values used to complete while silently meaning something else, e.g.
`MOSQ_LIFESPAN = 0.5`, `N_FOUNDERS = 12.5`, `MDA_TICK = 50.5`, `CHROMOSOMES = "1,1"`,
`EIP_DAYS = -5`, `VERBOSE = "FALSE"` (which ran verbose), `P_TRANS_VH = "0.5"` (which Eidos
compared as text) and `SEASONAL_AMP = -0.5`. `MOSQ_LIFESPAN` and `SEASONAL_PHASE` become doubles,
since a mean and a day of the year need not be whole, so `params.tsv` writes them as `10.0` and
`0.0`. `uptake_moments()` keeps its precision for a tiny lambda or a huge `k`, where it returned a
mean of Inf or below 1, and where lambda/k overflows, and its variance is never negative.
`uptake_calibrate()` matches the target mean exactly, returns `k = 0` for a target at or below the
Poisson variance and otherwise a `k` of at least 0.001 (and warns when even `k = 0.001` falls more
than 1% short of the target variance), and refuses a target mean above 1000. Run time is not
bounded: `N_HOSTS` and `N_MOSQ` are uncapped, `SIM_DAYS` is capped only at SLiM's 10^9, and a high
`RECOMB_RATE` or a dispersed uptake within the caps can be slow. Not settled here: a founder-free
run whose imports fail to establish by chance still stops at the model's liveness check, and some
impossible import schedules are caught only there.
**How to verify.** `tests/testthat/test-edge-values.R` runs the legitimate edge values to the end
and checks that each run reached the state its guard handles; the minimum delays are in
`tests/testthat/test-treatment.R`, and the direct `-d` checks in `tests/testthat/test-slim.R`.
`tests/testthat/test-params.R` refuses each impossible value with a message that names the
parameter, and stores the alternative spellings in SLiM's form.
`tests/testthat/test-model-invariants.R` checks the exact draw's mean and variance against
`uptake_moments()`. `tests/testthat/test-run.R` runs a set at the edge of the rounding allowance on
the caps, and a non-ASCII label in a C locale, reading it back unchanged from `params.tsv`, with
`run.log` complete.
---

## D32 — Treatment and recrudescence delays follow the 1/D convention · 2026-11-05 · Accepted
**Context.** D29 made every duration parameter the mean of a geometric duration with daily
probability 1/D, but two delays were left out. A course was scheduled
`max(1, rgeom(1, 1/(1 + TREATMENT_DELAY)))` days after the case, and a suppressed clone re-expanded
`max(1, rgeom(1, 1/(1 + MEAN_RECRUD_DELAY)))` days after its suppression. The one-day floor turns
the draw's zeros into ones, so the mean was D + 1/(1 + D): 7.125 days for the default
`TREATMENT_DELAY` of 7, and 10.59 for the default `MEAN_RECRUD_DELAY` of 10.5. `TREATMENT_DELAY = 0`
was accepted and meant the next day.
**Decision.** Both delays are `1 + rgeom(1, 1/D)`: each day after the case, or after the
suppression, the event happens with probability 1/D, so the delay is at least one day and its mean
is exactly D. Both parameters must be finite and at least 1; `validate_params()` refuses other
values, and so does the model's `initialize()` for runs started directly with SLiM.
`TREATMENT_DELAY = 1` means the next day, the fastest response the model can represent. The other
recrudescence shapes are unchanged: `"poisson"` is 1 + Poisson(D - 1), with mean D, and `"fixed"`
is round(D) days with halves rounded up, 11 days at the default 10.5.
**Supersedes.** In D24, the reading of `TREATMENT_DELAY = 0` as the next day, and the rejected
alternative "Forbid `TREATMENT_DELAY = 0`"; D24's same-tick fence in `applyTreatment()` stands. In
D26, the geometric draw and its verification signature; the 10.5-day mean and the selectable
shapes stand.
**Rationale.** Every duration parameter has one reading (D29), and both parameters now mean what
they say. A daily probability of 1/10.5 is Taylor et al. (2019)'s posterior median rate applied
per day.
**Alternatives rejected.** *Keep the draws and list them as exceptions* -- two parameters whose
mean is not their value, for no modelling reason. *Keep `TREATMENT_DELAY = 0` as another name
for 1* -- two values for one behaviour, where 1 already means the next day under the convention.
**Consequences.** Runs that schedule case-based courses (`BLOODCURE_COVERAGE > 0`) or suppress
clones (`P_BLOODCURE_SURVIVAL > 0` and `P_RECRUD_GIVEN_FAIL > 0`, after a course or an MDA round)
change, through slightly shorter delays and the random stream. All other runs are byte-identical,
the defaults and MDA rounds that suppress nothing included. The `TREATMENT_DELAY = 0` runs reported
in D20 and D24 were measured under earlier meanings of 0; the closest valid setting is now 1
(exactly one day), and their counts are not expected to reproduce. At the default 10.5 days the
recrudescence gap has median 7 days and about 94% fall within day 28 (93% under the old draw).
**How to verify.** In a run whose only cases are one pulse of detected imports, and where no
course clears anything (every clone still patent when its course lands is suppressed): at D = 1 no
course and no recrudescence comes later than the next day; at D = 4 the mean delay is 4 and a
quarter of the delays are one day, a share that any mass at zero would inflate
(`tests/testthat/test-treatment.R`). Values below 1 are refused by R and by the model
(`tests/testthat/test-params.R`, `tests/testthat/test-slim.R`). D26's recipe now gives the
geometric signature median 7 / sd 10.0 / P(<=7 d) 0.50 / P(<=28 d) 0.94.

---

## D31 — No age cap on hypnozoite batches · 2026-11-05 · Accepted; partly superseded by D35
**Context.** `MAX_HYPNOZOITE_AGE` (1,500 days) cleared a whole hypnozoite batch once it was older
than the cap. It began at 300 days as a "maximum hypnozoite survival" with no source, and was
raised to 1,500 when that truncated too many batches. At the defaults it binds for roughly one
batch in 35,000, and where it binds it produces an artificial spike of whole-batch clearances
at exactly the cap age. White et al. (2014) give hypnozoites exponential residence with no cap.
**Decision.** Remove `MAX_HYPNOZOITE_AGE`. A batch empties only through activation and death.
`validate_params()` requires `MU_HYPNOZOITE > 0`, so every batch is eventually cleared.
**Rationale.** A finite hypnozoite lifespan is expressed through the daily death probability
`MU_HYPNOZOITE`, not through a hard cutoff that distorts the late-relapse tail it was meant to
protect.
**Alternatives rejected.** *Keep it as a silent safety net* -- it rarely binds, but when it does it
bends exactly the long-latency tail the model needs.
**Consequences.** Runs shorter than 1,500 days are unchanged. `MU_HYPNOZOITE = 0` is no longer
accepted. Until hypnozoites can also die during latency, a temperate batch (`K_LATENCY > 0`)
leaves latency only at rate `DELTA_LATENCY`, so `validate_params()` requires `DELTA_LATENCY`
in (0, 1] whenever `K_LATENCY > 0`.
**How to verify.** `MAX_HYPNOZOITE_AGE` is an unknown parameter and `MU_HYPNOZOITE = 0` is
refused (`tests/testthat/test-params.R`); `events.tsv` and `state.tsv` are byte-identical to the
capped model for three seeded runs of 200-400 days.

---

## D30 — No age cap on blood-stage infections · 2026-11-05 · Accepted
**Context.** `MAX_BROOD_AGE` (180 days) force-cleared any blood-stage infection older than 180
days (the test was `> 180`, so on day 181). It is a leftover from an early version of the model
in which it was the infection duration itself, kept as a ceiling once clearance became a daily
draw. At the defaults it cut about 5% of infections, so the realised mean was about 57 days
rather than 60. Neither White et al. (2014) nor White et al. (2018) caps infection duration.
**Decision.** Remove `MAX_BROOD_AGE`. Every infection ends by the daily clearance draw (D29).
`validate_params()` requires `MEAN_BLOOD_STAGE_DURATION >= 1`, which the 1/D convention needs
anyway and which covers the only case the cap protected against, a duration that disables
clearance.
**Rationale.** With a daily clearance probability of 1/D every infection ends, after a geometric
time with mean D. A hard cap only truncates that distribution and makes the parameter mean
something other than what it says.
**Alternatives rejected.** *Raise the cap to a real safety net (e.g. 365 days)* -- it would still
truncate the tail in long-duration scenarios and guards against nothing validation does not
already catch.
**Consequences.** Infections are about 5% longer on average at the defaults, so prevalence rises
a little; the calibration bands are re-checked before the release. D25's open question of
whether a recrudescence can renew a clone's lifespan without a ceiling no longer arises: a
clone's daily chance of clearing does not depend on its age.
**How to verify.** With founders only, `MEAN_BLOOD_STAGE_DURATION = 1000` and 200 days, no founder
clears on exactly day 181; with the cap about 83% did (`tests/testthat/test-run.R`).

---

## D29 — Durations are geometric with daily probability 1/D · 2026-11-05 · Accepted
**Context.** Duration parameters were turned into daily probabilities in two different ways.
Blood-stage clearance used `1 - exp(-1/MEAN_BLOOD_STAGE_DURATION)`, a continuous-time rate made
discrete, while mosquito death used `1/MOSQ_LIFESPAN`, and biting, hypnozoite activation,
hypnozoite death and exit from latency were daily probabilities already. Under the exponential
form the mean blood-stage infection lasted 60.5 days rather than the 60 the parameter states
(ignoring the 180-day cap that D30 removes, which brought the realised mean to 57.5 days).
Separately, a mosquito bites at most once a day, so `BITING_RATE` times its seasonal and
resurgence multipliers above 1 was clipped at the peak and the realised mean fell below
`BITING_RATE`, without any warning.
**Decision.** Every duration parameter D is the mean of a geometric duration: the event happens
each day with probability 1/D. Blood-stage clearance becomes `1/MEAN_BLOOD_STAGE_DURATION`.
`validate_params()` warns when the daily biting rate the model will use, `BITING_RATE` times
its seasonal, intervention and resurgence multipliers, exceeds 1 on any day of the run. The
convention is stated once, in `?vivax_params` and `MODEL_DESCRIPTION.md` §10.
**Rationale.** One reading for every duration parameter, and each parameter means what it says.
A daily clearance probability of 1/60 is White et al. (2014)'s clearance rate (Table 1, after
Collins et al. 2003) applied per day.
**Alternatives rejected.** *Convert every duration with `1 - exp(-1/D)`* -- natural for an ODE
model, but this model steps in days and its other processes already use daily probabilities.
*Lift the one-bite-a-day cap* -- a mosquito takes at most one blood meal a day here.
**Consequences.** At the defaults the mean blood-stage infection shortens from 60.5 to 60 days,
or from 57.5 to 57.1 with the 180-day cap still in place, so every run's output changes; the
calibration bands are re-checked before the release. Mosquito death, biting, activation,
hypnozoite death and latency exit are unchanged.
**How to verify.** With founders only (no hypnozoites, no transmission) and
`MEAN_BLOOD_STAGE_DURATION = 2`, founder infections last 2 days on average; the exponential form
gives 2.54 (`tests/testthat/test-run.R`).

---

## D28 — MDA keeps its place and gains a coverage parameter · 2026-10-23 · Accepted; partly superseded by D44 (MDA reaches holders and clones patent that day)
**Context.** An early step map proposed withdrawing MDA from the pre-release. The strongest reason
given was that it reused `P_RAD_SURVIVAL`, whose meaning changed in D20, so the two paths would
silently disagree -- and that reason died when D20 routed MDA through `applyTreatment(h, T, "mda")`.
MDA now inherits host-level treatment, D24's patency floor, `RADCURE_MODE` and Step 4's
suppression for free, and it demonstrably works. The one real gap was that a round was
all-or-nothing: every host got the course.
**Decision.** Keep MDA and add `MDA_COVERAGE` (default 1.0), the fraction of the population a
round reaches. Drawn per host, and only when coverage is below 1.
**Rationale.** Coverage is the whole policy question for MDA -- it is what separates a feasible
campaign from an idealised one -- so without it the feature could not answer the question it
exists for. Removal would have cost more than it saved: six event types out of the **exported**
`event_types()` dictionary, `mda_bloodstage_clear` out of the episode vocabulary that span
reconstruction depends on, a metric branch, a shipped validation scenario, the structure diagram,
three passages in `MODEL_DESCRIPTION.md`, and two accepted ADRs needing supersession. There is
also a live mislabel downstream: `deploy_malariafwd_vivax.R` ships a scenario named
`E3_mda_imperfect` / "imperfect coverage" implemented as `P_RAD_SURVIVAL = 0.50`, which conflates
coverage with efficacy. `MDA_COVERAGE` is what that scenario actually wanted.
**Alternatives rejected.** *Remove MDA* -- above. *Leave it all-or-nothing and document the
limitation* -- cheaper, but it leaves the mislabelled scenario with nothing correct to point at.
**Consequences.** `MDA_COVERAGE = 1.0` is byte-identical to the previous behaviour: the draw is
skipped entirely at full coverage, so no existing MDA run changes. Also fixes a false-positive
fence found while testing this: `validate_params()` warned that `RADCURE_COVERAGE` has no effect
without `BLOODCURE_COVERAGE`, which is true of **case-based** treatment but not of MDA -- MDA
delivers radical cure to everyone it reaches. The fence was firing on exactly the scenario MDA
exists for, and now also checks that no usable MDA round is scheduled.
**How to verify.** At `MDA_COVERAGE = 1.0`, `events.tsv` must be identical to the pre-change model
(verified, 2 seeds). Across 0 / 0.25 / 0.5 / 0.75 / 1 the courses given must rise monotonically
(0 / 17 / 24 / 50 / 66) and the reservoir fall (5.39 / 4.38 / 4.12 / 2.30 / 1.38 against 5.44
before the round). The radical-cure fence must stay silent when MDA is scheduled and still fire
when it is not, when coverage is 0, or when the round falls after the run ends.

---

## D27 — The origin vocabulary and the travel window are argument-driven, never content-driven · 2026-10-23 · Accepted
**Context.** Two items were deferred from D23 as needing their own decision rather than being
folded into a schema fix. Both are now measured on the current model.
(a) `force_of_infection(split_imported = "auto")` relabelled relapses as
`relapse (imported)` / `relapse (local)` **only when the run contained imports**, so the returned
vocabulary depended on run content: `baseline` came back with `relapse`, `imports` with
`relapse (local)`, for the same quantity. Pooling a grid did not error -- it produced a table in
which one origin was split across three columns.
(b) `classify_who()`'s import predicate required `is_dormant == 0`, so a traveller who returned
carrying only hypnozoites never registered as having travelled: 43 of 91 import-exposed hosts on
the audit run.
**Decision.** (a) Split unconditionally (`split_imported = TRUE` by default; `"auto"` retired),
and return `origin` as a factor with **all** levels present, so the shape of the result follows
the argument and never the run. (b) Count any import establishment as a trip, dormant or not.
**Rationale.** The split is worth keeping -- an imported relapse is a distinct truth category,
which D20 went out of its way to expose -- so the fix is to stop it being conditional rather than
to remove it. For (b), a surveillance system records the journey, not which compartment the
parasite ended up in; requiring a patent blood-stage infection confuses the traveller with the
infection.
**Consequences, and an important negative result.** (b) is correct in principle but **changes
almost nothing**: flagged cases went 518 to 520, and imported relapses flagged stayed at 142 of
239. The dormant filter was never the binding constraint -- **the 60-day look-back window is**.
Measured, the gap from importation to an imported relapse has median 43 days but IQR 18-95 and
maximum 325, so 40% fall outside a 60-day window however the predicate is written. Widening the
window to 90 / 180 / 365 days recovers 73% / 92% / 100% of imported relapses while flagging
520 -> 685 / 935 / 1111 cases in total. That trade is now documented on `classify_who()`, because
it is the quantity the function exists to measure rather than a defect in it: no travel window
recovers a hypnozoite that activates a year after the trip. Anything previously computed with
`classify_who()` moves by two cases; anything that pooled `force_of_infection()` across runs with
and without imports was silently splitting one origin into three columns and should be re-run.
**How to verify.** `force_of_infection()` must return the same `levels(origin)` for a run with
imports and one without (5 split / 4 unsplit), with `is_acquisition` correct for the split labels;
a pooled `dcast` across both must produce one column per origin. `classify_who()` must flag a host
whose only import establishment was dormant. Suites: `test_dir()` 97/0,
`audit_events_and_ibd.R` 99/0, `R CMD check --as-cran --no-manual` 3 NOTEs.

---

## D26 — The recrudescence delay is geometric at mean 10.5 d, and the distribution is selectable · 2026-10-23 · Accepted; partly superseded by D32 and D42
**Context.** Step 4 shipped the suppression-to-recrudescence delay as a geometric at mean 21 d.
Two things needed settling: the plan (`TREATMENT_IMPLEMENTATION.md:32, :361`) specified
`1 + Poisson(MEAN_RECRUD_DELAY - 1)` instead, and the 21 d default was uncited. A literature and
reference-model review found that **none of the four reference models represents recrudescence
at all** -- White 2014 has no treatment-failure branch; the White 2018 IBM puts a failure
straight into clinical disease (`I_D`) with no delay; OpenMalaria's failure branch is literally
commented "No change in parasitological status"; and malariasimulation routes failures to the
untreated state and *declares* `late_parasitological_failure_probability` but never reads it.
So the reference models cannot arbitrate this, and the model goes beyond all four.
**Decision.** Keep the geometric and fix the mean: **`MEAN_RECRUD_DELAY` 21.0 -> 10.5**. Add
`RECRUD_DELAY_DIST` with `"geometric"` (default), `"poisson"` (the plan's shape, as a
tight-window sensitivity arm) and `"fixed"` (deterministic).
**Rationale.** Geometric is the discrete analogue of the only published fit for this quantity:
Taylor et al. (2019) model time-to-recrudescence as Exponential, and their own simulator draws
`round(rexp(1, rate = 1/10))` -- a rounded exponential on a day scale *is* a geometric, so the
implementation was already faithful to the reference. It is also the house style of every
reference model's state exits. 10.5 d is the reciprocal of Taylor's posterior median, replacing
an uncited number with a traceable one. The old 21 d was not merely unsourced but
self-contradicting: under a geometric it puts **26% of recrudescences beyond day 28**, outside a
standard WHO therapeutic-efficacy window, where a real study could never observe them -- while
the code comment claimed it "sits inside a 28-day follow-up window". Measured at the new default:
median gap 7 d, **92% inside day 28** (against 74% at 21 d).
**Alternatives rejected.** *The plan's Poisson* -- no source in any reference model, in Taylor,
or anywhere found, and it encodes an implausible claim: at mean 21 its SD is 4.5 d and
P(gap <= 7 d) is 0.0003, i.e. an almost-fixed delay, where the underlying process is subpatent
parasite regrowth from an unobserved starting density. Retained as a selectable arm rather than
deleted, because a tight-window contrast is a useful sensitivity. *A gamma with a shape
parameter* -- there is a real biological argument for a peaked delay with a hard floor (regrowth
is log-linear, so the delay is a logarithm of surviving biomass, and the one large empirical
hazard estimate for falciparum is non-monotonic), and a reviewer suggested mean 12 / min 5 /
shape 4. Not shipped: it needs a second parameter, and the family with a published fit behind it
should be the default. A clean extension if the biology argument wins later.
**Consequences.** Runs with `P_RECRUD_GIVEN_FAIL > 0` change -- the delay halves. Runs at the
shipped default are unaffected and remain byte-identical to the pre-Step-4 model (re-verified
after this change on 4 runs). `RECRUD_DELAY_DIST` joins `.malariafwd_string_params` and is fenced
on both sides. The plan documents' Poisson spec is now superseded and must be marked historical.
**How to verify.** Each option produces its own signature on a run with
`P_RECRUD_GIVEN_FAIL = 1`, pairing each recrudescence to its holder via `parent1_pedigree_id`
(the holder keeps its own id; the recurrence mints a new one): geometric median 7 / sd 11.8 /
P(<=7 d) 0.51 / P(<=28 d) 0.92; poisson median 10 / sd 3.3 / P(<=7 d) 0.19; fixed sd 0. An
invalid value must be refused on both the R and SLiM sides.

---

## D24 — Treatment cannot precede patency · 2026-10-23 · Accepted; partly superseded by D32 and D44 (MDA courses)
**Context.** `decideTreatment()` scheduled a course at `sim.cycle + rgeom(1, 1/(1+TREATMENT_DELAY))`.
A geometric with mean *D* returns 0 with probability `1/(1+D)` -- **12.5% at the shipped default
of 7**, and 100% at 0 -- so the course landed in the same tick's `early()`. `reproduction()` runs
before `early()`, so the genome that caused the case was minted moments earlier in the same tick,
and `applyTreatment()` cleared it on its own creation tick. Measured over 5 seeds on three
configurations: **796 episodes had zero patent days** (426 / 184 / 186), i.e. a parasite was
treated before it had been detectable for a single day, could cause symptoms, or could transmit.
This is not a treatment question, it is an internal inconsistency: natural clearance has enforced
the rule since the model was written (`malariafwd_vivax.slim`, "Cannot be cleared on same tick as
inoculation"); `applyTreatment()` simply never applied it.
**Decision.** Two changes, both needed. (1) `decideTreatment()` floors the delay at **one tick**:
`sim.cycle + max(1, rgeom(...))`, and `TREATMENT_DELAY = 0` now means "as fast as the model can
represent", i.e. next tick. (2) `applyTreatment()` skips any clone whose `bloodstage_start_tick`
is the current tick, so a course given for one case cannot reach a *different* clone that appeared
that same tick.
**Rationale.** A clone that became patent this tick cannot have caused symptoms, prompted
care-seeking or been diagnosed -- OpenMalaria and malariasimulation both place delays between
infection, symptoms and treatment. The floor must be at the decision site as well as the
application site: with only the `applyTreatment()` fence, a course whose every target was created
that tick is **consumed anyway** (`tx_pending` is cleared regardless of what was treated) and the
patient silently goes untreated. Measured, that overcorrection moved prevalence 24.4% -> 55.3% on
the `TREATMENT_DELAY = 0` run. With both changes the effect is small and in the expected direction.
**Alternatives rejected.** *Forbid `TREATMENT_DELAY = 0`* -- restricts the parameter space for no
reason; zero is a coherent request and now resolves to the fastest representable response.
*Fence only in `applyTreatment()`* -- wastes courses, as above. *Leave it and document* -- the
model would keep asserting two contradictory rules about the same tick.
**Consequences.** Every treated run changes. Over 5 seeds: zero-patent-day episodes 796 -> **0**;
prevalence 14.3+-8.0 -> 15.4+-4.0 (`treat_now`), 52.0+-10.7 -> 56.3+-6.2 (`cotrans`),
35.1+-5.4 -> 35.7+-5.8 (`mild`) -- every shift well inside the between-seed SD, so no calibration
conclusion moves. Treatment is now marginally *less* effective, which is the expected direction
once a parasite gets at least one patent day. Untreated runs are unaffected
(`BLOODCURE_COVERAGE = 0` is the shipped default, so `applyTreatment()` never runs).
Also makes the radical-cure efficacy check 3-seed. It is not a regression -- the ordering is
correct on both the old and new model, and *stronger* on the new one (6 of 6 seeds against 5 of 6)
-- but single-seed it was a near coin flip, so the shifted RNG stream flipped it with no change in
behaviour. Same fix as D20 applied to its sibling check; the two remaining single-seed directional
checks in `check_parameters.R` carry the same fragility.
**How to verify.** No episode may be cleared by `case_bloodstage_clear` or `mda_bloodstage_clear`
on the tick it was created, on any configuration (0 across `treat_now`, `cotrans`, `mild` at 5
seeds each, against 426/184/186 before). `tx_tick > sim.cycle` for every scheduled course. Suites:
`test_dir()` 71/0, `audit_events_and_ibd.R` 66/0, `check_parameters.R` 41/0,
`R CMD check --as-cran --no-manual` 3 NOTEs (the two known, plus "unable to verify current time",
which is the check failing to reach the network in this sandbox).

---

## D25 — Three-state recurrence: a recrudescence mints a new genome and is a case, not an acquisition · 2026-10-23 · Accepted; implemented except the subpatent_clearance event; relapse-fraction names by D65, the molFOB pair by D69; item 3's Champagne sentence superseded by D71; open item on MDA settled by D44
**Context.** Step 4 adds the third recurrence state of the Pv3Rs/Taylor taxonomy -- relapse,
reinfection, **recrudescence** -- by letting a clone that survives a blood-stage course be
suppressed below patency and later re-expand. A scratch prototype
(`/tmp/mfwchk/step4/proto.slim`) established the mechanics; a second prototype established the
representation. Both are throwaway and are not in the repo.
**Decision.**
1. **A recrudescence MINTS a new genome** -- `h.addCloned(holder)`, fresh `pedigree_id`, fresh
   tskit node at the recrudescence tick -- and the holder is retired. This is the relapse
   precedent exactly. Decision and mint share `reproduction()` (the holder already exists, so
   unlike D17 there is no deferral); retirement is the first thing in `early()`.
2. **It is a CASE but not an ACQUISITION.** `classify_cases()` gains a fourth `origin`, so case
   counts, incidence and the host timeline include it. The clone-acquisition metrics exclude it:
   FOI = imported + local reinfection, molFOB = FOI + relapse.
3. **Metrics that have both readings report both**, named by unit: `molfob_clones` against
   `molfob_episodes`, `relapse_frac_acquisitions` against `relapse_frac_episodes`, and
   `rc_by_origin(unit = c("episode","acquisition"))`. The published 0.55-0.85 band and the
   Champagne/Swiss TPH Rc comparison are both **acquisition**-level quantities.
4. `inoculation_id` and `oocyst_id` on the recrudescence row are **-1**, following the relapse.
5. `ends_episode`: `bloodcure_subpatent = TRUE`, `recrudescence = FALSE`,
   `subpatent_clearance = FALSE`.
6. `event_types()` gains an **`is_onset`** column so the onset vocabulary becomes derivable, as
   `ends_episode` already is.
7. **Patency, infectiousness and capacity are separated.** `getActiveLineageCount()` currently
   answers three different questions at three call sites, and subpatency splits them apart:
   `getActiveLineageCount()` = PATENT (the state writer, `moi_active`, episodes);
   `getTransmissibleLineageCount()` = can infect a mosquito (`h.tag`, and `human_MOI` in the bite
   loop); `getBloodStageLineageCount()` = present at all (the `MAX_MOI` capacity check). This is
   **required**, not optional: with them fused, `SUBPATENT_INFECTIVITY` cannot do anything, and
   the failure is invisible. **REVISIT AND RE-VERIFY THIS DURING IMPLEMENTATION** -- it was got
   wrong twice in prototype, see below.

   The subtlety that makes it a trap: splitting `h.tag` alone is NOT enough. The bite loop's
   `human_MOI` is drawn from the patent count and gates Event 2 with
   `if (human_MOI == 0) { "tagged infected but MOI=0; skipping" }` -- which reads as a defensive
   check against an impossible state, but becomes a reachable and legitimate state once a
   subpatent compartment exists, and silently blocks exactly the transmission the parameter is
   supposed to permit. Measured in prototype, uptakes sourced from a host with NO patent clone:
   **0 / 0 / 0** under the fused design at `SUBPATENT_INFECTIVITY` 0 / 0.5 / 1.0, and
   **0 / 68 / 233 / 918** (4.4% of all uptakes at the literature value of 0.15, 30% at 1.0)
   once all three roles are separated. The intermediate state -- `h.tag` split but `human_MOI`
   not -- looks exactly like "the subpatent compartment is too small to matter", which is the
   wrong conclusion and the one I drew first.
**Rationale.** Minting was chosen for pedigree uniqueness, not for node time. The node-time
argument is real but marginal -- staleness is median 19 ticks, and `arg_tmrca()` has zero callers
in the analyses repo. Reuse is what actually breaks things: `creation events mint unique
pedigree_ids` fails with 33/42/15 duplicates per run, `R/metrics.R:207`'s stated invariant becomes
false, and `rc_by_origin()` double-books 53-70% of all edges. Minting costs exactly +1 individual,
+2 nodes, +1 full-span edge, **+0 local trees** and 190 bytes (measured over 2035 relapse
`addCloned` children with zero variance) = 0.4% of individuals at the prototype rate. MOI,
prevalence, COI and the reservoir are *exactly* neutral, because the minted genome inherits
`lineage_id`; a fresh id would give a phantom +1 MOI. The case/acquisition split follows the unit
discipline D23 established: the difference between the two molFOB readings is itself the
quantity -- intermittent genotyping cannot tell a re-expanded clone from a newly acquired one, so
a real molFOB over-counts, and the model can measure by how much.
**Alternatives rejected.** *Reuse the Individual* -- above. *`treeSeqRememberIndividuals()` at the
recrudescence tick to re-date the node* -- killed by live probe: the node table comes back
identical and depth is still 0. *Mint without retiring the holder* -- one treatment failure
becomes a clone generator (2.21 recrudescences per holder, up to 6). *A single relapse-fraction
denominator* -- the two readings answer different questions and the band only applies to one.
**Consequences.** Nine origin-vocabulary sites in `R/plots.R` and `R/segments.R` need the fourth
level; `fill0()` at `R/plots.R:156` is the dangerous one, because it left-joins against a 3-origin
`CJ` and would drop recrudescence rows silently. Four audit assertions and one test change; only
`every clearance matches an onset` is new to minting, and (6) is its fix. The sentinel fix
(`181fa31`) and the patency fix (D24) were landed first so Step 4 starts from a clean base.
Literature anchors the defaults: recrudescence is **0.2-0.3% of recurrences** after chloroquine or
artesunate (Taylor 2019, n=1441), so the shipped default must be ~0 and a chloroquine-resistant arm
is needed for the state to have anything to classify; subpatent relative infectivity is
**0.10-0.20** (Koepfli 2021; White 2018 IBM), so neither 0 nor 1 is defensible.
**Open, to settle during implementation.** The plan's `1 + Poisson(MEAN_RECRUD_DELAY - 1)`
(`TREATMENT_IMPLEMENTATION.md:32, :361`, default 21 d) carries **no cited source**; it is
self-consistent with the plan's row 7 ("natural clearance while subpatent: off", which makes
`P_BLOODCURE_SURVIVAL` = P(recrudescence | treated) exactly), but Taylor's `RecurrentVivax` models
time-to-recrudescence as Exponential with mean 10 d (SD 10) against Poisson(20)'s SD 4.5. The
prototype used competing hazards, which is a real deviation from row 7 and breaks that identity.
Decide: pure delay with the clean parameter identity, or competing risks with a more mechanistic
subpatent exit. Also open: whether a later course may reach an already-subpatent holder; what MDA
does to a holder; whether repeated cycling may renew a clone's lifespan without a ceiling (minting
resets `bloodstage_start_tick`); and `P_DETECT_RECRUD`, where the plan's `= P_DETECT_PRIMARY`
is contradicted by Ross 2025 (recurrence detectability 0.13-0.49 of primary).
**How to verify.** `P_SUBPATENT_GIVEN_FAIL = 0` must stay byte-identical to the shipped model
(verified on the prototype). Every recrudescence must be IBD 1.000 to its holder and
`arg_tmrca` depth 0. `creation events mint unique pedigree_ids` must pass. MOI and prevalence must
be unchanged against a matched no-recrudescence run. The two molFOB readings must differ by exactly
the recrudescence count. For item 7 specifically, the test is mechanistic rather than
population-level: count mosquito uptakes sourced from a host with **no patent clone**, joining
each uptake at tick T to the state snapshot at T-1 (the bite loop runs in `reproduction()`, so it
sees the previous tick's end state). That count must be 0 at `SUBPATENT_INFECTIVITY = 0` and rise
monotonically above it. Do NOT use patent prevalence as the test -- at four seeds it came out
5.5 / 5.9 / 2.8 / 10.7% with SDs of 4.0 / 5.4 / 0.9 / 3.6, i.e. unresolved, while the uptake
channel was unambiguous.

---

## D23 — Metrics derive their vocabulary from classify_cases(), and output schemas do not depend on run content · 2026-10-22 · Accepted; the relapse_fraction() clause superseded by D65
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

## D22 — IBD counts all ancestry; recency is the caller's filter over returned segments · 2026-10-21 · Accepted; partly superseded by D63 (the shape of ibd_pairs()'s result) and D74 (the segment definition)
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

## D21 — Windowed series report the days they actually cover; onset filters derive from the event dictionary · 2026-10-20 · Accepted; partly superseded by D61 (the run end is SIM_DAYS) and D62 (every window; the first holds bin - 1 days)
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

## D19 — Case management gets its own window, separate from vector control · 2026-10-19 · Accepted; partly superseded by D46 (vector control does not scale imports)
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

## D2 — 14 chromosomes as ONE concatenated linkage group (Option B) · 2026-08-05 · Accepted; partly superseded by D34 (the concatenated map stands) and D74 (a shared stretch can run across a junction in the tree sequence; the IBD functions cut it there)
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

---

## Appendix — pre-release correctness fixes (moved out of NEWS.md)

These were fixed during development, before any public release. They were listed in `NEWS.md`
until the 0.1.0 rewrite, which took the view that a first-release changelog should describe what
the package **does**, not what changed on the way there — a new user never saw the broken
versions, so the fixes are not news to them. Kept here so the record is not lost.

Four of them stem from one Eidos scoping rule: a plain assignment to a global inside a callback
or function creates a *local*, so the write was never visible.

* **Hypnozoite batches never died.** Batches that were exhausted or hit the age cap kept
  activating with their old count — in one check run **64% of batches (74/116) drew more
  hypnozoites than they ever held** (636 excess), inflating relapses. Batches are now retired
  with `batch_count = 0` and killed in the `early()` cleanup, and a regression test asserts no
  batch over-draws. Relapse-related calibration from before this fix was re-run.
* **`inoculation_id` restarted every tick** (max id 3 across a whole run); now unique.
* **`event_id` was always `<tick>_1`** (6745 rows collapsing to 1073 ids); now unique.
* **Relapse-triggered radical cure was mostly undone** — it thinned the batch, then the relapse
  bookkeeping overwrote it. It now runs after the host's batches are updated.
* **MDA errored** when a treated host was left with no infections; that host is now skipped.
* **`ibd_index()` scored a genome against itself as 0.** tskit only stores distinct pairs, so a
  self-pair fell through to the below-detection default, which would make the diagonal of any
  pairwise matrix silently wrong. It now returns 1, and no longer carries node-map names on its
  return value.

The API was also renamed wholesale before release, to one prefix per family (`ibd_`, `arg_`,
`load_`, `plot_`, `slim_`, `uptake_`), with `run_` meaning only "execute". The old names
(`calc_ibd`, `run_params`, `run_metrics`, `who_classify`, `zt_moments`, ...) never shipped
publicly and are not aliased.
