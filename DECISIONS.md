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
**Context.** At the defaults about 85% of relapses start while a patent clone of the same lineage
is already in the host (70-82% in review M3's larger runs). Over four 3-year default runs that
clone was the primary from the same bite for 43% of them, an earlier relapse of the same batch for
45%, another batch's relapse for 10% and another inoculation of the lineage for 2%. Each such
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
defaults, 86% of 11,676 activations over eight 3-year runs are flagged.

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

## D38 — Imports go to uniformly random hosts, with the MAX_MOI check a bite gets · 2026-11-05 · Accepted
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
treatment of imported cases, are more common. Runs without imports, and a founder-free pulse into
a wholly susceptible population, are byte-identical; at `MAX_MOI = 0`, now refused, that pulse
was refused outright. The bite gate changes only runs with `SUBPATENT_INFECTIVITY = 0` in which a
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

## D37 — A run that goes extinct records the tick and says so · 2026-11-05 · Accepted
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
founder-free set whose only import falls outside the run, a pulse on day 0 for instance, is refused
before SLiM starts instead of stopping at the liveness check on tick 10. Every other run is
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
Tropical runs (`K_LATENCY = 0`, the default) are byte-identical. Since survival favours the batches
that move through quickly, a hypnozoite that comes through latency alive has been latent
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

## D33 — Legitimate edge values run to the end; values the model cannot run are refused first · 2026-11-05 · Accepted
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

## D28 — MDA keeps its place and gains a coverage parameter · 2026-10-23 · Accepted
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

## D26 — The recrudescence delay is geometric at mean 10.5 d, and the distribution is selectable · 2026-10-23 · Accepted; partly superseded by D32
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

## D24 — Treatment cannot precede patency · 2026-10-23 · Accepted; partly superseded by D32
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

## D25 — Three-state recurrence: a recrudescence mints a new genome and is a case, not an acquisition · 2026-10-23 · Accepted (design; implementation pending)
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

## D2 — 14 chromosomes as ONE concatenated linkage group (Option B) · 2026-08-05 · Accepted; partly superseded by D34 (the concatenated map stands)
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
