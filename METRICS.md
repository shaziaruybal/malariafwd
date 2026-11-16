# malariafwd v0.1.0 — derived metrics

What each number the package computes from a run means, and how it is counted. The output files
these are read from are described in section 14 of the
[model description](https://shaziaruybal.github.io/malariafwd/MODEL_DESCRIPTION.html), and the
[reading the outputs](https://shaziaruybal.github.io/malariafwd/articles/anatomy.html) article works
through them on a real run.

## Summary and case metrics

**Summary metrics.** `summarise_run()` gives each over one window of days, by default the second
half of `SIM_DAYS` (`burn_frac = 0.5`; `from` and `to` set it directly), whether or not the run went
extinct: a day after an extinction counts as a day with no parasite. Counts are clone onsets, as in
`classify_cases()`. The EIR, FOI and molFOB divide by the person-years of every host, `N_HOSTS` ×
the window's days ÷ 365; the relapse rate divides by the person-years at risk. Nothing is rounded,
and a ratio with nothing to divide is `NA`. `exclude_masked = TRUE` leaves the relapses into a
patent clone of their own lineage (`into_patent_lineage = 1`, section 5 of the model description)
out of every relapse column. To pool runs, sum the counts and the person-years rather than average the rates.

| metric | definition | `summarise_run()` column |
|---|---|---|
| prevalence | fraction of the window's host-days with `moi_active > 0` | `prevalence_frac` |
| COI / MOI | `moi_active` over the infected host-days; polyclonality is the share of them with 2 or more | `coi_mean`, `coi_median`, `polyclonal_frac` |
| reservoir | `hyp_reservoir_size` per host-day | `hypnozoites_per_host` |
| EIR | infectious bites per person-year: the sum of `infectious_bites` over hosts and days ÷ person-years (section 4.1 of the model description) | `eir_per_py` |
| FOI | new blood-stage clones from bites and imports per person-year — *not* the same as molFOB | `foi_per_py` |
| molFOB | new blood-stage **clone acquisitions** per person-year, relapses included (Koepfli et al. 2013). A recrudescence is a new episode but not a new clone, so only the episode reading counts it. `force_of_infection()` gives the same counts window by window: sum its `is_acquisition` rows | `molfob_clones_per_py`, `molfob_episodes_per_py` |
| relapses per person-year | relapses ÷ person-years at risk: each host from the day after its first establishment (a primary or a hypnozoite batch, from a bite or an import; founder seeding does not start it) to the window's end, also after an extinction | `relapses_per_py_at_risk` |
| relapse fraction | two readings, counted per clone onset (section 9a of the model description): relapses ÷ **acquisitions**, the unit of the relapse-fraction calibration target (section 13 of the model description), and relapses ÷ **every case**, recrudescence included; equal when a run has no recrudescence. Both include each host's first case, which is not a recurrence. ν was fitted with every relapse counted; leaving the masked ones out gives a reading below its target (section 13 of the model description) | `relapse_frac_acquisitions`, `relapse_frac_episodes` |
| imported share | imports among the new clones from bites and imports | `imported_frac` |
| Rc | the case reproduction number (Fraser 2007), by the case's true origin: secondary cases per case, through bites and the relapses they leave; each secondary case counts once, 1/2 to each parent of an outcrossed oocyst. Cases near the end are right-censored, above all through relapses; `from` and `followup` choose the onsets counted (`rc_by_origin()`) | — |
| time to first relapse | per hypnozoite batch, from the day it was established (for an import that also brought a blood-stage infection, the day after) to its first relapse; censored on the day it empties without relapsing, or at `SIM_DAYS` (`time_to_relapse()`, `emptied = "follow"` to follow emptied batches to the end) | — |

**Truth tables.** Read from the event log, so exact.

| table | definition |
|---|---|
| episode | the clones that became patent in one host on one day, named by `episode_id` ("<host>:<day>"). `classify_cases()` gives one row per clone onset, founder seeding excluded; `episodes()` gives every clone's span, founders included |
| pair IBD | `ibd_pairs()`: the fraction of the genome over which two genomes descend from a common ancestor, counting every shared stretch of any age and length |
| IBD segment | a maximal stretch over which two genomes descend from one common ancestor genome, cut at each chromosome junction, the pieces shorter than `min_span` left out after the cut; the same whatever other genomes are passed, from `ibd_pairs(segments = TRUE)` and `ibd_segments()` alike |
| segment age | `depth`: how long before the older of the two genomes the segment's common ancestor lived, in ticks |
| recurrence | every episode of a host after its first. `recurrences()` pairs each of its clones with each clone of the earlier episodes (or of the latest), never with a clone of its own episode |
| transmission link | a parent genome of an oocyst and a genome the oocyst's products infected: `direct` (a bite's blood stage) or `liver` (a relapse of the batch a bite left). Only local infections are targets; each parent of an outcrossed oocyst has weight 1/2. `transmission_links()` |
| import class (WHO) | each case's class by its true chain from an imported infection: `imported` (an import, a relapse of an import's batch, or a recrudescence of either), `introduced` (one link from one), `indigenous (from an import)` (two or more) or `indigenous (local lineage)` (no import in its chain), through both routes; an outcrossed oocyst's shorter chain counts. `classify_chain()` |
| travel-history class | what a rule that knows only travel concludes, to be scored against the import class: `imported (recent travel)` when the host had an import established, as a blood stage or as hypnozoites only, in the `window` days up to the case, and `local (no recent travel)` otherwise. `classify_travel()` |
| transmission chain | a connected component of the links, each taken in either direction, a recrudescence counted as the clone it came from: `chain_id` and `chain_size` on `transmission_links()`. A genome with no link is in no chain |
