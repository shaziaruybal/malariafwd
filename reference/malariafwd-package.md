# malariafwd: concepts and terms

malariafwd simulates *Plasmodium vivax* transmission and genetics
forward in time in SLiM, one host, mosquito and parasite genome at a
time, and reads what a run records from R. This page defines the terms
the help pages share.

## Start here

[`example_run()`](https://shaziaruybal.github.io/malariafwd/reference/example_run.md)
returns a small run that ships with the package, so every function that
reads a run can be tried without SLiM.
[`vignette("malariafwd")`](https://shaziaruybal.github.io/malariafwd/articles/malariafwd.md)
shows how to run a simulation and reads the bundled run,
[`vignette("parameters")`](https://shaziaruybal.github.io/malariafwd/articles/parameters.md)
gives the parameters that are easy to get wrong and how to choose them,
and
[`vignette("anatomy")`](https://shaziaruybal.github.io/malariafwd/articles/anatomy.md)
walks through what a run records (all three are also on the
[website](https://shaziaruybal.github.io/malariafwd/));
[`event_types()`](https://shaziaruybal.github.io/malariafwd/reference/event_types.md)
lists every event a run logs. The model is described in full in
[MODEL_DESCRIPTION](https://shaziaruybal.github.io/malariafwd/MODEL_DESCRIPTION.html),
the reasons for its choices are in [design
decisions](https://shaziaruybal.github.io/malariafwd/DECISIONS.html),
and the scripts that validate the model are installed in
`system.file("validation", package = "malariafwd")`.

## Runs and their tables

A *run* is one simulation:
[`run_sim()`](https://shaziaruybal.github.io/malariafwd/reference/run_sim.md)
writes it to a directory, and
[`load_run()`](https://shaziaruybal.github.io/malariafwd/reference/load_run.md)
reads that directory back as a `vivax_run`. Every function that reads a
run takes either. A run's outputs are four files, each read by its own
function: `params.tsv`, the parameter set it was made from
([`load_params()`](https://shaziaruybal.github.io/malariafwd/reference/load_params.md));
`events.tsv`, one row per event
([`load_events()`](https://shaziaruybal.github.io/malariafwd/reference/load_events.md),
described in
[events.tsv](https://shaziaruybal.github.io/malariafwd/reference/malariafwd-outputs.md));
`state.tsv`, every host and mosquito on every day
([`load_state()`](https://shaziaruybal.github.io/malariafwd/reference/load_state.md),
hosts by default); and `plasmo.trees`, the genealogy of every genome as
a tree sequence
([`trees_path()`](https://shaziaruybal.github.io/malariafwd/reference/trees_path.md)).
An argument named `events` takes an events table or a run, one named
`state` a state table or a run, and a plot takes a run or the table its
function returns.

## Time

A *tick* is one day. A run's days are 1 to `SIM_DAYS`, and every column
named `*_tick` is a day. Most durations are geometric: a state with a
mean duration of `D` days ends on each day with probability `1/D`.
`EIP_DAYS` and `GAMETOCYTE_DELAY` are fixed delays, `RECRUD_DELAY_DIST`
sets the shape of a recrudescence's delay, and with `K_LATENCY` above 0
a batch's latency is a sequence of geometric stages. A year is 365 days,
so a rate per person-year divides by person-days over 365.

## Genomes and identifiers

Each parasite genome is one individual of the simulation. A bite, an
import, a relapse and a recrudescence each make a new genome, and a
batch of hypnozoites is one genome that stands for its count.
`pedigree_id` names a genome, uniquely within a run; it is 0 or more,
and -1 means none. `lineage_id` names a clonal group: a genome keeps its
parent's lineage through clonal and selfed descent, and a recombinant
starts a new one. In this package *lineage* always means this group; the
genealogy itself is the tree sequence. `oocyst_id` names one meiosis in
a mosquito, `human_id` and `mosq_id` a host and a mosquito, and
`entity_id` either of them in `state.tsv`.

## Cases, episodes and acquisitions

A *case* is a clone onset, one row of
[`classify_cases()`](https://shaziaruybal.github.io/malariafwd/reference/classify_cases.md):
a bite or an import that gave a blood-stage infection, a relapse, or a
recrudescence. Founder seeding, where a run starts, is not a case. An
*episode* is the clones that became patent in one host on one day, named
by `episode_id`, `"<human_id>:<tick>"`, so a bite that delivered two
genotypes is two cases of one episode
([`episodes()`](https://shaziaruybal.github.io/malariafwd/reference/episodes.md),
[`recurrences()`](https://shaziaruybal.github.io/malariafwd/reference/recurrences.md)).
An *acquisition* is a case that began a new blood-stage clone in its
host: an import, a local reinfection or a relapse. A recrudescence is a
case but not an acquisition, since it re-expands a clone already
counted. A relapse into a patent clone of its own lineage is flagged by
`into_patent_lineage`, and `exclude_masked = TRUE` leaves such relapses
out of the counts.

## COI and MOI

Here the complexity of infection (COI) and the multiplicity of infection
(MOI) are one quantity: the distinct patent blood-stage lineages in a
host, `moi_active` in `state.tsv`. `MAX_MOI` gates new inoculations: a
host already carrying that many blood-stage lineages, patent or
subpatent, takes no new bite or import. Relapses bypass it, and a bite
that delivers several genotypes is admitted whole, so a host's COI can
exceed it.

## Transmission

A *transmission link* joins a parent genome of an oocyst to a genome its
products infected: `direct`, through a bite's blood stage, or `liver`,
through a relapse of the batch a bite left.
[`transmission_links()`](https://shaziaruybal.github.io/malariafwd/reference/transmission_links.md)
lists them, each with its *chain*, a connected group of links. *Rc*, the
case reproduction number, is the mean number of secondary cases per
case, by true origin
([`rc_by_origin()`](https://shaziaruybal.github.io/malariafwd/reference/rc_by_origin.md)).
[`classify_chain()`](https://shaziaruybal.github.io/malariafwd/reference/classify_chain.md)
classes each case by its true chain from an import, and
[`classify_travel()`](https://shaziaruybal.github.io/malariafwd/reference/classify_travel.md)
by what a rule that knows only travel history would conclude.

## Relatedness

Every genome descends from root genomes, which are unrelated: the
founder templates, one per founder lineage (`N_FOUNDER_LINEAGES`), and
the imports. The founder hosts of one lineage carry clonal copies of its
template. So the genealogy is a forest. Identity by descent (IBD) is
read from it exactly:
[`ibd_pairs()`](https://shaziaruybal.github.io/malariafwd/reference/ibd_pairs.md)
for pairs,
[`ibd_segments()`](https://shaziaruybal.github.io/malariafwd/reference/ibd_segments.md)
for the stretches behind them, and
[`ibd_lookup()`](https://shaziaruybal.github.io/malariafwd/reference/ibd_lookup.md)
to read a pair.
[`arg_tmrca()`](https://shaziaruybal.github.io/malariafwd/reference/arg_tmrca.md)
gives how far back two genomes share an ancestor, and
[`arg_summary()`](https://shaziaruybal.github.io/malariafwd/reference/arg_summary.md)
the shape of the forest.

## See also

Useful links:

- <https://github.com/shaziaruybal/malariafwd>

- <https://shaziaruybal.github.io/malariafwd/>

- Report bugs at <https://github.com/shaziaruybal/malariafwd/issues>

## Author

**Maintainer**: Shazia Ruybal-Pesántez <shazia.ruybal@isglobal.org>
([ORCID](https://orcid.org/0000-0002-0495-179X)) \[copyright holder\]

Authors:

- Shazia Ruybal-Pesántez <shazia.ruybal@isglobal.org>
  ([ORCID](https://orcid.org/0000-0002-0495-179X)) \[copyright holder\]
