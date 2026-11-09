#' malariafwd: concepts and terms
#'
#' malariafwd simulates *Plasmodium vivax* transmission and genetics forward in time in SLiM, one
#' host, mosquito and parasite genome at a time, and reads what a run records from R. This page
#' defines the terms the help pages share.
#'
#' @section Start here:
#' [example_run()] returns a small run that ships with the package, so every reader can be tried
#' without SLiM. `vignette("malariafwd")` runs a simulation and reads it, and `vignette("anatomy")`
#' walks through what a run records (both are also on the
#' [website](https://shaziaruybal.github.io/malariafwd/)); [event_types()] lists every event a run
#' logs. The model is described in full in
#' [MODEL_DESCRIPTION](https://shaziaruybal.github.io/malariafwd/MODEL_DESCRIPTION.html), the
#' reasons for its choices are in
#' [DECISIONS](https://shaziaruybal.github.io/malariafwd/DECISIONS.html), and the scripts that
#' validate the model are installed in `system.file("validation", package = "malariafwd")`.
#'
#' @section Runs and their tables:
#' A *run* is one simulation: [run_sim()] writes it to a directory, and [load_run()] reads that
#' directory back as a `vivax_run`. Every function that reads a run takes either. A run's outputs
#' are four files, each with its reader: `params.tsv`, the parameter set it was made from
#' ([load_params()]); `events.tsv`, one row per event ([load_events()], described in
#' [events.tsv][malariafwd-outputs]); `state.tsv`, every host and mosquito on every day
#' ([load_state()], hosts by default); and `plasmo.trees`, the genealogy of every genome as a tree
#' sequence ([trees_path()]). An argument named `events` takes an events table or a run, one named
#' `state` a state table or a run, and a plot takes a run or the table its function returns.
#'
#' @section Time:
#' A *tick* is one day. A run's days are 1 to `SIM_DAYS`, and every column named `*_tick` is a day.
#' Most durations are geometric: a state with a mean duration of `D` days ends on each day with
#' probability `1/D`. `EIP_DAYS` and `GAMETOCYTE_DELAY` are fixed delays, `RECRUD_DELAY_DIST` sets
#' the shape of a recrudescence's delay, and with `K_LATENCY` above 0 a batch's latency is a
#' sequence of geometric stages. A year is 365 days, so a rate per person-year divides by
#' person-days over 365.
#'
#' @section Genomes and identifiers:
#' Each parasite genome is one individual of the simulation. A bite, an import, a relapse and a
#' recrudescence each make a new genome, and a batch of hypnozoites is one genome that stands for
#' its count. `pedigree_id` names a genome, uniquely within a run; it is 0 or more, and -1 means
#' none. `lineage_id` names a clonal group: a genome keeps its parent's lineage through clonal and
#' selfed descent, and a recombinant starts a new one. In this package *lineage* always means this
#' group; the genealogy itself is the tree sequence. `oocyst_id` names one meiosis in a mosquito,
#' `human_id` and `mosq_id` a host and a mosquito, and `entity_id` either of them in `state.tsv`.
#'
#' @section Cases, episodes and acquisitions:
#' A *case* is a clone onset, one row of [classify_cases()]: a bite or an import that gave a
#' blood-stage infection, a relapse, or a recrudescence. Founder seeding, where a run starts, is not
#' a case. An *episode* is the clones that became patent in one host on one day, named by
#' `episode_id`, `"<human_id>:<tick>"`, so a bite that delivered two genotypes is two cases of one
#' episode ([episodes()], [recurrences()]). An *acquisition* is a case that began a new blood-stage
#' clone in its host: an import, a local reinfection or a relapse. A recrudescence is a case but not
#' an acquisition, since it re-expands a clone already counted. A relapse into a patent clone of its
#' own lineage is flagged by `into_patent_lineage`, and `exclude_masked = TRUE` leaves such relapses
#' out of a reader's counts.
#'
#' @section COI and MOI:
#' Here the complexity of infection (COI) and the multiplicity of infection (MOI) are one quantity:
#' the distinct patent blood-stage lineages in a host, `moi_active` in `state.tsv`. `MAX_MOI` gates
#' new inoculations: a host already carrying that many blood-stage lineages, patent or subpatent,
#' takes no new bite or import. Relapses bypass it, and a bite that delivers several genotypes is
#' admitted whole, so a host's COI can exceed it.
#'
#' @section Transmission:
#' A *transmission link* joins a parent genome of an oocyst to a genome its products infected:
#' `direct`, through a bite's blood stage, or `liver`, through a relapse of the batch a bite left.
#' [transmission_links()] lists them, each with its *chain*, a connected group of links. *Rc*, the
#' case reproduction number, is the mean number of secondary cases per case, by true origin
#' ([rc_by_origin()]). [classify_chain()] classes each case by its true chain from an import, and
#' [classify_travel()] by what a rule that knows only travel history would conclude.
#'
#' @section Relatedness:
#' Every genome descends from root genomes, which are unrelated: the founder templates, one per
#' founder lineage (`N_FOUNDER_LINEAGES`), and the imports. The founder hosts of one lineage carry
#' clonal copies of its template. So the genealogy is a forest. Identity by descent (IBD) is read
#' from it exactly: [ibd_pairs()] for pairs, [ibd_segments()] for the stretches behind them, and
#' [ibd_lookup()] to read a pair. [arg_tmrca()] gives how far back two genomes share an ancestor,
#' and [arg_summary()] the shape of the forest.
#'
#' @import data.table
#' @importFrom stats setNames
#' @importFrom utils combn
"_PACKAGE"
