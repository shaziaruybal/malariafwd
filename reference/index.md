# Package index

## Start here

The terms used across the package, a run to try them on, and where to
read more.

- [`malariafwd`](https://shaziaruybal.github.io/malariafwd/reference/malariafwd-package.md)
  [`malariafwd-package`](https://shaziaruybal.github.io/malariafwd/reference/malariafwd-package.md)
  : malariafwd: concepts and terms
- [`example_run()`](https://shaziaruybal.github.io/malariafwd/reference/example_run.md)
  : A small run that ships with the package

## Set up

malariafwd runs the external program SLiM, and reads the tree sequence
with Python’s tskit. Check what the package finds before anything else.

- [`check_setup()`](https://shaziaruybal.github.io/malariafwd/reference/check_setup.md)
  : Check what malariafwd needs: R, SLiM, the model, and Python with
  tskit
- [`slim_available()`](https://shaziaruybal.github.io/malariafwd/reference/slim_available.md)
  : Is a SLiM binary available?
- [`slim_model()`](https://shaziaruybal.github.io/malariafwd/reference/slim_model.md)
  : Path to the bundled SLiM model template
- [`tskit_available()`](https://shaziaruybal.github.io/malariafwd/reference/tskit_available.md)
  : Can the tree-sequence functions reach Python with tskit?

## Configure and run

Build a parameter set, run one simulation, or run a grid of parameter
sets and seeds in parallel. Time is in days throughout.

- [`vivax_params()`](https://shaziaruybal.github.io/malariafwd/reference/vivax_params.md)
  [`as.data.frame(`*`<vivax_params>`*`)`](https://shaziaruybal.github.io/malariafwd/reference/vivax_params.md)
  : Build a validated parameter set for the vivax model
- [`param_info()`](https://shaziaruybal.github.io/malariafwd/reference/param_info.md)
  : What each parameter means, its unit and the values it may take
- [`run_sim()`](https://shaziaruybal.github.io/malariafwd/reference/run_sim.md)
  : Run one simulation
- [`run_grid()`](https://shaziaruybal.github.io/malariafwd/reference/run_grid.md)
  : Run a grid of parameter sets x seeds in parallel
- [`load_grid()`](https://shaziaruybal.github.io/malariafwd/reference/load_grid.md)
  : Load the record of a grid

## Set the oocyst and sporozoite draws

LAMBDA_OOCYSTS and LAMBDA_SPOROZOITES are rates before zero-truncation,
not the mean counts a feeding study reports: convert between the two
before setting them.

- [`uptake_moments()`](https://shaziaruybal.github.io/malariafwd/reference/uptake_moments.md)
  : Realised mean and variance of a zero-truncated count distribution
- [`uptake_calibrate()`](https://shaziaruybal.github.io/malariafwd/reference/uptake_calibrate.md)
  : Convert a target realised mean (and variance) into lambda and k

## Read a run’s outputs

The files a run writes: the event log, the daily state of every host and
mosquito, the parameters it was made from, and its tree sequence.

- [`load_run()`](https://shaziaruybal.github.io/malariafwd/reference/load_run.md)
  :

  Load a completed run as a `vivax_run` object

- [`load_events()`](https://shaziaruybal.github.io/malariafwd/reference/load_events.md)
  : Load the event log of a run

- [`load_state()`](https://shaziaruybal.github.io/malariafwd/reference/load_state.md)
  : Load the per-tick state table of a run

- [`load_params()`](https://shaziaruybal.github.io/malariafwd/reference/load_params.md)
  : The parameter set a run was made from

- [`trees_path()`](https://shaziaruybal.github.io/malariafwd/reference/trees_path.md)
  : Path to a run's tree-sequence file

- [`event_types()`](https://shaziaruybal.github.io/malariafwd/reference/event_types.md)
  : The malariafwd event dictionary

- [`malariafwd-outputs`](https://shaziaruybal.github.io/malariafwd/reference/malariafwd-outputs.md)
  [`events.tsv`](https://shaziaruybal.github.io/malariafwd/reference/malariafwd-outputs.md)
  [`state.tsv`](https://shaziaruybal.github.io/malariafwd/reference/malariafwd-outputs.md)
  [`params.tsv`](https://shaziaruybal.github.io/malariafwd/reference/malariafwd-outputs.md)
  [`plasmo.trees`](https://shaziaruybal.github.io/malariafwd/reference/malariafwd-outputs.md)
  : The files a run writes, and their columns

## Check a run

Summary numbers and a dashboard, to see whether a run went extinct,
saturated or settled.

- [`summarise_run()`](https://shaziaruybal.github.io/malariafwd/reference/summarise_run.md)
  : A run's epidemiology over one window of days
- [`plot_diagnostics()`](https://shaziaruybal.github.io/malariafwd/reference/plot_diagnostics.md)
  : One-glance run-health dashboard

## Cases and their origin

Every blood-stage clone and recurrence, labelled by its true origin
(relapse, local reinfection, import or recrudescence), beside the
classes surveillance would give it, and incidence and force of infection
by origin.

- [`episodes()`](https://shaziaruybal.github.io/malariafwd/reference/episodes.md)
  : Every blood-stage clone of a run, from onset to end
- [`recurrences()`](https://shaziaruybal.github.io/malariafwd/reference/recurrences.md)
  : Each recurrence, against the earlier episodes of its host
- [`classify_cases()`](https://shaziaruybal.github.io/malariafwd/reference/classify_cases.md)
  : Classify each blood-stage clone onset ("case") by its true origin
- [`classify_chain()`](https://shaziaruybal.github.io/malariafwd/reference/classify_chain.md)
  : Classify each case by its transmission chain from an import (the WHO
  classes)
- [`classify_travel()`](https://shaziaruybal.github.io/malariafwd/reference/classify_travel.md)
  : Each case as a travel rule classes it: recent travel or not
- [`force_of_infection()`](https://shaziaruybal.github.io/malariafwd/reference/force_of_infection.md)
  : Incidence and force of blood-stage infection, by origin
- [`plot_force_of_infection()`](https://shaziaruybal.github.io/malariafwd/reference/plot_force_of_infection.md)
  : Force of blood-stage infection: molFOB vs FOI

## Hosts over time

Infection states across the population, complexity of infection, and one
host’s history.

- [`host_state_prevalence()`](https://shaziaruybal.github.io/malariafwd/reference/host_state_prevalence.md)
  : Host-state prevalence and the hypnozoite reservoir over time
- [`plot_host_state_prevalence()`](https://shaziaruybal.github.io/malariafwd/reference/plot_host_state_prevalence.md)
  : Host-state prevalence over time (stacked area)
- [`coi_distribution()`](https://shaziaruybal.github.io/malariafwd/reference/coi_distribution.md)
  : Distribution of complexity of infection over the infected host-days
- [`plot_host_trajectory()`](https://shaziaruybal.github.io/malariafwd/reference/plot_host_trajectory.md)
  : Within-host trajectory for one host

## Hypnozoites and relapse

The dormant liver-stage reservoir and the relapses it drives.

- [`relapse_fraction()`](https://shaziaruybal.github.io/malariafwd/reference/relapse_fraction.md)
  : Relapse fraction, by episodes and by acquisitions
- [`dormancy_durations()`](https://shaziaruybal.github.io/malariafwd/reference/dormancy_durations.md)
  : Hypnozoite dormancy durations, from a batch to each of its relapses
- [`time_to_relapse()`](https://shaziaruybal.github.io/malariafwd/reference/time_to_relapse.md)
  : Time from an infection to the first relapse of its hypnozoite batch
- [`plot_relapse_km()`](https://shaziaruybal.github.io/malariafwd/reference/plot_relapse_km.md)
  : Kaplan-Meier curve of the time to first relapse, per hypnozoite
  batch

## Transmission chains

Who infected whom, by a bite or through a relapse, and onward
transmission by origin.

- [`transmission_links()`](https://shaziaruybal.github.io/malariafwd/reference/transmission_links.md)
  : Who infected whom: the ground-truth transmission links
- [`rc_by_origin()`](https://shaziaruybal.github.io/malariafwd/reference/rc_by_origin.md)
  : Case reproduction number (Rc) by true origin

## Relatedness (IBD)

Identity by descent read exactly from the recorded genealogy, not
inferred from markers: pairwise IBD, the segments behind it, and pairs
of known relationship for testing an IBD method. Needs Python with
tskit.

- [`ibd_pairs()`](https://shaziaruybal.github.io/malariafwd/reference/ibd_pairs.md)
  : Ground-truth IBD from a run's tree sequence
- [`ibd_lookup()`](https://shaziaruybal.github.io/malariafwd/reference/ibd_lookup.md)
  [`as.matrix(`*`<vivax_ibd>`*`)`](https://shaziaruybal.github.io/malariafwd/reference/ibd_lookup.md)
  [`as.data.frame(`*`<vivax_ibd>`*`)`](https://shaziaruybal.github.io/malariafwd/reference/ibd_lookup.md)
  : The IBD of pairs of genomes, from ibd_pairs()
- [`ibd_segments()`](https://shaziaruybal.github.io/malariafwd/reference/ibd_segments.md)
  : Pairwise IBD segments along the genome
- [`chrom_map()`](https://shaziaruybal.github.io/malariafwd/reference/chrom_map.md)
  : Chromosome layout of a run's concatenated genome
- [`plot_ibd_pair()`](https://shaziaruybal.github.io/malariafwd/reference/plot_ibd_pair.md)
  : Plot the IBD "karyotype" of one or more pairs
- [`ibd_truth_pairs()`](https://shaziaruybal.github.io/malariafwd/reference/ibd_truth_pairs.md)
  : Pairs of genomes with a known true relationship, and their IBD
- [`plot_relatedness_ladder()`](https://shaziaruybal.github.io/malariafwd/reference/plot_relatedness_ladder.md)
  : Plot the relatedness ladder
- [`new_infection_diversity()`](https://shaziaruybal.github.io/malariafwd/reference/new_infection_diversity.md)
  : Genetic diversity of new infections, over the run or by window

## Genotype data

Marker genotypes painted on the run’s genealogy, and their VCF.

- [`simulate_genotypes()`](https://shaziaruybal.github.io/malariafwd/reference/simulate_genotypes.md)
  [`print(`*`<vivax_genotypes>`*`)`](https://shaziaruybal.github.io/malariafwd/reference/simulate_genotypes.md)
  : Simulate marker genotypes from a run's genealogy
- [`write_vcf()`](https://shaziaruybal.github.io/malariafwd/reference/write_vcf.md)
  : Write simulated genotypes as a VCF file

## Genealogy summaries

The recorded genealogy (the ancestral recombination graph, ARG) is a
forest of unconnected founder lineages: read ?arg_summary first, and use
arg_tmrca() for relatedness.

- [`arg_summary()`](https://shaziaruybal.github.io/malariafwd/reference/arg_summary.md)
  : Structure of the recorded ancestral recombination graph (ARG)
- [`arg_tmrca()`](https://shaziaruybal.github.io/malariafwd/reference/arg_tmrca.md)
  : Pairwise time to the most recent common ancestor, the relatedness
  measure of the ancestral recombination graph
