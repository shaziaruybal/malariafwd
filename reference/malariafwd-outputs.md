# The files a run writes, and their columns

A run's directory holds:

- `events.tsv`, one row per event, read by
  [`load_events()`](https://shaziaruybal.github.io/malariafwd/reference/load_events.md);

- `state.tsv`, one row per host and per mosquito per day, read by
  [`load_state()`](https://shaziaruybal.github.io/malariafwd/reference/load_state.md);

- `params.tsv`, the parameter set and the run's record, read by
  [`load_params()`](https://shaziaruybal.github.io/malariafwd/reference/load_params.md)
  and
  [`load_run()`](https://shaziaruybal.github.io/malariafwd/reference/load_run.md);

- `plasmo.trees`, the tree sequence, found by
  [`trees_path()`](https://shaziaruybal.github.io/malariafwd/reference/trees_path.md)
  and read by
  [`ibd_pairs()`](https://shaziaruybal.github.io/malariafwd/reference/ibd_pairs.md),
  [`ibd_segments()`](https://shaziaruybal.github.io/malariafwd/reference/ibd_segments.md),
  [`arg_summary()`](https://shaziaruybal.github.io/malariafwd/reference/arg_summary.md)
  and
  [`simulate_genotypes()`](https://shaziaruybal.github.io/malariafwd/reference/simulate_genotypes.md);

- `model.slim`, `cmdline.txt`, `run.log` and `DONE`: the model run, the
  command, its log, and the record that the run finished, with the size
  of every file.

## Details

`params.tsv` has a `param` and a `value` column: the record rows
`ACTUAL_SEED`, `SLIM_VERSION` and `MALARIAFWD_VERSION`, then every
parameter with `SEED`, with the derived `GENOME_LENGTH`, `N_CHROMOSOMES`
and `RECOMB_WITHIN`, and last `EXTINCT_TICK` or `NEVER_STARTED` when the
run went extinct or never started. Whole numbers are written as such,
other numbers to 15 significant digits with `.0` on a whole one, flags
as `T` or `F`, and text as given.

`plasmo.trees` holds one node per genome, keyed by its pedigree id, with
node times in ticks before the end of the run; SLiM's second, vacant
node of each genome is not a genome. A run with
`RETAIN_FULL_HISTORY = FALSE` keeps the genomes alive at its end, and
most earlier ones leave it.

## See also

[`load_run()`](https://shaziaruybal.github.io/malariafwd/reference/load_run.md),
[`load_events()`](https://shaziaruybal.github.io/malariafwd/reference/load_events.md),
[`load_state()`](https://shaziaruybal.github.io/malariafwd/reference/load_state.md),
[`load_params()`](https://shaziaruybal.github.io/malariafwd/reference/load_params.md),
[`trees_path()`](https://shaziaruybal.github.io/malariafwd/reference/trees_path.md),
[`event_types()`](https://shaziaruybal.github.io/malariafwd/reference/event_types.md)

## events.tsv

One row per event, in the order the events were logged. -1 means not
applicable and NA no lineage. A run with no event gives no row, and a
column a file lacks reads as NA.

- `tick` (integer):

  the simulation day

- `event_id` (character):

  `"<tick>_<counter>"`, unique

- `event_type` (character):

  the kind of event;
  [`event_types()`](https://shaziaruybal.github.io/malariafwd/reference/event_types.md)
  lists them

- `inoculation_id` (integer):

  the bite or import that brought the parasite, on bite rows
  (`vh_inoculation` and its no-take form) and import rows (an import's
  batch row, logged the next day, carries its primary's id); -1
  elsewhere

- `mosq_id` (integer):

  the mosquito, on uptake, bite and `mosq_death` rows; -1 elsewhere

- `human_id` (integer):

  the host; -1 only on `mosq_death` rows

- `lineage_id` (integer):

  the parasite's lineage, its clonal group; NA on recombinant uptakes,
  on course rows about the host (`bloodcure_given`, `radcure_given`,
  `radcure_failure` and their MDA forms) and on
  `vh_inoculation_import_no_take`

- `parent1_lineage_id` (integer):

  on uptake rows, the lineage of the first gamete; NA elsewhere

- `parent2_lineage_id` (integer):

  on recombinant uptake rows, the lineage of the second gamete; NA
  elsewhere

- `is_recombinant` (integer):

  1 on an outcrossed uptake (`hv_uptake_recombinant`), 0 elsewhere

- `p_outcross` (double):

  on uptake rows, the Simpson diversity of the gametocyte pool the
  mosquito drew from, to 6 significant digits, a diagnostic; 0 elsewhere

- `recipient` (character):

  `"mosquito"` on uptake and `mosq_death` rows, `"human"` on the others

- `is_dormant` (integer):

  1 on a hypnozoite batch's rows (its establishment, deaths and drug
  clearance), 0 on blood-stage and host rows, -1 on no-take rows

- `is_imported` (integer):

  1 when the parasite was acquired by importation, carried within the
  host to its batches, relapses and recrudescences; a local infection by
  an import's offspring is 0

- `pedigree_id` (integer):

  the genome, the key to the tree sequence (pedigree 0 is a real
  genome); -1 on uptake, course and no-take rows

- `parent1_pedigree_id` (integer):

  the genome a row's genome came from: a gamete of an uptake, the
  mosquito product of a bite, the batch of a relapse, the suppressed
  clone of a recrudescence; -1 where none

- `parent2_pedigree_id` (integer):

  the second gamete of a recombinant uptake; -1 elsewhere

- `n_hypnozoites` (integer):

  the hypnozoites a row stands for: a batch's size when established, the
  number activating, dying or cleared, or the reservoir surviving a
  radical-cure failure; 1 on other rows

- `oocyst_id` (integer):

  on an uptake row, the oocyst the mosquito made; on a local bite's
  blood-stage and batch rows, every relapse of that batch and the drug
  clearances of its clones, the oocyst of the infecting product; -1 for
  founders, imports and their batches and relapses, recrudescences, and
  rows that do not record it (`human_clearance`, `hypnozoite_clearance`,
  `mosq_death`, course rows)

- `detected` (integer):

  on the blood-stage onset rows of bites, imports, relapses and
  recrudescences, 1 if the observation layer reported the case and 0 if
  not; -1 elsewhere, founder seeding included

- `into_patent_lineage` (integer):

  on `hypnozoite_activation` rows, 1 if a patent clone of the same
  lineage, born on an earlier day, was still in the host (one that ends
  that day counts), 0 if not; -1 elsewhere

## state.tsv

One row per host and per mosquito per day.
[`load_state()`](https://shaziaruybal.github.io/malariafwd/reference/load_state.md)
gives the host rows, and `entity = "mosquito"` or `"all"` the
mosquitoes' or both.

- `tick` (integer):

  the simulation day

- `entity_type` (character):

  `"human"` or `"mosquito"`

- `entity_id` (integer):

  the host's or the mosquito's id

- `tag` (integer):

  hosts: 1 if the host carries a patent blood-stage clone, or a
  subpatent one when `SUBPATENT_INFECTIVITY > 0`, whether or not it
  carries gametocytes; 0 otherwise, hypnozoites alone included;
  mosquitoes: 1 if it carries a parasite product, inside or past the EIP

- `moi_active` (integer):

  hosts: the distinct patent blood-stage lineages (the COI), suppressed
  clones excluded; mosquitoes: the distinct lineages it carries, inside
  the EIP included

- `moi_subpatent` (integer):

  hosts: the distinct lineages suppressed below patency, present but
  undetectable; 0 unless `P_RECRUD_GIVEN_FAIL > 0`

- `moi_dormant` (integer):

  hosts: the distinct genotypes in the hypnozoite reservoir

- `hyp_reservoir_size` (integer):

  hosts: the hypnozoites, the sum of the batches' counts

- `lineage_ids` (character):

  the lineages present, `|`-separated and de-duplicated, even for one;
  NA when none

- `n_transmissible` (integer):

  mosquitoes: the distinct lineages past the EIP; NA for hosts

- `mosq_age_days` (integer):

  mosquitoes: the age in days; NA for hosts

- `infectious_bites` (integer):

  hosts: the bites received that day from a mosquito carrying a product
  past the EIP, whatever the bite then did, the EIR's numerator; NA for
  mosquitoes
