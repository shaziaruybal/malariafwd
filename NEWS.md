# malariafwd 0.1.0

The first release. malariafwd runs a forward, individual-based simulation of *Plasmodium vivax*
transmission and genetics from R, and reads its recorded genealogy as ground truth: who infected
whom, whether each recurrence is a relapse, a reinfection or a recrudescence, and how related any
two parasite genomes are. It is a pre-release: the interface and the model may change before
1.0.0. It needs SLiM 5.2 or newer, and Python with tskit for the relatedness functions;
`check_setup()` reports both.

## What it contains

* **The model**, in SLiM with tskit tree-sequence recording: day-by-day human-mosquito
  transmission with superinfection; a hypnozoite reservoir whose batches activate (relapse) or
  die; case management, radical cure, mass drug administration and recrudescence; four-strand
  meiosis in the mosquito on the 14 chromosomes of the PvP01 reference genome; importation,
  vector control, seasonality and resurgence. Everything beyond the core is off by default.
* **Running it**: `vivax_params()`, with every parameter described in `?vivax_params` and
  `param_info()`; `run_sim()`; `run_grid()`, for parallel sweeps that resume where they stopped.
* **An example run**: `example_run()` returns a small run that ships with the package, so every
  function that reads a run can be tried without SLiM.
* **Ground truth**: `classify_cases()`, `episodes()`, `recurrences()`, `transmission_links()`,
  `classify_chain()` and `rc_by_origin()`, with `classify_travel()` for what a travel-history
  rule concludes instead.
* **Relatedness**: `ibd_pairs()` and `ibd_segments()` read identity-by-descent from the tree
  sequence; `simulate_genotypes()` and `write_vcf()` paint marker genotypes on the genealogy, with
  the truth behind every call; `arg_tmrca()` and `arg_summary()` describe the genealogy itself.
* **Summaries and plots**: `summarise_run()`, `force_of_infection()`, `relapse_fraction()`,
  `new_infection_diversity()`, `plot_diagnostics()` and the other `plot_*()` functions.
* **Validation**: the scripts in `inst/validation`, installed under
  `system.file("validation", package = "malariafwd")`, audit the event log, the host states, the
  tree sequence and the IBD against what the model must produce, and check that each parameter
  moves the outputs the way it should.

## Before drawing conclusions

* `NU_HYPNOZOITE_BATCH = 3` is fitted to four calibration targets, ranges taken from published
  *P. vivax* studies, at moderate transmission (`inst/validation/calibrate_nu.R` repeats the fit). The
  other defaults are an illustrative preview, 40 hosts for a year, not a calibrated setting.
* Imported parasites are unrelated to local ones, so an import-detection result is an upper bound.
* No study has measured how many genotypes a natural bite delivers. At the defaults
  (`LAMBDA_OOCYSTS` and `LAMBDA_SPOROZOITES` of 2) a bite can co-transmit several: sweep them.
* There are no within-host densities, no post-treatment prophylaxis and no spatial structure.

The README lists these limitations with the others. The
[model description](https://shaziaruybal.github.io/malariafwd/MODEL_DESCRIPTION.html) gives the
model in full, with its assumptions and its validation, and the
[design decisions](https://shaziaruybal.github.io/malariafwd/DECISIONS.html) give the reasons
for each modelling choice.

## In development

* **Each hypnozoite genotyped independently of the bite's primary infections**, so that a relapse
  can differ from the primary of its own bite. This will change every default run.
* **A structured source population for imports**, so that imported parasites carry real ancestry
  and can be related to one another.
* **Post-treatment prophylaxis**: protection that wanes after each course.
* **More extensive calibration of the treatment and recurrence machinery**, including a
  chloroquine-resistant setting, in which recrudescence is material rather than negligible.

The roadmap lists every planned extension with its motivation.
