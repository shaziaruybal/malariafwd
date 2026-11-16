# Design decisions

Why the model and the R package are built the way they are. Each entry states a choice, the reason
for it, the alternatives considered and set aside, and the tests or validation scripts that check
it. The [model description](https://shaziaruybal.github.io/malariafwd/MODEL_DESCRIPTION.html) says
what the model does; this page says why.

The entries are grouped by topic and numbered D1 to D54. The model file
(`inst/slim/malariafwd_vivax.slim`), the R code and the model description cite them by number, as in
"(D21)". *How to check* names tests in `tests/testthat/`, which `devtools::test()` runs, and
validation scripts in `inst/validation/`, whose last results are in `inst/validation/RESULTS.md`.

- The simulator, a run and its record (D1–D8)
- Hosts, blood stage and superinfection (D9–D12)
- The mosquito stage, meiosis and the genetics it makes (D13–D20)
- Hypnozoites, relapse and the reservoir (D21–D28)
- Detection, treatment, recrudescence and MDA (D29–D35)
- Seasonality, vector control and importation (D36–D41)
- Parameters, outputs and the R interface (D42–D54)

## The simulator, a run and its record

What a run is and how it is made: the model in SLiM, one day a tick, launched from R as a
subprocess, with a record that repeats it, the outcomes it writes down, and a grid that resumes from
what is on disk.

### D1 — The package is the engine and the readers of its ground truth

**Decision.** malariafwd is the SLiM model, the R functions that run it one simulation or a whole
grid at a time, and the readers of what a run records: loaders, epidemiological and genetic
summaries, ground-truth IBD, marker genotypes painted on the genealogy, and plots, one function to a
question (D54). Research analyses built on it sit outside the package. The package and its model
share the name (`malariafwd_vivax.slim`), which leaves room for species other than *P. vivax*.

**Why.** The package holds what every use of the model needs, the engine and the readers of its
truth, and nothing that belongs to one study, so it does not move with research in progress.

**Alternatives considered.** The engine alone would be too thin to use: every user would rebuild the
IBD and the metrics a run is read with. The engine with a study's analysis framework would tie the
package to that research.

**How to check.** `tests/testthat/test-exports.R` lists the exports exactly.

### D2 — SLiM runs as a subprocess, and Python is reached through one check

**Decision.** `run_sim()` starts `slim` as a process through processx, with no shell, and gives it
every parameter as a `-d NAME=value` constant; the model writes its own output files, and SLiM's
console output goes whole to `run.log`. The binary is the one `MALARIAFWD_SLIM` names, or `slim` on
the `PATH`, and must be SLiM 5.2 or newer: the model is tested on 5.2, and 5.0 and 5.1 give
different output for the same seed. The same model file runs without R, with `slim -d`, and its
`initialize()` then makes the checks R makes (D43). Its per-event narration prints only under
`VERBOSE = TRUE`: a quiet log holds little beyond SLiM's start-up, the genome line and the closing
banner, and a run's events, states and genealogy are the same either way. Python serves only to read
tree sequences, and reticulate is suggested rather than required, so running the model and reading a
run's tables need no Python. The tree-sequence functions reach tskit (>= 1.0) and numpy through
reticulate behind one check (`R/python.R`), which declares the modules before Python starts: a
Python the user chose is used as it is, and with none chosen reticulate builds its own environment
with them. A check that finds everything holds for the session. A function that cannot reach the
modules stops, naming itself, with class `malariafwd_no_python` and one paragraph: the Python in use
and how it came to be used, what is missing, and what to do. `tskit_available()` answers without a
condition, and `check_setup()` reports SLiM, the model and the Python side.

**Why.** Launching a process costs little beside a run, which takes tens of seconds for 150 hosts
over three years (MODEL_DESCRIPTION 16.6), and a model that writes its own files needs nothing from
R while it runs. Without a shell, the quotes of an Eidos string, such as `OUTDIR='.'` or a `SIM_ID`
with an apostrophe, reach SLiM as written. One model file cannot drift from a quiet copy of itself;
the quiet default keeps a sweep's logs small, and a single run can still be narrated to debug it.
Declaring the modules changes nothing for a Python the user chose and is what reticulate needs to
build its own, and one check gives every function the same modules and the same message.

**Alternatives considered.** slimr adds a layer that parses SLiM's output. `system2()` passes the
arguments through a shell, which strips the quotes of an Eidos string. A separate quiet model file
would drift from the narrated one, and removing the narration would lose the debugging mode.
Declaring the modules when the package loads would add tskit to reticulate's environment for a user
who never reads a tree sequence, and `py_module_available()` alone is `FALSE` until Python has
started.

**How to check.** `tests/testthat/test-run.R` (an apostrophe and a space in the output path and the
`SIM_ID`, and a non-ASCII `SIM_ID` in a C locale), `tests/testthat/test-slim.R` (a SLiM older than
5.2 is refused), `tests/testthat/test-python.R` (every branch of the Python check, and
`check_setup()`) and `tests/testthat/test-run-dir.R` (a run's key ignores `VERBOSE`).

### D3 — Time advances in daily ticks, and a duration parameter is the mean of a geometric duration

**Decision.** One tick is one day, in SLiM's nonWF order. `reproduction()` draws the day's
hypnozoite activations, as relapses, and deaths; clones the hypnozoite batch of an import that
arrived the day before; runs the bites, with uptake and meiosis in the mosquito and inoculation in
the host; and re-expands the suppressed clones whose day has come, as recrudescences. `early()`
removes the emptied hypnozoite batches and the clones the day's recrudescences replaced, which SLiM
cannot kill in `reproduction()`; gives the courses that are due; clears blood-stage clones; kills
and replaces mosquitoes; gives an MDA round on `MDA_TICK`; brings the day's imports; and ends the
run if no parasite is left and none can still come (D5). `late()` writes the day's rows of
`state.tsv`, and on the last day the tree sequence. A duration parameter D, such as
`MEAN_BLOOD_STAGE_DURATION` or `MOSQ_LIFESPAN`, is the mean of a geometric duration: the event
happens each day with probability 1/D. `TREATMENT_DELAY` and, by default, `MEAN_RECRUD_DELAY` follow
the same rule from the day after the case or the suppression, so each lasts at least a day (D32).
`BITING_RATE`, `ALPHA_ACTIVATE`, `MU_HYPNOZOITE` and `DELTA_LATENCY` are daily probabilities, so the
waiting times they set are geometric too; `EIP_DAYS` and `GAMETOCYTE_DELAY` are delays of a set
number of days. A mosquito bites at most once a day, so `vivax_params()` warns when `BITING_RATE`
times its seasonal, intervention and resurgence multipliers passes 1 on a day of the run.

**Why.** Every duration has one reading, and each parameter means what it says: at the default of 60
days, blood-stage clearance is White et al. (2014)'s rate of 1/60 a day applied as a daily
probability, and an infection lasts 60 days on average. Two steps sit where SLiM requires: imports
arrive in `early()`, since `reproduction()` runs only while a parasite exists and an import must be
able to seed an empty population (D38), and courses are given in `early()`, the first point at which
a genome made in that day's `reproduction()` can be seen, so radical cure reaches the batch the
case's own bite seeded (D31).

**Alternatives considered.** Converting each duration with 1 - exp(-1/D), natural in a
continuous-time model, would not fit a model that steps in days and whose other processes are daily
probabilities already, and it would give a parameter of D days a mean of about D + 1/2 (2.54 days at
D = 2).

**How to check.** `tests/testthat/test-run.R` (founder infections with
`MEAN_BLOOD_STAGE_DURATION = 2` last 2 days on average), `tests/testthat/test-treatment.R` (the
treatment and recrudescence delays have mean D and last at least a day) and
`tests/testthat/test-params.R` (the biting warning follows the model's own seasonal, intervention
and boost timeline).

### D4 — Every genome is remembered in the tree sequence before it is killed

**Decision.** Every path that removes a genome first calls
`treeSeqRememberIndividuals(permanent = RETAIN_FULL_HISTORY)`: natural blood-stage clearance; a case
course or an MDA round that clears a clone, a subpatent one included; the removal of a hypnozoite
batch emptied by activation and death, by radical cure or by MDA; the death of a mosquito, with the
products it carries; and the removal of a clone that a recrudescence has replaced.
`RETAIN_FULL_HISTORY` is `TRUE` by default, so `plasmo.trees` keeps every genome a run made. With
`FALSE` the genomes are remembered without the permanent flag and most earlier ones leave the tree
sequence: `ibd_pairs()` lists a requested id it does not hold in `missing`, warns once and looks it
up as `NA`, and asking only for such ids is an error.

**Why.** IBD between a recurrence and an infection that has already cleared, which `recurrences()`
reports and every relapse-to-origin comparison rests on, needs the cleared genome's ancestry, and a
genome killed without being remembered can be simplified out of the tree, its ancestry with it.
`FALSE` is for runs that will never ask about earlier infections.

**How to check.** `inst/validation/audit_events_and_ibd.R` (the logged cases are in the tree
sequence), `tests/testthat/test-run.R` (every relapse is found, at IBD 1, with the batch it woke
from) and `tests/testthat/test-no-parasite.R` (with `RETAIN_FULL_HISTORY = FALSE`, the readers give
their `NA` rows, and the error names the setting).

### D5 — A run that dies out, and one that never starts, are results the run writes down

**Decision.** A run goes extinct on the first tick on which no parasite is left in any host or
mosquito, hypnozoites included, and no importation is still to come. It ends on that tick, and
records it even when that tick is its last day: it writes its tree sequence, its `events.tsv` and
`state.tsv` stop there, and its `params.tsv` ends with an `EXTINCT_TICK` row, which `load_run()`
reads into `extinct_tick` (`NA` for a run that did not go extinct) and `run_sim()` reports. A
founder-free run whose imports all fail to establish a parasite, by chance (a Poisson rate that
draws none, imports that take in neither compartment) or because its settings let none take
(`P_PRIMARY = 0` with a hypnozoite batch mean of 0, an import rate of 0 throughout its window),
never starts: it ends on the first tick after its last import day but not before tick 10, and in any
case on its last day, writes its outputs as an extinct run does, and ends `params.tsv` with a
`NEVER_STARTED` row, which `load_run()` reads as `never_started = TRUE`. `run_grid()` gives each
outcome a column and counts its cells in a message of its own; a cell that never started is `ok`. An
import is still to come until the run passes `LAST_IMPORT_DAY`, the last day of the run on which one
can arrive (0 if none): a window inside the run counts to its last day even while its rate is 0, and
a pulse or window outside the run brings nothing. The model's start-up check and its end-of-day
check read that one constant, and R works out the same days. A founder-free set with no import path
inside the run is refused before SLiM starts (D42). Every reader takes a run with no parasite: its
tables have their usual columns and no rows, its counts are 0, and a quantity that needs a parasite
is `NA`; the summaries of a window count the days after an extinction as days with no parasite
(D49).

**Why.** Fade-out and a failure of importation are results: how often imports fail to seed an
outbreak is itself worth counting. Every later tick would be empty, so the run ends there. Recording
an extinction on the last day too means `NA` says only that a run did not go extinct, which a grid
counting extinctions needs, and the two outcomes are facts of the run rather than parameters, so
every run's parameter set has the same names. One constant for the import days keeps the checks from
disagreeing about whether an import can still come.

**Alternatives considered.** Stopping a run that never started with an error would lose a sweep's
cell to chance. Recording it as an extinction would confuse a run that never had a parasite with one
whose parasites died out. Refusing a founder-free set whose imports cannot take would remove a
legitimate null arm of a grid, as a closed border is, whose honest outcome is that it never started.

**How to check.** `tests/testthat/test-run.R` (a run that dies out, also with imports set outside
it; founder-free runs whose imports never take, and the tick each import schedule ends on),
`tests/testthat/test-params.R` and `tests/testthat/test-slim.R` (a founder-free set with no import
path inside the run is refused, by R and by the model run directly), `tests/testthat/test-load.R`
and `tests/testthat/test-no-parasite.R` (extinct and never-started runs load with the usual columns,
and every reader takes them).

### D6 — A run directory holds one whole run or what it held before

**Decision.** `run_sim()` writes a run into a new directory beside `outdir`, starts SLiM there
(`OUTDIR` is `'.'`) on the run's own copy of the model, `model.slim`, and moves the directory to
`outdir` only once SLiM has exited 0, the four outputs exist and `DONE` has been written, last. A
run that fails, times out or is interrupted goes to `<outdir>.failed` with its `run.log`, which
holds SLiM's error output too; SLiM is stopped if R dies, which leaves an `<outdir>.partial-*`
directory to delete by hand. `DONE` records when the run finished, how long it took, its seed, a key
of its parameters and seed (the md5 of the sorted arguments SLiM is given, without `OUTDIR` and
`VERBOSE`), the md5 of the model, the versions and the size of every file. A directory without
`DONE`, or whose files have changed since, loads with `complete = FALSE`, and every reader given it
warns once per directory per session (class `malariafwd_incomplete_run`). An `outdir` that holds
this run is returned with `skip_existing = TRUE`; one that holds any run, finished or not, is
replaced with `overwrite = TRUE` once the new run is in place, and is otherwise refused (class
`malariafwd_outdir_in_use`); one that holds other files is never used, and only a run's own files
are ever deleted.

**Why.** A run directory is ground truth: it should hold one whole run or none, say which, and not
be replaced by accident. Writing beside `outdir` makes the final move a rename within one file
system.

**Alternatives considered.** Deleting `outdir` first would let a failed run destroy the one before
it. Writing in place with `DONE` alone would still leave a failed run's files mixed with those of
the run before it. Loaders that refuse a directory without `DONE` would lock out runs started with
`slim` directly, which have none. Hashing every output into `DONE` would cost time, where the sizes
catch a truncated copy at no cost. A marker written by the model cannot know that R holds the whole
run. Replacing a run by default would let a script run again replace ground truth without a word.

**How to check.** `tests/testthat/test-run-dir.R` (a run stopped part-way leaves `outdir` as it is;
`outdir` appears only when the run has finished; SLiM's error output reaches `run.log`; a run is
reused or replaced only when asked; other files are never touched; `DONE`'s rows and key; an
incomplete directory warns once).

### D7 — A run's record rebuilds the run, and one seed gives one set of tables

**Decision.** `params.tsv` holds every resolved parameter, with `ACTUAL_SEED`, the SLiM and package
versions, the derived `GENOME_LENGTH`, `N_CHROMOSOMES` and `RECOMB_WITHIN`, and the outcome rows of
D5; floats are written to the 15 significant digits `run_sim()` passes SLiM, a whole one with its
decimal point (`10.0`) and an exponent form as printf writes it (`7.4e-07`). `vivax_params()` stores
each value in its default's type and as SLiM is given it: a count or a scheduling day as an integer,
and a double as the value of its 15-digit text. `load_params()` reads the file back into the
`vivax_params()` set that made the run, every parameter of its type and none validated, so
`run_sim(load_params(run), seed = run$seed, outdir = tempfile())` repeats a run started from R
whenever each of its non-zero floats lies between 1e-130 and 1e15 in magnitude. The seed, the
versions, the derived genome and the outcomes are typed fields of the `load_run()` object instead. A
`params.tsv` that lacks a parameter, or holds a row that is neither a parameter nor a record row,
loads with the default or without the row and warns once per directory per session (class
`malariafwd_params_mismatch`). `run_sim()` gives SLiM the seed both as `-s` and as `SEED`, which the
model sets, so the tree sequence's provenance and `params.tsv` agree, and `cmdline.txt` holds SLiM's
arguments, so `slim` run with them beside `model.slim` repeats the run. Every genome's x, y and z
are set to 0 where it is made, in `setHostGenome()` or `tagMosqProduct()`, so two runs of one seed
give the same tree-sequence tables, all but the provenance row, where SLiM records the run's date,
its duration and the memory it used.

**Why.** A run's record should be the run, and a parameter set should hold what the run is given;
the record rows describe a run rather than set it, so they belong to the run object. A default
filled in for a missing parameter is not a setting of that run, hence the warning. SLiM writes each
individual's coordinates into the tree sequence even in a non-spatial model and leaves them
undefined until a model sets them, and a new individual can take over the memory of one that died,
so unset coordinates would make two runs of one seed differ in their tables while their `events.tsv`
and `state.tsv` agree.

**Alternatives considered.** Seventeen digits, or hexadecimal floats, would carry any double
exactly, but Eidos cannot read hexadecimal floats, and seventeen digits print values with digits
nobody set (0.1 as 0.10000000000000001); at 15, the set, its record and what SLiM is given still
hold one value. The record rows in an attribute of the set would be lost to most list operations.
Validating on load would let a file's content stop a reader; `run_sim()` checks a set before it
runs. Setting the coordinates once a day would record a genome made and remembered on one day before
they are set.

**How to check.** `tests/testthat/test-params.R` (a double is stored as SLiM is given it, each value
has its default's type, a rebuilt set is the set), `tests/testthat/test-slim.R` (`params.tsv`
records floats to 15 significant digits, from R and from a direct run), `tests/testthat/test-run.R`
(the set read from a run is the set that made it, a run from it gives the same files, and the record
is typed), `tests/testthat/test-load.R` (types, and a file that does not match warns once),
`tests/testthat/test-run-dir.R` (the directory repeats its own run),
`tests/testthat/test-model-invariants.R` (coordinates are 0 where SLiM reuses a dead genome's
memory, and two runs of one seed give the same tables) and `inst/validation/check_parameters.R` (one
seed gives a byte-identical event log).

### D8 — A grid checks every cell first, runs one future per cell, resumes from disk and records itself

**Decision.** `run_grid()` rebuilds every parameter set and checks every cell's directory against
what `run_sim()` would do with it before any cell starts: a refused set, or a cell that holds
another run, one that did not finish or files that are not a run's, stops the grid with them named.
A set's name names its cells' directories (`<set>_seed<seed>`), so names must differ, also ignoring
case, and be usable as directory names; seeds must be distinct whole numbers. Each cell is one
`future::future()`, labelled with the cell so that an error from a lost worker names it, and runs
under the caller's plan or, with `workers`, a multisession plan for the call; it takes no seed,
since the cells draw no random numbers in R, and the caller's random-number stream is left as it is.
A cell that already holds its run is not run again (`skip_existing = TRUE`), so a second call
resumes a grid that stopped part-way, and `overwrite = TRUE` runs the cells of an edited set again.
A cell that fails does not stop the grid: its directory is left as it is, its files go to
`<cell>.failed`, and the grid warns, with the reason in the `error` column. The calling session
writes `grid.tsv` at the end, with paths relative to the grid's directory: each cell's outcome,
versions and key, and every parameter as SLiM is given it, which `load_grid()` reads. The result is
a `data.table`, one row per cell, with the parameters that vary between the sets.

**Why.** A grid is launched, left and launched again after a failure: it should cost only the cells
that are missing, say what it holds, and not run half its cells before saying what is wrong. A
future plan runs the same call on a laptop's cores or a cluster's nodes, and cells run by other
means, such as one `run_sim()` per job of a job array, are recorded by calling `run_grid()` with the
same sets and seeds afterwards, which runs nothing.

**Alternatives considered.** A grid class would misreport subsets, since data.table keeps a subclass
through subsetting. Keying a cell on the SLiM version would let an upgrade invalidate every finished
grid; the versions are recorded instead. Workers writing `grid.tsv` would race on one file.
`parallel::mclapply()` forks on Unix only and has no path to a cluster. furrr or future.apply would
each add a dependency for what one future per cell does.

**How to check.** `tests/testthat/test-grid.R` (names, seeds, the caller's stream, resuming, an
edited set, `grid.tsv` and `load_grid()` across a moved directory, a failed cell), the grid tests in
`tests/testthat/test-run.R`, and part 3 of `inst/validation/check_parameters.R`, which runs a grid
with two workers on the installed package.

## Hosts, blood stage and superinfection

The host side: how long a clone lasts, how many a host can take, which of its parasites a mosquito
can see, and when a clone can pass parasites on.

### D9 — A blood-stage clone clears on a daily draw, with no age ceiling

**Decision.** Each patent blood-stage clone clears on its own each day, in `early()`, with
probability `1/MEAN_BLOOD_STAGE_DURATION`, so a host carrying several clones loses them
independently and its MOI decays rather than collapsing at once. A clone cannot clear on the day its
blood stage began, so its blood stage lasts a geometric number of days, at least one, with mean
`MEAN_BLOOD_STAGE_DURATION` (D3). No ceiling clears a clone by age. A subpatent clone makes no draw,
and a recrudescence is a new clone (D33); since no draw depends on age, a recrudescence clears as a
new infection does, and a clone that cycles through suppression has no limit on its total
persistence. `MEAN_BLOOD_STAGE_DURATION` must be finite and at least 1, in `vivax_params()` and in
the model's `initialize()`.

**Why.** The default of 60 days is White et al. (2014)'s clearance rate, 1/60 a day, used as a daily
probability. A ceiling would cut off the geometric tail, so the realised mean would fall short of
`MEAN_BLOOD_STAGE_DURATION` and the parameter would mean something other than what it says. Nor
would it guard against anything: the only means that stop clearance, a negative or an infinite one,
are refused by the requirement that the mean be finite and at least 1, which the daily probability
needs anyway, since below 1 it would exceed 1.

**Alternatives considered.** A high ceiling as a safety net, 365 days say: it would still cut the
tail wherever `MEAN_BLOOD_STAGE_DURATION` is long, and guard against nothing the parameter check
does not already refuse.

**How to check.** `tests/testthat/test-run.R` checks, with founders only, that clones last a
geometric time with the parameter's mean and that none is forced to clear at a set age.
`tests/testthat/test-params.R` refuses a mean below 1 or an infinite one, and
`tests/testthat/test-slim.R` a mean below 1 in a run started directly with `slim -d`.

### D10 — One superinfection gate for every new inoculation, counting the clones a host gained that day

**Decision.** A host already carrying `MAX_MOI` blood-stage lineages, patent or subpatent, takes no
new inoculation, a local bite or an import alike: the whole inoculation is skipped, blood-stage
primary and hypnozoite batch, so a full host's reservoir is not topped up either. A bite is checked
whatever the host's infectious state, before its `P_TRANS_VH` draw, so a refused bite makes no such
draw; the refusal stops only the host's inoculation, and on the same bite an infectious host can
still infect the mosquito. Each host keeps a list of the lineages of the blood-stage clones it
gained earlier in the day, by a relapse or a bite, emptied at the start of each day's
`reproduction()`, and a bite's check counts them with the clones the host carries; an import arrives
in `early()`, after the day's bites, so its check finds them among the host's clones. A bite on a
host whose only clones arrived that day is inoculated like any other: only the mosquito's uptake
needs a clone the host already shows, and with none it is skipped without a draw. Relapses bypass
the gate, so realised MOI can exceed `MAX_MOI`; so can a bite that delivers several genotypes, which
is admitted whole, and a founder host whose genome is a copy of another founder's, since the copy
arrives on day 2 after that day's bites and so counts against neither a day-1 import (D40) nor a
day-2 bite. A refused import is logged as `vh_inoculation_import_no_take`; a refused bite leaves no
row. `MAX_MOI` is at least 1, in `vivax_params()` and in the model's `initialize()`.

**Why.** `MAX_MOI` is about parasite burden, not about what a microscope or a mosquito can see, so
subpatent clones count (D11), and it is one rule for every inoculation, whatever its route. The
order of the bites within a day is an artefact of the loop and should decide neither whether a bite
counts nor whether the cap holds; SLiM adds a clone made in `reproduction()` to its host only when
that callback ends, so without the day's list the clones a host gained earlier that day would be
hidden from its later bites. `MAX_MOI = 0` would refuse every inoculation, an empty host's included,
where 1 already means no superinfection. At the default of 20 the gate is a guard that seldom binds,
COI being typically 1 to 3; set low, it becomes the quantity a run measures, as section 14 of the
model description warns.

**Alternatives considered.** Gating a bite only when the host has an infectious clone: two gates for
one rule about parasite burden, under which, at `SUBPATENT_INFECTIVITY = 0`, a host whose clones are
all subpatent would take bites past the cap while refusing imports. Counting only the clones a host
shows: a second bite on the same day would get a second clone in even at `MAX_MOI = 1`.

**How to check.** `tests/testthat/test-model-invariants.R`: at `MAX_MOI = 2`, with no relapses and
one genotype per bite, no host ever carries more than 2 blood-stage lineages; at `P_TRANS_VH = 1`,
with case management and no host full, each host's infectious bites a day are its inoculations that
day. `tests/testthat/test-run.R`: at `MAX_MOI = 1` no import lands in a host that still carries its
clone, and the refusals are logged; at `SUBPATENT_INFECTIVITY = 0`, with every course failing into
subpatency, no bite reaches a host once its clone is suppressed. `tests/testthat/test-params.R` and
`tests/testthat/test-slim.R` refuse `MAX_MOI = 0`, the second in a run started directly with
`slim -d`. `inst/validation/check_parameters.R` checks that `MAX_MOI = 2` lowers mean COI.

### D11 — Patency, infectiousness and capacity are three separate counts

**Decision.** The model asks three questions of a host's blood stage and answers each with its own
count of distinct lineages; none counts hypnozoite batches. `getActiveLineageCount()` counts the
patent lineages: the COI, `moi_active` in `state.tsv`, which leaves out a subpatent clone, since it
is undetectable and in no episode. `getTransmissibleLineageCount()` counts the lineages that can
infect a mosquito, subpatent ones included whenever `SUBPATENT_INFECTIVITY > 0`; it sets the host's
infectious tag, the `tag` column of `state.tsv`, and is the test the bite loop makes before a
mosquito can take up parasites. Which of those clones a mosquito draws from is settled at each
uptake: each subpatent clone joins with probability `SUBPATENT_INFECTIVITY`, and the gametocyte gate
applies (D12). `getBloodStageLineageCount()` counts every blood-stage lineage present, patent or
subpatent, which is the `MAX_MOI` question for an import; for a bite, `getCapacityLineageCount()`
adds the lineages the host gained earlier that day (D10). Each call site asks the question it means.
At the shipped defaults no clone is ever subpatent and the three agree; they part once a failed
course can leave a clone suppressed (D33).

**Why.** Subpatency separates three properties that coincide without it. A suppressed clone is
undetectable, so it is in no episode and no part of the COI. It can still infect a mosquito, if less
well than a patent clone: in the feeding studies Koepfli et al. (2021) compile, microscopy-positive
carriers are 5 to 10 times as infective as submicroscopic ones. And it still occupies its host,
since `MAX_MOI` is about parasite burden (D10).

**Alternatives considered.** The patent count for all three questions, or a separate infectious tag
with the bite loop still testing the patent count. Either way the bite loop passes over a host whose
clones are all subpatent, so that host never infects a mosquito whatever `SUBPATENT_INFECTIVITY`
says, and nothing shows it: the subpatent compartment is usually small, so a missing route looks
like a minor one. The patent count would also let a host whose clones are subpatent take bites past
`MAX_MOI`.

**How to check.** `tests/testthat/test-run.R` refuses a bite on a host at `MAX_MOI = 1` whose clone
is suppressed at `SUBPATENT_INFECTIVITY = 0`, a host the infectiousness count leaves out;
`tests/testthat/test-edge-values.R` checks that a host with only hypnozoite batches is not
infectious, and `tests/testthat/test-load.R` that every host tagged infectious carries a patent or
subpatent clone. The uptake route itself has no test of its own; to check it, join each
`hv_uptake_*` row on day T to its host's `state.tsv` row on day T − 1, in a run whose courses fail
into subpatency: uptakes from hosts with `moi_active = 0` number 0 at `SUBPATENT_INFECTIVITY = 0`
and rise with it.

### D12 — Gametocyte gating is in the model and acts as the identity at the shipped defaults

**Decision.** A blood-stage clone can infect a mosquito only from `gametocyte_ready_tick`, its onset
plus `GAMETOCYTE_DELAY` days, and only if it carries gametocytes, which a Bernoulli draw at
`P_GAMETOCYTE` decides once for the clone and keeps. Every new clone, a relapse's and a
recrudescence's included, starts its own delay and makes its own draw. A clone first meets the bite
loop the day after its onset, since a clone made in `reproduction()` joins its host only when that
callback ends and an import arrives in `early()`, after the bites; so a `GAMETOCYTE_DELAY` of 0 acts
as 1. A mosquito that bites a host none of whose clones pass the gate takes up nothing. The gate
applies after the subpatent clones have been thinned at `SUBPATENT_INFECTIVITY` (D11); at the
shipped `GAMETOCYTE_DELAY = 0` and `P_GAMETOCYTE = 1` it passes every clone and draws no random
number, so a default run is identical to one without the gate.

**Why.** Vivax gametocytes appear early, at or near the onset of symptoms rather than after the 10
to 12 days falciparum gametocytes take to mature (Bousema & Drakeley 2011), so infectiousness runs
close to the blood stage, and a realistic setting, a delay of 2 to 3 days with carriage of about
0.9, lies near the identity. The parameters are there for the sensitivity analysis that section 14
of the model description calls for. The delay is one value for every clone: at a few days against a
60-day mean blood stage (D9), its variance can do little to the dynamics.

**Alternatives considered.** Removing the two parameters: they are the knob a sensitivity analysis
of gametocyte timing and carriage needs. Shipping the realistic setting as the default: it would
move the transmission at which `NU_HYPNOZOITE_BATCH` is fitted (D21) for little gain, the setting
being close to the identity. A delay drawn for each infection, from a Poisson or a uniform over 0 to
3 days: simple to add, but it would matter only to a question about the delay's tail, such as the
share of early transmitters or how the timing of treatment meets infectiousness.

**How to check.** `inst/validation/check_parameters.R` runs a configuration with
`GAMETOCYTE_DELAY = 3` and `P_GAMETOCYTE = 0.8` through its structural checks; no check measures the
gate's effect on uptake.

## The mosquito stage, meiosis and the genetics it makes

Where the genetics come from: how much a mosquito takes up, what meiosis does to it, the genome it
happens on, and how relatedness is read back out of the recorded genealogy.

### D13 — Uptake and inoculation counts are zero-truncated, and their lambdas are pre-truncation rates

**Decision.** A successful blood meal forms a zero-truncated Poisson(`LAMBDA_OOCYSTS`) number of
oocysts, each its own meiosis (D14), negative binomial with variance λ + λ²/K when `K_OOCYSTS > 0`.
An infectious bite samples a zero-truncated `LAMBDA_SPOROZOITES` number of a mosquito's post-EIP
products, with replacement and weighted by `spz_weight`, and repeat draws of one lineage collapse to
a single genotype. A lambda of 0 gives exactly one. Both default to 2, so an uptake forms 2.31
oocysts and a bite samples 2.31 products on average, and a bite delivers several genotypes,
co-transmitted, whenever it draws more than one lineage the mosquito carries. Both lambdas are
pre-truncation rates, not mean counts: `uptake_moments(lambda, k)` gives the moments a setting
produces and `uptake_calibrate(mean, var)` converts a target mean, and optionally a variance, into
the lambda and `K` to set. An uptake whose 10,000 gamma-Poisson tries all give 0 draws its count
exactly, by inversion of the zero-truncated negative binomial.

**Why.** Truncation, not inflation, is the right family: `P_TRANS_HV` and `P_TRANS_VH` already
decide whether the event happens, so the zeros of a feeding dataset are represented upstream and
inflating here would count one failure twice. Between-donor gametocytaemia variation is likewise
emergent, so `K_OOCYSTS` carries only the residual within-feed overdispersion that vivax membrane
feeds show (a zero-inflated negative binomial describes 96 direct membrane-feeding assays of 22,236
mosquitoes, Miura et al. 2020); the sampling scheme follows Wong et al. (2018). Conditioning away
the zeros pulls the realised mean up — λ = 1 gives 1.58, λ = 2 gives 2.31, and λ = 2 with K = 1
gives 3.0 — so a lambda cannot be read off a reported mean oocyst count, and the two helpers make
that conversion explicit instead of leaving it to the reader. The exact draw is reached by 37% of
uptakes where a positive count is one draw in 10,000 (`LAMBDA_OOCYSTS` about 1e-4, at any K) and by
fewer than one in a million above about 0.003. No study has measured how many genotypes a natural
bite delivers, so the defaults of 2 are an assumption: they let co-transmission happen in every
default run, at a modest rate, and section 14 of the model description lists both lambdas for
sensitivity analysis.

**Alternatives considered.** Settling for a single oocyst when the tries run out; it would make the
count less dispersed than `uptake_moments()` reports. Drawing the count by inversion wherever a
positive count is rare; at the smallest `K_OOCYSTS` the model accepts, 0.001, ordinary settings
already need about 110 tries per positive count, so the exact draw would displace the rejection draw
where rejection works well. Both lambdas at 0 by default, one oocyst per uptake and one genotype per
bite: co-transmission, and the sibling genotypes it brings into one host, would be absent from every
default run.

**How to check.** `tests/testthat/test-model-invariants.R` (the exact draw's mean and variance
against `uptake_moments()`); `tests/testthat/test-edge-values.R` (edge settings, among them the
small `K_OOCYSTS` that `uptake_calibrate()` can return, run to the end);
`inst/validation/calibrate_nu.R` (the share of bites that deliver more than one genotype, at each
batch mean it tries).

### D14 — An outcrossed oocyst gives a full tetrad by per-chromosome four-strand meiosis, with meiosis I pairing at the centromeres

**Decision.** Each oocyst draws two gametes independently and uniformly, with replacement, from the
host's gametocyte pool, and outcrosses exactly when the two differ, so selfing and outcrossing
emerge; the pool's Simpson index is logged as `p_outcross` and drives nothing. An outcrossed
zygote's two genomes replicate into four chromatids per chromosome; chiasmata fall at the union of
two gamete-rate breakpoint draws within each chromosome, each joining a random chromatid of each
parent; meiosis I sends the two chromatids carrying the same parent at the chromosome's annotated
centromere to one pole, drawn independently for each chromosome; meiosis II gives one of them to
each product of that pole. The four products therefore form two pairs of meiosis II sisters on every
chromosome, and every product carries the full map of crossovers. A product with no switch at all is
a whole copy of its starting parent and keeps that parent's lineage. A selfed oocyst's four
identical products are stored as one object with `spz_weight = 4`, which the sampling at a bite
corrects for (D13).

**Why.** This is canonical four-strand meiosis without chromatid interference, and independent
uniform gamete draws give both the marginal P(outcross) = 1 − Σpᵢ² and the correct pairing P(i,j) ∝
pᵢpⱼ (Wong et al. 2018), which rolling outcrossing first and then picking two distinct lineages
would not. The centromeres are the centres of the features `PVP01_CEN01`–`PVP01_CEN14` of the PvP01
annotation (Sanger/GeneDB, released with Auburn et al. 2016), held as fractions of each chromosome
so that `GENOME_MODE = "core"` keeps them (D16); the split barely depends on where they sit
(chromosome midpoints give sisters 0.614 and the other pairs 0.193, the annotated positions 0.603
and 0.198). A sequenced *P. falciparum* tetrad supports the pairing: read against the 3D7 centromere
annotation, its four products form two pairs that carry the same parent at all 13 centromeres that
can be read, sharing 0.67 of the genome against the other pairs' 0.17, where the model gives 0.60
and 0.20; six further pairs of cells in the same study's sibling groups carry the same parent at
every readable centromere (Nkhoma et al. 2020, Figure 6 and Data S3; the centromere readings are
ours, from the published figures). Since IBD and relatedness methods are scored against this
genealogy as ground truth, the shape of same-oocyst sharing is itself a modelling target.

**Alternatives considered.** A random deal of each chromosome's four chromatids to the four
products, as in the simulator of Wong et al. (2018); it drops meiosis I and makes every same-oocyst
pair look alike, around 1/3 with sd about 0.08, which would flatter the methods scored against the
model. Pairing the products without centromeres, by a deal fixed at each chromosome's start; it
gives sisters 0.50 and the other pairs 0.25, so the centromere positions matter. Gene conversion,
which the same tetrad shows in two tracts (872 bp and 181 kb); one tetrad gives no rate to fit. An
obligate chiasma per bivalent; no source supports one for *P. vivax*, and it would move sisters to
about 0.58 and the other pairs to about 0.21, about twice the effect of the centromere positions.

**How to check.** `tests/testthat/test-model-invariants.R` runs the model's own `seedTetrad()` for
400 tetrads with a probe at every centromere: 2:2 segregation at every probe, two sister pairs per
oocyst carrying the same parent at all 14 centromeres, sister and other pairs' sharing against 1/3 +
(2/3)e^(−3rt) and 1/3 − (1/3)e^(−3rt) at map distance rt from the centromere, each product's
crossovers against the map, the reciprocal pairing on chromosomes with one crossover, and the same
centromere fractions in `"core"` mode and with a `CHROMOSOMES` subset.
`inst/validation/audit_arg_recombination.R` and `inst/validation/audit_arg_junctions.R` check the
same machinery in a run's recorded genealogy.

### D15 — The 14 chromosomes are one concatenated linkage group

**Decision.** The genome is one SLiM chromosome holding the PvP01 chromosomes end to end, with the
within-chromosome rate along each and r = 0.5 at each of the 13 junctions. `seedTetrad()` discards
the junction breakpoint draws and assorts the chromosomes itself at the centromeres (D14), which
gives every product exactly r = 0.5 across a junction. The model holds the chromosome starts on the
concatenated map, and `chrom_map()` recovers the layout in R from the genome parameters the run
recorded, so the segment functions can cut at the junctions (D18).

**Why.** One concatenated linkage group is genetically equivalent to explicit chromosomes for every
IBD quantity the package reports, and it keeps one `.trees` file per run. SLiM's explicit
multi-chromosome mode writes a directory of per-chromosome tree sequences, so the whole IBD pipeline
would have to aggregate across 14 files; its one unique payoff, per-chromosome trees for
chromosome-level selection scans, is of no use in a model with no mutation and no selection (D17).

**Alternatives considered.** Explicit SLiM chromosomes — that payoff against an IBD pipeline spread
over 14 files per run.

**How to check.** `tests/testthat/test-model-invariants.R` (a product changes parent at a junction
with probability 1/2); `inst/validation/audit_arg_recombination.R` (breakpoint spikes at all 13
junctions on a uniform within-chromosome background) and `inst/validation/audit_arg_junctions.R` (a
product's parent decorrelates across a junction, φ ≈ 0, while staying linked over 200 kb within a
chromosome, φ ≈ 0.71); `tests/testthat/test-run.R` (the `chrom_map()` lengths sum to the genome
length the run recorded).

### D16 — GENOME_MODE is a coordinate switch, and one uniform crossover rate serves both modes

**Decision.** `GENOME_MODE = "full"`, the default, gives the exact PvP01 chromosome lengths, 24.2
Mb; `"core"` gives 90% of each, a proxy for the subtelomere-masked core genome of the MalariaGEN Pv4
release, about 21.8 Mb. Both modes use one uniform within-chromosome crossover rate of 7.4e-7 per
bp. `RECOMB_RATE` overrides that rate for sensitivity analysis, and `CHROMOSOMES` chooses which of
the 14 chromosomes to simulate.

**Why.** 7.4e-7 per bp (13.5 kb/cM) is the euchromatic core rate of Jiang et al. (2011), a *P.
falciparum* estimate used for want of a comparable published vivax genome-wide map; the genome-wide
figure (9.6 kb/cM) is higher only because it averages in hyper-recombinogenic subtelomeres, which a
uniform map does not place spatially, so spreading it everywhere would overstate recombination
exactly where the analysable genome and the IBD signal lie. `"full"` is the default because the
genome carries no mutation and all genetics come from exact recorded ancestry (D17): subtelomere
masking exists in real data because those regions are hard to genotype and align, which a simulation
with exact ancestry is not subject to, so the complete physical scaffold is the honest
representation and `"core"` is for deliberately mimicking a masked panel. One rate for both modes
leaves them directly comparable, differing only in extent, so the coordinate choice can be isolated:
full gives about 17.9 Morgans of within-chromosome map, core about 16.1.

**Alternatives considered.** Masking the subtelomeres by default — it imports a genotyping
limitation the simulator does not have. A rate per mode, each the self-consistent uniform rate that
reproduces its own extent's empirical map length — it confounds extent with recombination and leaves
the two modes incomparable. The exact Pv4 core-region coordinates in place of the proxy — they are
not shipped with the package, and the proxy affects only the mode that is not the default.

**How to check.** `tests/testthat/test-params.R` (only `"full"` and `"core"` are accepted, whatever
the spelling); `inst/validation/check_parameters.R` (`"core"` shortens the genome);
`tests/testthat/test-model-invariants.R` (`"core"` and a `CHROMOSOMES` subset place each centromere
at the same fraction of its chromosome).

### D17 — The genome carries no mutation, so relatedness is identity by descent read from a forest

**Decision.** The mutation rate is 0 and all genetics come from the recorded genealogy. Every
founder template and every import is an independent root, so the genealogy is a forest and most
pairs never coalesce; tskit's branch statistics, `diversity()` and `divergence()` among them,
therefore measure no relatedness here. The package's two measures are `ibd_pairs()`, where a pair's
IBD counts every stretch over which the two genomes descend from a common ancestor genome, of any
age and any length, and `arg_tmrca()`, which gives the span-weighted mean depth over the part of the
genome where a common ancestor exists, the fraction of the genome where one does, and that
ancestor's raw node time alongside. `min_span` does not change a pair's IBD, and a caller who wants
a recency rule writes it as a filter on the returned segments' `depth`.

**Why.** Independent roots make a forest — about 93% of pairs never coalesce — so a statistic that
scores a non-coalescing span as zero divergence is not a relatedness measure. A depth from the older
of the two genomes makes one relationship give one number wherever in the run it occurred: a clonal
parent-child pair reads 0 whether it arose on day 10 or day 390, while a raw node time is dominated
by when the pair was sampled. Keeping the recency rule in the caller's filter keeps the time
convention visible in the analysis that chooses it rather than hidden behind an argument, and the
data support no single obvious convention: about a third of segments join genomes sampled 100 or
more ticks apart, and in 0.3% the older genome is itself the common ancestor, where depth from the
older is 0 and depth from the younger is not.

**Alternatives considered.** A cutoff on the common ancestor's absolute node time, as tskit applies
one: node times run backwards from the final tick, so the cutoff is a window pinned to the end of
the run, unrelated to the pair being asked about, and it erases exactly the long-dormant relapses
whose ancestry is old by definition — scored that way, relapses rank below local reinfections for
relatedness to earlier cases, inverting the comparison the measure exists to make. A depth argument
on `ibd_pairs()` itself — it presumes one depth convention.

**How to check.** `tests/testthat/test-model-invariants.R` (`arg_tmrca()` gives one number for one
fixed relationship whatever the sampling day, and never a negative depth; the segments reconstruct
each pair's IBD); `tests/testthat/test-ibd.R` and `tests/testthat/test-arg.R`; the ARG and IBD
checks of `inst/validation/audit_events_and_ibd.R` (the genealogy is a forest, relapse-to-origin IBD
is exactly 1, distinct founder lineages are unrelated).

### D18 — An IBD segment is a maximal stretch from one common ancestor, cut at the chromosome junctions

**Decision.** A segment is a maximal stretch over which two genomes descend from one common ancestor
genome, cut at each chromosome junction: abutting rows with the same ancestor genome merge, a change
of ancestor starts a new segment even with no gap, and `min_span` drops the pieces shorter than it
after the cut, a piece of exactly `min_span` being kept. A pair gets the same segments whatever
other genomes are passed, and `ibd_pairs(segments = TRUE)` and `ibd_segments()` give the same rows,
each with the chromosome, both coordinate systems, the length, the common ancestor and its age. A
pair's IBD fraction counts every shared stretch, of any age and any length, so `min_span` does not
move it (D17).

**Why.** A segment whose extent depends on which other genomes were asked for, or that runs across
two chromosomes, cannot be compared, drawn or counted. tskit's rows are pieces of the paths through
the trees, so a pair's ancestry splits into rows differently depending on the genomes the tree was
simplified to, and a stretch continues across a junction whenever both genomes carry the same parent
on both sides of it (D15). Merging on the ancestor genome rather than on its time keeps two
ancestors born on one tick distinct. At `min_span = 0` the merged and cut segments add up to each
pair's IBD, which ties the two views together.

**Alternatives considered.** Taking tskit's rows as they come — path-dependent, and a row can run
across a junction. Merging on the ancestor's time — two ancestors born on one tick would merge.
Carrying the length floor into the pair fraction — a floor applied before the cut depends on the
other genomes passed, and applying it after the cut costs four to five times the time and memory.

**How to check.** `tests/testthat/test-ibd.R`: merging, cutting and flooring on a synthetic table;
on a run, the same segments from any set of genomes and from both functions; no segment across a
junction or below `min_span`; the segments adding up to each pair's IBD; one schema for every pair;
a genome with itself; pairs counted from 1.

### D19 — IBD comes back as a pedigree-keyed object, and an id the tree sequence does not hold reads NA

**Decision.** `ibd_pairs()` returns a `vivax_ibd` object keyed by what callers hold, pedigree ids:
`pairs` (every pair of the genomes found, `a < b`, zeros included), `samples` (each genome's node
time and the day it arose), `missing`, optional `segments` (with the common ancestor's pedigree id
and its depth), the sequence length, the floor used and the run's `end_tick`. `ibd_lookup()` reads
one pair and gives 1 for a genome with itself; `as.matrix()` and `as.data.frame()` give the table. A
requested id the tree sequence does not hold warns, is listed in `missing` and reads `NA`; asking
only for such ids is an error that gives the reason the run's record holds, whether it never started
or ran with `RETAIN_FULL_HISTORY = FALSE`, and the readers built on it then give their `NA` row or
empty table. The Python side is one private module, called with arguments and returning arrays.

**Why.** An object keyed by pedigree ids needs no translation from the ids callers already have, and
`NA` for a genome the run does not record cannot be mistaken for an unrelated pair, which a 0 can.
Returning each genome's node time is what makes a depth filter writable by the caller (D17). A
module called with arguments pastes nothing into code and passes no table through a file, which
costs one array conversion rather than a round trip through text.

**Alternatives considered.** A sparse table with implicit zeros — every reader would have to fill
them in. 0 for an absent id — the confusion `NA` avoids. Returning the simplified tree itself — a
live Python object cannot be saved, so an R object's validity would depend on a live session,
breaking `saveRDS()` and grid workers.

**How to check.** `tests/testthat/test-ibd.R` (the object read back; its shape on a run; lookups of
relapses against their origin; segment ancestors and ages; ids not in the tree sequence; a run with
no parasite) and the IBD checks of `inst/validation/audit_events_and_ibd.R`.

### D20 — Marker genotypes are painted on the root genomes from a panel's allele frequencies

**Decision.** `simulate_genotypes()` takes biallelic loci as chromosome, 1-based position within
that chromosome and PLAF, the ALT frequency. Each root genome — a founder template or an import
(D39) — draws its allele at each locus independently, every genome carries the allele of the root
genome it descends from there, and `roots` records that root genome for every genome and locus, the
truth behind every call. Samples are clones, one genome each, or host samples grouped by `sample_id`
(0, 1, or 0.5 where both alleles are present); `error` can replace a call and `missing` follows the
errors. `write_vcf()` writes VCF 4.2 with one contig per simulated chromosome, haploid clone calls,
unphased diploid host calls and the painting frequency as `PLAF` in `INFO`. The root genomes'
alleles are drawn for every root genome of the run before any sample is assembled.

**Why.** In a neutral model with no mutation (D17), painting the roots is exact: every genome
descends at each position from exactly one root genome and root genomes are unrelated, so two
genomes share a painted allele exactly where they share a common ancestor (D18). PLAF is what users
have from a panel. Drawing for every root genome of the run first makes a genome's alleles
independent of the company it is asked for in, of the grouping into samples, and of the error and
missingness laid over them.

**Alternatives considered.** Coalescent mutations over the recorded genealogy — a mutation never
falls on a root, so no site separates two root genomes, and the frequencies would come from a
coalescent prior rather than from the panel. Recapitation before mutation — the same prior, and it
needs a demography the run does not provide. Painting in tskit and reading its genotype matrix — a
genome alone under its root reads as missing wherever that root carries REF, and tskit's VCF writer
gives one contig in concatenated coordinates rather than one per chromosome. Drawing alleles only
for the root genomes behind the genomes asked for — a genome's alleles would depend on which others
were asked for.

**How to check.** `tests/testthat/test-genotypes.R`: every genome carries its root genome's allele,
by the bookkeeping, by painting the root nodes in tskit and by walking the unsimplified genealogy;
root genomes are their own roots; the same root genome exactly where `ibd_segments()` finds a common
ancestor; the root genomes and their days; the seed's invariance; host samples; the VCF; refusals;
the schema. `inst/validation/audit_events_and_ibd.R` repeats four of those checks on each audit run.

## Hypnozoites, relapse and the reservoir

The reservoir the model exists for: how a bite seeds hypnozoites, how they wake or die, and how
relapse is counted once a relapse can hide inside an infection already there.

### D21 — A bite seeds one geometric hypnozoite batch, split over the genotypes it delivers, at a fitted mean

**Decision.** Each genotype a bite delivers gives a blood-stage primary with probability `P_PRIMARY`
(1 by default). Independently, the bite draws one hypnozoite batch, geometric on 0, 1, 2, ... with
mean ν = `NU_HYPNOZOITE_BATCH` (times `NU_HYPNOZOITE_BATCH_BOOST` from `RESURGENCE_TRIGGER_TICK`),
and each of its hypnozoites goes to one of the delivered genotypes uniformly at random; a genotype
whose share is above 0 lays down a batch that is a clone of it, as its primary is, so the two are
IBD 1. A genotype can give both, either or neither, and one that seeds only hypnozoites makes no
case (D29). An import draws its primary and its batch the same way (D38); founder hosts start with a
blood-stage clone and no batch. The mean ships at 3, the value fitted at moderate transmission to
four calibration targets, which section 13 of the model description gives with their sources; the
fit counts every relapse, masked ones included (D26). `?vivax_params` says the 40-host defaults are
an illustrative preview, not a calibrated setting, and that the mean should be re-fitted when
transmission or scale change.

**Why.** The batch is the one White et al. (2014) define, with mean ν per bite, and their base model
treats a primary as certain given a successful inoculation, hence `P_PRIMARY = 1`. Drawing the batch
once per bite and splitting it keeps ν that quantity and keeps co-transmission from inflating the
reservoir. A default should be the value the documentation and its figures use, and the one the fit
supports.

**Alternatives considered.** A placeholder default apart from the fitted value was considered: no
documented figure uses one, and the fit does not support it.

**How to check.** `inst/validation/calibrate_nu.R` repeats the fit and checks that the shipped value
sits inside all four targets; `tests/testthat/test-params.R` keeps the default of `vivax_params()`
equal to the model's own; `tests/testthat/test-model-invariants.R` checks that an inoculation's
primary and its batch are the same genome, for bites and imports alike;
`inst/validation/check_parameters.R` checks that a batch mean of 0 gives no relapse and that the
batch boost enlarges the reservoir.

### D22 — Activation and death are competing daily risks per hypnozoite, clamped so a boost stays valid

**Decision.** Each day a batch not in latency (D25) draws its activations, Binomial(`batch_count`,
α), and then its deaths among the rest, Binomial(remaining, μ/(1 − α)), with α = `ALPHA_ACTIVATE`
and μ = `MU_HYPNOZOITE`. Each hypnozoite thus activates with probability α and dies with probability
μ that day, so its fate is settled after a geometric time with mean 1/(α + μ), and it activates
rather than dies with probability α/(α + μ). From `RESURGENCE_TRIGGER_TICK` α is multiplied by
`ALPHA_ACTIVATE_BOOST`; on every day the effective α is clamped at 1 − μ, and the death probability
is held at 1 at most, since at α + μ = 1 the ratio can round just above it. An extreme reactivation
boost therefore saturates rather than failing.

**Why.** The two binomials are each hypnozoite's daily choice between activating, dying and staying,
drawn for the whole batch at once, with activation and death competing as in White et al. (2014) and
read as daily probabilities, as every rate in the model is (D3). The clamp keeps both probabilities
valid, α below 1 and μ/(1 − α) at most 1: no more than a whole batch can activate in a day. The two
rates are literature values, not fitted: α = 1/170, an average carriage of about 170 days (Mehra et
al. 2024b; White et al. 2014 estimate 1/332), and μ = 1/442 (Mehra et al. 2024a, who credit White et
al. 2014, whose table gives 1/425).

**How to check.** `tests/testthat/test-run.R` checks that the share of resolved hypnozoites that
activated is α/(α + μ); `tests/testthat/test-edge-values.R` runs α + μ = 1, the edge of the clamp,
to the end, every batch resolving on its first dormant day; `tests/testthat/test-params.R` accepts
an activation boost of 1000, which the model clamps; `inst/validation/check_parameters.R` checks
that an activation boost raises the relapse fraction.

### D23 — One relapse per activation event, clonal to the batch it woke from

**Decision.** However many of a batch's hypnozoites activate on one day, they seed one blood-stage
genome, cloned from the batch, so relapse-to-origin IBD is exactly 1. The batch loses all of them;
the one `hypnozoite_activation` row names the batch in `parent1_pedigree_id` and carries the count
in `n_hypnozoites`, and one detection draw (`P_DETECT_RELAPSE`) and at most one course follow.
Batches that wake on the same day give a relapse each, and in one host those relapses are clones of
one episode (D46). Relapses therefore grow less than in proportion to the batch mean: a larger batch
has more hypnozoites waking together.

**Why.** The hypnozoites of a batch come from one inoculation and are genetically identical, so
those that wake together are one recurrent infection.

**Alternatives considered.** A genome for each activating hypnozoite was considered: it would
over-count relapses, keep the lineage in the blood until every copy cleared, draw detection and
radical cure once per hypnozoite, and weight that lineage by its count in the gametocyte pool a
mosquito draws from.

**How to check.** `tests/testthat/test-run.R` and `inst/validation/audit_events_and_ibd.R` check
that every relapse has IBD 1 with its batch; the audit also checks that every activation names an
established batch and that a relapse shares every root genome with it;
`tests/testthat/test-model-invariants.R` checks that a batch's activations and deaths, counted in
hypnozoites, never exceed what it held.

### D24 — A hypnozoite batch empties only by activation and death

**Decision.** No cap clears a batch by its age: outside radical cure, by a case course (D31) or an
MDA round (D35), a batch ends when activation and death have taken its last hypnozoite.
`MU_HYPNOZOITE` must be above 0, in `vivax_params()` and in the model itself, so every batch ends. A
batch that runs out in `reproduction()` has its own count set to 0 and is killed by the cleanup at
the top of the same day's `early()`, since `killIndividuals()` cannot run in `reproduction()`; one
emptied by radical cure, which runs later in `early()`, goes at the next day's cleanup and is
skipped until then. The mark sits on the batch itself rather than in a global list, which avoids the
Eidos scoping trap of section 18 of the model description.

**Why.** White et al. (2014) give hypnozoites an exponential residence with no cap. A finite
lifespan belongs in the daily death probability, not in a cutoff, which would clear whole batches at
exactly the cap age and cut the long late-relapse tail the model is built to carry.

**Alternatives considered.** An age cap kept as a silent safety net was considered: wherever it
binds it bends exactly that tail, and `MU_HYPNOZOITE > 0` already ends every batch.

**How to check.** `tests/testthat/test-model-invariants.R`, `inst/validation/audit_events_and_ibd.R`
and `inst/validation/check_parameters.R` check that no batch draws more hypnozoites than it held,
which a batch left in place after running out would do; the audit also checks that the hypnozoites
left unresolved equal the final reservoir in `state.tsv`; `tests/testthat/test-params.R` and
`tests/testthat/test-slim.R` refuse `MU_HYPNOZOITE = 0`.

### D25 — Latent hypnozoites die, and a batch moves through latency as one unit

**Decision.** With `K_LATENCY` = k above 0, every new batch starts k latency stages away from
activation. Each latent day it first loses Binomial(`batch_count`, μ) hypnozoites, logged as
`hypnozoite_clearance` like any other death, and then, if any are left, moves on one stage with
probability δ = `DELTA_LATENCY`; it cannot activate while latent, nor on the day it leaves latency.
A hypnozoite therefore relapses with probability [δ(1 − μ)/(1 − (1 − δ)(1 − μ))]^k × α/(α + μ),
0.480 at k = 2 and the shipped rates against 0.722 without latency, and because survival favours the
batches that move through quickly, one that leaves latency alive has been latent k/(1 − (1 − δ)(1 −
μ)) days on average (163 at k = 2), not k/δ. `DELTA_LATENCY` must lie in (0, 1] when `K_LATENCY` is
above 0.

**Why.** It is the model the documentation cites: White et al. (2014) apply death in every latent
compartment and track latency per infection, which is why the batch moves as one unit. At
`DELTA_LATENCY = 0` a temperate batch could never leave latency, and so could never relapse.

**How to check.** `tests/testthat/test-run.R` imports 300 infections at `K_LATENCY = 2`, with no
onward transmission, runs until every batch has resolved, and checks the share of hypnozoites that
activated against the formula above; it also checks that latent deaths in imported batches carry
`is_imported = 1`. `tests/testthat/test-params.R` refuses `DELTA_LATENCY = 0` at `K_LATENCY = 2`,
and `inst/validation/check_parameters.R` checks that temperate latency delays each batch's first
relapse.

### D26 — A relapse into a patent clone of its own lineage is counted like any other, and flagged

**Decision.** A relapse can begin while a patent clone of its own lineage is in the host, most often
the primary of the same bite or an earlier relapse of the same batch; at the shipped defaults about
80% of relapses do. Such a relapse is a new genome, counted as a relapse and drawn for detection
like any other, so it can be reported and treated. Its `hypnozoite_activation` row has
`into_patent_lineage = 1` (0 on other activations, -1 on every other row), and
`exclude_masked = TRUE` leaves these masked relapses out of `classify_cases()`,
`relapse_fraction()`, `force_of_infection()`, `recurrences()` and `summarise_run()`, and passes over
them in `time_to_relapse()` (D27). The activation runs in `reproduction()` and sees only clones born
on earlier days (one that ends that day counts), so a relapse that starts on the same day as another
clone of its lineage, from another batch waking, a recrudescence or a bite, is not flagged for that
clone: two batches of one lineage that wake together, with no earlier clone of it in the host, are
both 0.

**Why.** A study that genotypes the host would see only the infection already there, not a
recurrence, while hypnozoites did wake and the model's observation layer can report the relapse;
which count an analysis needs depends on its question, and the flag gives it both.

**Alternatives considered.** Merging such an activation into the clone already there was considered:
it would build one reading into the model.

**How to check.** `tests/testthat/test-run.R` rebuilds each host's patent clones from `events.tsv`
in a run with treatment failures, recrudescence and an MDA round, and the flag matches on every
activation; `tests/testthat/test-episodes.R` checks that `exclude_masked` removes exactly the
flagged activations, and `tests/testthat/test-summary.R` that it leaves them out of every relapse
column of `summarise_run()` and changes no other column.

### D27 — Time to first relapse is measured per hypnozoite batch

**Decision.** `time_to_relapse()` gives one row per hypnozoite batch, timed from the day it is laid
down (for an import that also brings a blood-stage infection, the day after it) to its first
relapse, which the activation's `parent1_pedigree_id` ties to the batch. A batch with no relapse is
censored on the day its last hypnozoite leaves it, by death or radical cure (`emptied = "censor"`,
the default), or at `SIM_DAYS`; `emptied = "follow"` follows an emptied batch to `SIM_DAYS`, the
curve a trial reports, where radical cure shows as relapses prevented. With `exclude_masked = TRUE`
the clock stops at the first relapse a study could see as a recurrence, although a masked relapse's
hypnozoites still leave the batch. `plot_relapse_km()` draws the Kaplan–Meier fraction of batches
not yet relapsed.

**Why.** A relapse comes from a batch, and the time from the infection that laid a batch down to its
first relapse is the quantity cohort and trial data report; censoring at emptying gives the time to
relapse of the hypnozoites that live.

**Alternatives considered.** A clock per host was considered: a host's first relapse can come from
the batch of a later bite, which mixes the wait for reinfection into the measure, and since founders
carry no batch, a clock started at founder seeding makes the curve's shape depend on how many hosts
are founders.

**How to check.** `tests/testthat/test-relapse-metrics.R`: a fixture's five batches, censored,
followed and with a masked relapse passed over; the curve's steps; a run with no batch; on an import
cohort, every batch timed from its establishment, and the batches radical cure reached first ending
as emptied; a batch emptied in a run that goes extinct; every liver-stage row that names a genome
takes hypnozoites out of its batch. `inst/validation/check_parameters.R` checks that temperate
latency delays each batch's first relapse.

### D28 — The relapse rate is per person-year at risk, and the relapse fraction has two readings

**Decision.** `summarise_run()` (D50) gives the relapse rate as `relapses_per_py_at_risk`: every
relapse in its window over `person_years_at_risk`, the person-time of the hosts ever infected. A
host's time starts the day after its first establishment, a primary or a hypnozoite batch from a
bite or an import, and runs to the window's end, also after an extinction; founder seeding does not
put a host at risk. The relapse fraction has two readings, both counted per clone onset and both
including each host's first case: `relapse_frac_acquisitions`, relapses over imported, local
reinfection and relapse cases, the unit of the relapse-fraction target in section 13 of the model
description, and `relapse_frac_episodes`, relapses over every case, recrudescence included.
`relapse_fraction()` and `summarise_run()` give both, each `NA` when its denominator is 0, and the
two are equal on a run with no recrudescence, as at the defaults. Masked relapses (D26) are counted
unless `exclude_masked = TRUE`; `?relapse_fraction` says the batch mean is fitted with every relapse
counted, and what leaving the masked ones out gives.

**Why.** A host never infected cannot relapse, while one whose first infection left only a batch
can, and the run's length, not its last event, is the period observed. A recrudescence is a new
episode but not a newly acquired parasite (D33), so once recrudescence is on both readings are
needed, and names that say their unit keep them from being swapped. A fraction with nothing to
divide is not 0.

**Alternatives considered.** Dividing by the person-time of every host was considered: it counts
hosts never at risk, so the rate falls with the share never bitten. Counting every ever-infected
host over the window's whole length was considered too: it counts a host's days before its first
infection.

**How to check.** `tests/testthat/test-relapse-metrics.R`: both fractions, with and without the
masked relapse, and a run with no case. `tests/testthat/test-summary.R`: a fixture's relapses and
person-time, with and without a window, and a run with no infection; on a run, every activation
counted and person-time to `SIM_DAYS` after an extinction; `summarise_run()`'s columns against
`relapse_fraction()` and the person-time over several windows, and `exclude_masked` in every relapse
column. `tests/testthat/test-treatment.R`: on a run with recrudescence, each reading equals the
share of relapses among its cases, and the acquisitions reading is the larger.

## Detection, treatment, recrudescence and MDA

The observation layer and the drugs: what a surveillance system sees, what a course does to a host
and its reservoir, what survives one, and what a mass round reaches that a case course cannot.

### D29 — Detection is drawn once per bite, import, relapse or recrudescence, and an infection with no blood stage is no case

**Decision.** Whether a case is reported is drawn once per bite, for every genotype the bite
delivers, at `P_DETECT_PRIMARY`; once per import, at `P_DETECT_PRIMARY` too; once per activation,
which makes one relapse (D23), at `P_DETECT_RELAPSE`; and once per recrudescence, at
`P_DETECT_RECRUD`. The outcome is the `detected` flag on the case's onset rows. A detected case can
be prescribed a course (D30, D31), and the reported cases, `detected == 1` in `classify_cases()` and
`force_of_infection(detected_only = TRUE)`, are the observed line list; the true one holds every
case. A bite or an import that establishes only a hypnozoite batch, which `P_PRIMARY` below 1
allows, makes no blood stage, so it is no case: no row of it carries a detection flag, and it gets
no course. Detection depends on origin only through the three probabilities, so active and passive
surveillance cannot be told apart.

**Why.** The genotypes one bite delivers are separate cases of one episode (D46), but one person
presents with them, so they are reported and treated together, never for one genotype and not the
other. Detection differs by origin because relapses are often subpatent, so a realistic
`P_DETECT_RELAPSE` is much lower than a primary's (0.15, say), and Ross et al. (2025) put recurrence
detectability at 0.13–0.49 of a primary's, so `P_DETECT_RECRUD` should usually sit well below
`P_DETECT_PRIMARY`. A batch in the liver cannot have caused symptoms, and a course given for it
would clear the host's other clones and, with radical cure, the batch the same bite had just seeded.

**How to check.** `tests/testthat/test-treatment.R`: at `P_PRIMARY = 0`, with full coverage and
undetected relapses, no course is given. `tests/testthat/test-load.R` and
`inst/validation/audit_events_and_ibd.R`: the detection flag is set only on the blood-stage onset
rows of bites, imports, relapses and recrudescences, never on founder seeding or a batch.
`inst/validation/check_parameters.R`: lower detection probabilities thin the observed line list, and
at the defaults every case is reported. `tests/testthat/test-windows.R`:
`force_of_infection(detected_only = TRUE)` counts the reported cases.

### D30 — Case management has its own window, separate from vector control

**Decision.** `TREATMENT_START` and `TREATMENT_END` bound case management: a case detected on days
`TREATMENT_START` to `TREATMENT_END - 1` can be prescribed a course. They default to the whole run,
so a coverage alone is the switch: set `BLOODCURE_COVERAGE` and detected cases are treated. One
helper, `caseManagementOn()`, gates every case-based course; it is asked in `decideTreatment()`,
which every case goes through. `INTERVENTION_*` is vector control and scales biting only (D37): it
does not touch treatment, and the treatment window does not touch biting. The window gates the
decision, not the delivery, so a case detected on the programme's last day completes its course
after the window closes. A coverage set against a window that holds no day of the run, one that is
closed, ends by day 1 or opens after the run ends, is refused, by `validate_params()` and by the
model itself (D43).

**Why.** Case management is drugs given to diagnosed patients and vector control is bed nets and
spraying; nothing should make one imply the other, and a scenario can treat cases without changing
biting. Expressing the gate once keeps the sites where cases arise from drifting apart. A patient
who presents on the programme's last day is not refused treatment. A coverage against a closed
window would silently give no treatment at all.

**Alternatives considered.** One window for both programmes, with a warning when a coverage is set
without it: it keeps a coupling that has no modelling reason, and forces a change in biting on
anyone who wants treatment. A treatment efficacy of its own within the shared window: the timing of
the two programmes would stay tied.

**How to check.** `tests/testthat/test-treatment.R`: a coverage alone gives courses; opening the
vector-control window alone gives none; a bounded window gives courses only from its first day, with
a short tail after it closes; and a coverage against a closed window is refused. The window test
asserts when courses come rather than how many, because a window open throughout a run can give
fewer courses than a late one, by preventing the cases it would treat. `tests/testthat/test-slim.R`:
a run started directly with SLiM refuses a coverage against a closed window too.

### D31 — Treatment is decided where the case arises and applied to the host later; radical cure is nested inside it

**Decision.** `decideTreatment()` runs where a detected case arises, at a bite, an import, a relapse
or a recrudescence. Inside the treatment window (D30), and when the host has no course pending, it
draws a blood-stage course at `BLOODCURE_COVERAGE` and, nested inside it, radical cure at
`RADCURE_COVERAGE`, and schedules the course (D32). `applyTreatment()` gives the course in `early()`
on the day it lands, before the day's natural clearance, to the host rather than to one lineage:
every patent blood-stage clone the host carries, except one patent only since that day (D32), clears
unless it survives, at `P_BLOODCURE_SURVIVAL` (D33), and radical cure, if drawn, acts on every
hypnozoite batch the host carries. Radical cure is liver-stage only and never clears the current
blood stage, and `RADCURE_COVERAGE` is the probability of radical cure given a course, not a
coverage of its own. Its efficacy is per person by default (`RADCURE_MODE = "person"`): one
all-or-nothing draw, in which the whole reservoir survives with probability `P_RAD_SURVIVAL` (0.29,
about 71% cure). `RADCURE_MODE = "hypnozoite"` draws each hypnozoite's survival independently, so a
course can clear part of a reservoir; there `P_RAD_SURVIVAL` is a per-hypnozoite probability, and
the same number is not the same efficacy. A host takes one course at a time: a second case while a
course is pending schedules none, and a treated case stays infectious until the drugs land.

**Why.** `early()` is the first point at which the genomes made in `reproduction()` are visible, so
a course applied there reaches the batch seeded by the infection behind the case, an import's
included, whose batch is cloned the day after it arrives (D38). Drugs clear the blood-stage
parasites in a person, not one clone, so the course treats the patient, and a detected relapse gets
a schizonticide like any other case. Radical cure is the 8-aminoquinoline class (primaquine,
tafenoquine), the one lever on the dormant liver reservoir and a class apart from the blood
schizonticides; primaquine is not given without a schizonticide, hence the nesting. Per-person
efficacy is how clinical efficacy is reported, as recurrence prevented per patient, and 0.29 sits
near the 70% effectiveness White et al. (2018) assume for a 14-day primaquine course, given poor
adherence. The per-hypnozoite mode is the only one in which some of a host's genotypes survive
radical cure, which is what reshapes relapse diversity after treatment.

**Alternatives considered.** A treatment flag on the new lineage, set where the case arises: it
treats one clone rather than the patient, and cannot reach the batch the same bite seeds, which SLiM
keeps out of sight until `reproduction()` ends. Per-person efficacy alone: it cannot clear part of a
reservoir, so the per-hypnozoite mode is offered beside it.

**How to check.** `tests/testthat/test-treatment.R`: a blood-stage coverage alone gives courses, and
radical cure nested in them clears hypnozoites. `tests/testthat/test-relapse-metrics.R`: on a cohort
of imports, a case's radical cure empties the batch its own import seeded before any of it
activates, and with `RADCURE_COVERAGE = 0` none is emptied; in the per-hypnozoite mode, case and MDA
radical cure take hypnozoites out of the batches they name. `inst/validation/check_parameters.R`:
radical cure shrinks the reservoir, a weaker one leaves a bigger reservoir, and blood-stage
treatment lowers prevalence.

### D32 — A course lands at least a day after its case, and the treatment and recrudescence delays follow the 1/D rule

**Decision.** A course is scheduled `1 + rgeom(1, 1/TREATMENT_DELAY)` days after its case: each day
after the case it arrives with probability 1/`TREATMENT_DELAY`, so the delay is at least one day and
its mean is exactly `TREATMENT_DELAY`, 7 by default. On the day it lands, a case course passes over
any clone whose blood stage began that day. `TREATMENT_DELAY` is at least 1, where 1 means the next
day, the fastest response the model can represent. The delay from suppression to recrudescence,
`MEAN_RECRUD_DELAY`, follows the same rule under its default shape (D34). Both parameters are
finite, from 1 to 10^6 days; `validate_params()` and the model itself (D43) refuse other values.

**Why.** A clone that became patent today cannot yet have caused symptoms, prompted care-seeking or
been diagnosed, and natural clearance follows the same rule: a clone cannot clear on the day its
blood stage began. The rule is needed at both ends. At the decision, because a course landing on the
day of its case would pass over the clone behind the case and be spent all the same, leaving the
patient untreated; at the application, because a course given for one case must not reach a
different clone that appeared on its landing day. Under the 1/D rule every duration parameter means
what it says (D3). MDA waits for no symptoms, so a round is exempt from the application rule (D35).

**Alternatives considered.** A geometric draw of mean `TREATMENT_DELAY` that can return 0, floored
at one day: the zeros pile onto the first day, so the mean exceeds the parameter, 7.125 days for 7.
Accepting `TREATMENT_DELAY = 0` as another name for the next day: two values for one behaviour,
where 1 already means the next day.

**How to check.** `tests/testthat/test-treatment.R`: with a pulse of detected imports as the only
cases, at a delay of 1 every course and every recrudescence comes the next day, and at 4 both delays
have mean 4 with a quarter of them one day long, a share that any mass at zero would inflate.
`tests/testthat/test-params.R` and `tests/testthat/test-param-dictionary.R`: values below 1,
infinite or above 10^6 are refused. `tests/testthat/test-slim.R`: a run started directly with SLiM
refuses values below 1 or infinite.

### D33 — A recrudescence mints a new genome clonal to its holder, and is a case but not an acquisition

**Decision.** A clone that survives a course carries on at full density (`bloodcure_persist`) or,
with probability `P_RECRUD_GIVEN_FAIL`, is suppressed below patency (`bloodcure_subpatent`, which
ends its patent span). A suppressed clone, a holder, is not patent: it is out of `moi_active`
(`state.tsv` counts it in `moi_subpatent`), cannot be detected, and no later case course reaches it.
It is present, counting toward `MAX_MOI`, and transmissible: the bite loop asks whether a host is
infectious, not whether it is patent, and a holder joins the gametocyte pool a feeding mosquito
draws from with probability `SUBPATENT_INFECTIVITY`, 0.15 by default. A holder does not clear
naturally. After its delay (D34) it re-expands as a `recrudescence`: a new genome cloned from the
holder in `reproduction()`, with a pedigree id of its own and its own place in the tree sequence,
dated at the recurrence, while the holder itself is removed in that day's `early()`, before the
courses and clearance. The new genome carries the holder's lineage, and its row gives -1 for
`inoculation_id` and `oocyst_id` and names the holder in `parent1_pedigree_id`. A recrudescence is a
case, detected at `P_DETECT_RECRUD` and treated like any other (D29, D31), and the fourth origin of
`classify_cases()`, so it counts in case counts, incidence and host timelines. It is not an
acquisition (`is_acquisition = FALSE`), so FOI leaves it out, and the metrics that have both
readings, molFOB, the relapse fraction and Rc, report both, by unit (D28, D50, D53).
`P_RECRUD_GIVEN_FAIL` defaults to 0, so a default run has no recrudescence.

**Why.** A recurrence has three causes, relapse, reinfection and recrudescence (Taylor et al. 2019),
and the model represents all three. Minting gives every case a genome of its own: the recurrence's
genome is dated at the recurrence rather than at the original inoculation, and its IBD of 1 to the
earlier episode is a measurement between two genomes, not a genome paired with itself. Because the
new genome inherits the lineage, removing the holder and adding it leaves MOI, prevalence, COI and
the reservoir exactly as they are. Patency, infectiousness and presence are counted apart because a
holder separates them: counted as patent it would make an episode no one can see, and a patent count
in the bite loop would silently block the transmission `SUBPATENT_INFECTIVITY` permits. The
mosquito-feeding studies Koepfli et al. (2021) compile put submicroscopic carriers at 0.10–0.20 of
the infectivity of microscopy-positive ones, so neither 0 nor 1 is defensible. With natural
clearance off while subpatent, every holder recrudesces unless an MDA round clears it (D35), so
without MDA `P_RECRUD_GIVEN_FAIL` is exactly P(recrudescence | the course failed), comparable with
the failure rate of a therapeutic-efficacy study. No new parasite enters the host, so a
recrudescence is no acquisition; a field molFOB, which compares genotypes between visits, often
counts it as a new clone, because the clone is undetectable at the earlier visit, and the gap
between the two readings estimates that bias. Its ids are -1 because a recrudescence is no bite, and
inherited ones would put it in another case's co-transmission group. The default of 0 matches
recrudescence at 0.2–0.3% of recurrences after chloroquine or artesunate (Taylor et al. 2019, n =
1441); the state matters in a chloroquine-resistant setting.

**Alternatives considered.** Reusing the holder's individual: two onsets would share one pedigree
id, against the rule that creation events mint unique ids, and the recurrence's onward transmission
could not be told from the earlier infection's. Re-dating the holder's node at the recurrence: a
node keeps the time its genome is made. Minting without removing the holder: one failed course would
recrudesce again and again.

**How to check.** `tests/testthat/test-treatment.R`: at `P_RECRUD_GIVEN_FAIL = 0` nothing is
suppressed; at 1 every failure suppresses, and each recrudescence has an id of its own, a holder as
its parent, -1 ids and a later day than its suppression, is a case of its own origin and no
acquisition, and nothing clears twice; a recrudescence is IBD 1 to its holder at a TMRCA of 0, in
the holder's lineage. `inst/validation/audit_events_and_ibd.R`, on its configuration with
recrudescence: creation events mint unique pedigree ids, every clearance matches an onset, and none
repeats. `tests/testthat/test-run.R`: a holder counts toward `MAX_MOI` even at
`SUBPATENT_INFECTIVITY = 0`.

### D34 — The recrudescence delay is geometric with a mean of 10.5 days, selectable, and a literature assumption

**Decision.** Each day after its suppression a holder re-expands with probability
1/`MEAN_RECRUD_DELAY`, so the delay is at least one day with that mean, 10.5 days by default: a
median of 7 days, and about 94% within day 28. `RECRUD_DELAY_DIST` sets the shape: `"geometric"`,
the default; `"poisson"`, 1 + Poisson(D − 1) with mean D, a tight-window sensitivity arm; or
`"fixed"`, round(D) days with halves rounded up, 11 at 10.5. Every shape is at least one day. 10.5
days is a literature assumption, pending parameterisation, which section 14 of the model description
lists for sensitivity analysis.

**Why.** The only published fit for this quantity models the time to recrudescence as
exponential (Taylor et al. 2019), and the ceiling of an exponential delay on a day scale is a
geometric one. Every other duration in the model is geometric too (D3). 10.5 days is the
reciprocal of their posterior median rate, 0.095 a day, but they hold the rate with a strongly
informative Gamma(100, 1000) prior, mean 0.1 a day, for identifiability, and the posterior
barely moves from it, so the value is close to the prior rather than an estimate from data. At
this mean nearly all recrudescences fall within the standard 28-day follow-up of a
therapeutic-efficacy study, where a study can see them.

**Alternatives considered.** A Poisson delay as the default: no source supports it, and its spread,
a standard deviation of about 3 days at this mean, is too narrow for regrowth from an unobserved
density; it is offered as the `"poisson"` arm, a tight-window contrast. A gamma delay with a shape
parameter: regrowth is log-linear, so the delay is a logarithm of the surviving biomass, which
argues for a peaked delay with a floor; but it needs a second parameter, and the family with a
published fit behind it is the default.

**How to check.** `tests/testthat/test-treatment.R`: at a mean of 1 every shape gives a one-day
delay, and at 4 the geometric delay has mean 4 with a quarter of delays one day long.
`tests/testthat/test-params.R` and `tests/testthat/test-param-dictionary.R`: a mean below 1,
infinite or above 10^6 days, and a shape outside the three, are refused.
`tests/testthat/test-slim.R`: a run started directly with SLiM refuses a mean below 1 or infinite,
and runs each shape.

### D35 — MDA is one round of the same course, drawn per host, and reaches the clones a case course cannot

**Decision.** `MDA_TICK` gives one round of mass drug administration through the same
`applyTreatment()` as case management: blood-stage drugs and radical cure for every host it reaches,
those carrying only hypnozoites included, at the efficacies `P_BLOODCURE_SURVIVAL`, `P_RAD_SURVIVAL`
and `RADCURE_MODE` set. `RADCURE_COVERAGE`, a probability given a case course, plays no part in it.
The round reaches each host with probability `MDA_COVERAGE`, drawn only when the coverage is below
1, so a round that reaches everyone draws no number. Only hosts not empty at the round are drawn;
the draw does not depend on infection, so the coverage is that of the whole population, and
`mda_bloodcure_given` and `mda_radcure_given` count the courses to those hosts. The round comes in
`early()`, after the day's natural clearance and before its imports, and the treatment window does
not gate it. Because it waits for no symptoms, it reaches two kinds of clone a case course does not:
a clone that became patent earlier the same day, and a holder (D33), which it clears with
probability 1 − `P_BLOODCURE_SURVIVAL`, cancelling its recrudescence. A cleared holder is logged as
`mda_subpatent_clear`, which ends no episode, since its patent span ends at its suppression; a
holder that survives keeps its schedule and is not suppressed again. A patent clone that survives
the round persists or is suppressed as under a case course. `MDA_TICK = 1` is refused when a founder
host carries a copy of another founder's genome, which is seeded on day 2 (D40): a day-1 round could
clear the template before the copy is made.

**Why.** One course for both paths means the two cannot diverge. Coverage is the whole policy
question for MDA, what separates a feasible campaign from an idealised one; without it an imperfect
campaign could only be mimicked through efficacy, which confounds the two. Neither filter a case
course keeps has its reason under MDA, which is given to everyone, symptoms or not: the one-day rule
stands for the delay from patency to symptoms and care (D32), and the holder rule for a course
prescribed for a case other than the holder's. Drugs taken on the day a clone reaches the blood act
on it too.

**Alternatives considered.** A round blind to holders: it would leave the very reservoir it exists
to clear. Drawing a surviving holder's suppression again: it would reset the holder's schedule and
log a second end to a patent span that has already ended.

**How to check.** `tests/testthat/test-treatment.R`: at a round, each holder is cleared
(`mda_subpatent_clear`, never to recrudesce) or left on its schedule with no row that day, a clone
the round has just suppressed is not drawn again, and nothing clears twice; with
`P_BLOODCURE_SURVIVAL = 0`, a day-2 round clears every founder copy seeded that morning and leaves
no blood-stage clone, and the host trajectory plot draws such a one-day clone.
`tests/testthat/test-params.R` and `tests/testthat/test-slim.R`: `MDA_TICK = 1` is refused only when
a founder host carries a copy. `tests/testthat/test-relapse-metrics.R`: in the per-hypnozoite mode,
a round's radical cure takes hypnozoites out of the batches it names.

## Seasonality, vector control and importation

What comes from outside and what scales it: the season, local vector control, the three import
routes, and the two ways a case can be called imported.

### D36 — The seasonal factor peaks on the day its parameter names, and scales biting and the Poisson import rate

**Decision.** On day t of the run a mosquito bites with probability `BITING_RATE` × max(0, 1 +
`SEASONAL_AMP` cos(2π(t − `SEASONAL_PHASE`)/365)), times the vector-control factor inside its window
(D37) and `BITING_RATE_BOOST` from the resurgence trigger. The cosine puts the maximum exactly on
day `SEASONAL_PHASE`, and at every 365 days from it, so the parameter is what its table says it is:
the day of the year of the peak. `SEASONAL_AMP` is at least 0, since a negative amplitude would put
a trough on that day (a trough on day x is `SEASONAL_PHASE` = x + 182.5); above 1, the factor is 0
for part of each year, a legitimate edge value that runs (D42). The same factor multiplies the
Poisson import rate, `IMPORT_RATE_PER_TICK`; the pulse and the daily count take none (D37).

**Why.** A parameter should mean what it says. With a sine in place of the cosine the peak would
fall about 91 days, a quarter of a year, after `SEASONAL_PHASE`, and anyone setting it to a real
transmission-season day would be a quarter of a year out. The import rate follows the season because
travel, or transmission where the travellers pick up their infections, can follow the same season as
local biting.

**How to check.** `tests/testthat/test-params.R`: `vivax_params()`, which works the factor out as
the model does, puts the import-rate and biting peaks on `SEASONAL_PHASE` (day 100 at
`SEASONAL_PHASE = 100`, and day 3020, eight years on, once a boost starts on day 3000), and refuses
a negative `SEASONAL_AMP`, as the model does in a run started directly
(`tests/testthat/test-slim.R`); `tests/testthat/test-edge-values.R`: an amplitude of 1.5 runs to its
last day.

### D37 — Vector control scales biting only

**Decision.** On days `INTERVENTION_START` to `INTERVENTION_END` − 1 the daily biting probability is
multiplied by 1 − `INTERVENTION_EFFICACY`, together with the seasonal factor and the resurgence
boost (D36), and the factor scales nothing else: not importation, and not treatment, which has its
own window, `TREATMENT_START`/`TREATMENT_END` (D30). The Poisson import rate takes the seasonal
factor and `IMPORT_RATE_BOOST` from `RESURGENCE_TRIGGER_TICK`; the pulse (`N_IMPORTS`) and the daily
count (`N_IMPORTS_PER_TICK`) take no factor. A border measure is `IMPORT_RATE_BOOST`, which at 0
closes the border from the trigger, or the import window itself. `vivax_params()` works out the
import rate's peak for its warning (D40) with the same factors.

**Why.** Imports are infections acquired elsewhere, which local bed nets and spraying do not reach;
what limits them is a border measure, which the boost and the window express.

**Alternatives considered.** Scaling the Poisson import rate by the vector-control factor would tie
a border measure to local vector control, so an arm with vector control and open borders could not
be expressed at all. Scaling the pulse and the daily count by the season or the boost as well: they
are values the user sets directly, and a user who wants them to vary can vary them.

**How to check.** `tests/testthat/test-edge-values.R`: under full vector control from day 10 no bite
or uptake happens from day 10 while imports keep arriving, and `IMPORT_RATE_BOOST = 0` from the
trigger stops the imports; `tests/testthat/test-params.R`: an import rate of 40 a day into 30 hosts
warns at its full value under 90% vector control; `tests/testthat/test-treatment.R`: opening the
vector-control window gives no course; `inst/validation/check_parameters.R`: vector control lowers
prevalence.

### D38 — An import is minted in a scratch subpopulation in early(), and its second compartment is cloned the next day

**Decision.** Imports arrive in `early()`. Each import's outcomes, the same as a bite's, are rolled
before any genome exists: a blood-stage primary with probability `P_PRIMARY` and a hypnozoite batch
of geometric size with a bite's mean; an import that takes in neither is logged as
`vh_inoculation_import_no_take` and mints nothing. The day's imports that take are minted together
in one scratch subpopulation with `addSubpop()`, moved into their hosts with `takeMigrants()`, and
the scratch subpopulation is removed in the same callback. The migrated genome becomes whichever
compartment the import made, so an import that makes only hypnozoites arrives as its batch and is no
case (D29). An import that takes in both queues its batch, which the next day's `reproduction()`
clones from the primary and logs under the same `inoculation_id`. Both rows log
`parent1_pedigree_id` as -1, while the tree sequence records the primary as the batch's parent. An
imported case is detected and decided for treatment on arrival like any other (D31), and its course
lands no earlier than the next day's `early()` (D32), after the batch is cloned, so a radical cure
given for it reaches that batch. In the tree sequence the migrated genome's node keeps the scratch
subpopulation as its population, so a node's population does not name an import's host; nothing in
the package reads node populations.

**Why.** SLiM calls `reproduction()` once for each existing individual, so with no parasite alive it
never runs: an import seeded there could never reach an empty parasite population or re-seed a run
that reached elimination. Outside `reproduction()`, `addSubpop()` is the only way to mint a genome,
the idiom that also mints the founder lineages on day 1, and `takeMigrants()` keeps the pedigree id,
the tag and every stored value. SLiM refuses to clone a genome minted the same day, and two genomes
minted by `addSubpop()` are unrelated roots, so a clone on the next day is the only way to make an
import's batch IBD 1 to its primary: one inoculation delivers one sporozoite population, and a
relapse of an imported infection is then clonal to the infection that caused it, as a local relapse
is. Rolling the outcomes first leaves no genome with nowhere to go. The logged parent stays -1
because nothing inside the simulation is an import's parent, and `inoculation_id` and `lineage_id`
already say that the two rows are one inoculation; logging the primary as the batch's parent would
let anyone reconstructing transmission from the event log infer a blood-to-liver transmission within
the host that never happened.

**Alternatives considered.** A permanent sentinel individual that keeps `reproduction()` running in
an empty population: alive at the end of every run, it would be a tree-sequence sample with a
pedigree id in every all-samples analysis, and would need guards through the model's busiest
callback. A transient inoculum genome cloned into both compartments, which would give an import the
sibling topology of a local bite: it costs a one-day delay on every imported blood-stage infection
and an extra individual, for a topology that a structured source population (D39) gives as a
by-product.

**How to check.** `tests/testthat/test-model-invariants.R`: an inoculation's primary and its batch
are IBD 1 in the tree sequence, imports included; `inst/validation/audit_events_and_ibd.R`: each
`inoculation_id` is one host, a local inoculation's rows share a day, and an import's span at most
one; `tests/testthat/test-run.R`, among others, runs founder-free populations seeded only by
imports.

### D39 — Imported genotypes are unrelated founders

**Decision.** Every import enters as a root of the genealogy, with no ancestry: the limit of an
infinitely large, infinitely diverse source population, in which no two imports are related and no
import is related to the parasites it joins. `simulate_genotypes()` paints each import, as it paints
each founder template, independently from the same allele frequencies (D20), so an import's genotype
differs from local ones through its genealogy, never through differentiated allele frequencies.

**Why.** Unrelated founders add no parameter to the model, while a structured source would add the
source's diversity as a calibration axis and change what an import is. The cost is stated rather
than left implicit: imports are maximally distinguishable by construction, so any result on how
identifiable imported cases are is an upper bound, conditional on this assumption. The README's
limitations and sections 8, 11 and 15 of the model description say so.

**Alternatives considered.** A structured source deme with migration, the population-genetics
convention for immigration, under which imported genomes carry real ancestry and coalesce with one
another in the source: a planned extension (section 15 of the model description). Under it an
import's inoculum would be a real genome, so its two compartments would be siblings cloned from it,
as a local bite's are, in place of the primary-to-batch clone of D38.

**How to check.** `tests/testthat/test-genotypes.R`: the root genomes are the founder templates and
the imports, and distinct root genomes share no IBD segment.

### D40 — Imports land in uniformly random hosts under the bite's gate, at most one a host a day

**Decision.** A day's imports go to distinct hosts drawn uniformly at random, as bites do, so an
import can land in a host that already carries a clone. Each meets the gate every new inoculation
meets, bite or import: a host already carrying `MAX_MOI` blood-stage lineages, patent or subpatent,
counting those it gained earlier that day (D10), takes none, whether or not it is infectious, and a
refused import is logged as `vh_inoculation_import_no_take`. `MAX_MOI` is at least 1, where 1 means
no superinfection; relapses bypass the gate. A host takes at most one import a day, so a day's
imports beyond `N_HOSTS` are dropped, and `vivax_params()` warns when `N_IMPORTS`,
`N_IMPORTS_PER_TICK` or the two together would bring more than `N_HOSTS` on a day of the run, or
when the Poisson rate at its seasonal and boosted peak, with the daily count, would draw more than
`N_HOSTS` on a day with probability above 0.01. A pulse on day 1 lands like a pulse on any other
day, alongside the founders, with one difference: a founder host carrying a copy of another
founder's genome (`N_FOUNDER_LINEAGES < N_FOUNDER_HOSTS`) receives the copy on day 2, so on day 1 it
is still empty, takes an import even at `MAX_MOI = 1`, and from day 2 carries both. A naive
population is `N_FOUNDER_HOSTS = 0` with an import path inside the run.

**Why.** Random placement makes an import as likely as a bite to share its host with a local clone,
which is the setting a study detecting or classifying imported cases faces. `MAX_MOI` is a rule
about parasite burden, not about what a microscope or a mosquito can see, so one gate on every
blood-stage lineage present serves every new inoculation. One import a host a day keeps the largest
import rate the parameters allow, 1e6 a day, cheap to run, and the warning makes the dropped excess
visible. A day-1 pulse keeps the founders because the parameters ask for them; a naive population is
expressible directly.

**Alternatives considered.** Filling a day's imports from hosts without a blood-stage infection
first: an import would almost never meet a local clone, which makes imports easier to detect and
classify than random placement does. Placing every import with replacement: at the cap of 1e6
imports a day, each of the default 40 hosts would take about 25,000 a day. Gating a bite only on a
host with an infectious clone: two gates for one rule, since at `SUBPATENT_INFECTIVITY = 0` a host
whose clones are all subpatent would take bites past the cap while refusing imports. Reading a day-1
pulse as a naive population: it would silently drop the founders the parameters ask for.

**How to check.** `tests/testthat/test-run.R`: 100 imports into 200 hosts, half of them founders,
put about half into founder hosts; at `MAX_MOI = 1` none lands in a founder host still carrying its
clone, and the refusals are logged; a pulse of 50 into 30 hosts brings exactly 30, one a host, with
the warning; a bite never reaches a host at `MAX_MOI` whose clone is subpatent and not infectious; a
day-1 pulse keeps the founders. `tests/testthat/test-params.R`: `MAX_MOI = 0` is refused, as it is
in a run started directly (`tests/testthat/test-slim.R`), and a day's counts above `N_HOSTS`, or a
rate that often draws more, warn.

### D41 — An import class comes from the true chain, and a travel rule is scored against it

**Decision.** `classify_chain()` gives every case, a row of `classify_cases()`, its generation: the
fewest transmission links from a clone acquired by importation, through both routes (D52), found
breadth first from every import at once. A recrudescence counts as the clone it came from, and of an
outcrossed oocyst's two parents the shorter chain counts. The class follows the WHO malaria
terminology (World Health Organization, 2021): `imported` (generation 0: an import, a relapse of an
import's batch, or a recrudescence of either), `introduced` (1), `indigenous (from an import)` (2 or
more) and `indigenous (local lineage)` (no import in the chain); `import_pedigree_id` names the
import row the chain starts from, the batch for an imported relapse and the smaller id on a tie.
`classify_travel()` answers the other question, what a rule that knows only travel history
concludes: a case is `imported (recent travel)` when its host received an import that established,
as a blood stage or as hypnozoites only, on one of the `window` days that end on the case's day (60
by default), and `local (no recent travel)` otherwise. Each import is one trip, on the day it
arrived, although an import that takes in both compartments logs its batch the next day (D38); an
import that established nothing is no trip.

**Why.** The true class is what a travel rule, or a genetic classifier, is scored against. A
surveillance system records the journey, not which compartment the parasite ended in, so a traveller
who brings back only hypnozoites has travelled. The look-back decides most of the result, and for
vivax no window is right: a hypnozoite can wake many months after the trip, so a short window misses
many imported relapses, while one long enough to catch a hypnozoite that activates a year after the
trip also flags what its host acquires locally in that year; that disagreement is the quantity
`classify_travel()` exists to measure, what a travel rule costs, not a defect of it.

**Alternatives considered.** Walking each case's parents back one case at a time gives the same
classes, more slowly, and loses the chain at a recrudescence unless the walk takes it back to its
clone. Counting only the imports that established a blood stage as trips confuses the traveller with
the infection. Counting an import that established nothing as a trip: the model's only travel is the
import of a parasite, and travellers who bring none back are not modelled, so it would add a few of
them and leave out the rest.

**How to check.** `tests/testthat/test-imports.R`: a fixture's classes, generations and import rows,
with an imported relapse traced from its batch and an import's recrudescence; on a run, the
`imported` class is exactly the imported cases and each generation is one link from the one before;
the travel window's bounds, one trip per import on its arrival day, a hypnozoite-only import counted
and one that established nothing not; `tests/testthat/test-example-run.R`: the example run holds all
four classes.

## Parameters, outputs and the R interface

The surface a user meets: what is refused and what runs, how a run's tables read the same way every
time, the truth tables built from the event log, and the one function per question that reads them.

### D42 — Values the model cannot run are refused before SLiM starts; legitimate edge values run to the end

**Decision.** `validate_params()` refuses, before SLiM starts and with a message that names the
parameter, every value the model cannot run or would run as something else: a value that is missing,
not a single value, or of another kind than its default (a number, `TRUE`/`FALSE` or a string), or a
number with a class, such as a `bit64::integer64` or a lubridate duration; a count, `SIM_DAYS`, a
scheduling day or a fixed delay that is not whole; a value outside the bounds that `?vivax_params`
and section 10 of the model description give it, such as a probability outside [0, 1], a mean
duration below one day or a negative count; a number that is not finite, or not 0 and below 1e-300
in size; text that is not UTF-8 or holds control characters; a `CHROMOSOMES` that is unreadable,
names a chromosome outside 1 to 14 or twice, or runs a range high-low; and a few sizes, capped far
above anything realistic: a mean of 1000 oocysts per successful uptake, with `K_OOCYSTS` 0 or at
least 0.001, `LAMBDA_SPOROZOITES` at 1000, the hypnozoite batch mean at 1e6 (times its boost once
the trigger falls inside the run), the import rate at 1e6 a day, `TREATMENT_DELAY` and
`MEAN_RECRUD_DELAY` at 1e6 days, and `SIM_DAYS` at 1e9, SLiM's last tick. A run that reaches day
99999, the day the defaults use for never, must move those days past its end. A set in which no
parasite can ever arise, such as one with no founder and no import on a day of the run, and a
coverage set against a closed treatment window (D30) are refused too. The legitimate edges run at
their natural limit instead: no imports at a rate of 0, a one-day delay at the minimum, a
hypnozoite's death probability capped at 1 when activation and death together reach 1, no
infectiousness for a host without a blood stage, a zero batch from a mean of 0, a zero oocyst count
from a gamma draw that underflows to 0, and an exact draw, by inversion, of the oocyst count of an
uptake whose 10,000 gamma-Poisson tries all find none; a founder-free run whose imports never take
ends cleanly as never started (D5). Spellings that mean the same value, such as a named number, a
value in `I()`, a whole number carrying floating-point noise (`1.1 * 100`), a number or a factor
from `expand.grid()` for a text parameter, text in another encoding, or spaces around the commas and
dashes of `CHROMOSOMES`, are stored in the form SLiM is given (D7). The model makes the same checks
for a run started without R (D43).

**Why.** A sweep must not lose cells to engine errors, and a refusal that names the parameter beats
a stack trace from tick 200; each refused value would otherwise stop the run part-way or run it as
something else: a string `P_TRANS_VH` is compared as text, 12.5 founder hosts seed 13, and a
lubridate duration counts seconds. SLiM, and the provenance it writes into the tree sequence, need
UTF-8 text, and SLiM's command line cannot read a number between 0 and about 2.2e-308. A set in
which no parasite can arise would end with an empty event log, which in a sweep cannot be told from
a cell that saw no transmission. The caps sit far above any realistic setting and well below where
SLiM's draws fail: binomial sizes are 32-bit, Poisson draws wrap above 2^32 and the geometric delay
saturates near 4e9 days. The oocyst cap is on the mean count per uptake rather than on
`LAMBDA_OOCYSTS`, since a small `K_OOCYSTS` gives a tail whose mean is many times lambda; at
`K_OOCYSTS = 0.001` the count already comes within 0.5% of the largest variance a zero-truncated
negative binomial reaches at its mean, while the tries an uptake takes grow about as 1/K below it. A
cap allows a relative rounding of 1e-12, so a fitted value a few ulps above it (`uptake_calibrate()`
at a mean of 1000) runs. Day 99999 cannot be used in a run that reaches it because `run_sim()`
passes every parameter explicitly, so a day set to 99999 on purpose cannot be told from the default.
The exact draw keeps the count as dispersed as `uptake_moments()` says; an uptake reaches it with
probability P(0)^10000, 37% where a positive count is 1 draw in 10,000 (`LAMBDA_OOCYSTS` about 1e-4)
and below one in a million above a `LAMBDA_OOCYSTS` of about 0.003.

**Alternatives considered.** Clamping out-of-range values in the model: they would run as something
else, unannounced. Refusing the legitimate edges too: a closed border, a null arm without subpatent
infectivity and the minimum delay are real scenarios. Refusing the other spellings of a value: they
mean that value, and refusing them would cost sweep cells for nothing. Capping the import rate at
`N_HOSTS` a day: the model runs such rates, and the cap would refuse only small populations, so the
rate is capped where its draw would fail, and a day likely to bring more than `N_HOSTS` imports
warns.

**How to check.** `tests/testthat/test-params.R` refuses each impossible value with a message that
names the parameter, and stores the other spellings in SLiM's form;
`tests/testthat/test-edge-values.R` runs each legitimate edge to the end and checks that the run
reached the state its guard handles, with the minimum delays in `tests/testthat/test-treatment.R`;
`tests/testthat/test-model-invariants.R` checks the exact draw's mean and variance against
`uptake_moments()`; `tests/testthat/test-run.R` runs a set at the edge of the rounding allowance on
the caps.

### D43 — A run started with slim -d makes the same checks and takes either spelling of a number

**Decision.** The model's `initialize()` makes the checks of `validate_params()` (D42) again for a
run started directly with `slim -d`: each parameter's type and that it is a single value, whole
numbers within ±2^53 for the counts, days and seed, finite numbers, the probability ranges, the
founder bounds and the other single-value bounds, `GENOME_MODE`, `RADCURE_MODE`, `RECRUD_DELAY_DIST`
and `CHROMOSOMES`, a coverage against a closed treatment window, an overflowing biting rate, the
size caps, the days that mean never, and a set with no founder and no import path. Its caps allow a
relative rounding of 1e-11 rather than 1e-12, since the values arrive rounded to 15 digits, so a set
R accepts is never refused there; its messages give a float to 15 significant digits and quote a
string. Not repeated: a misspelt name, which SLiM defines as a constant the model never reads; the
limits on what text can carry (UTF-8, no control characters) and on numbers below 1e-300; and a
`SIM_DAYS` outside 2 to 1e9, not finite or not a single value, which SLiM refuses itself while it
reads the tick ranges, before `initialize()` runs, in a message that does not name `SIM_DAYS`.
Because `-d X=7` makes an Eidos integer where `-d X=7.0` makes a float, and Eidos will not mix the
two in `rbinom()`, `min()`, `seqLen()`, `setSeed()`, subsetting or a tick range, the model converts
each parameter to the type it needs (`asInteger()`, `asFloat()`) at each such use, so either
spelling gives the same run.

**Why.** A direct run is the one place where a typo cannot be caught in R, and each of these values
would otherwise start a run that means something else, as `RADCURE_MODE = 'Person'` would run the
per-hypnozoite mode and a probability above 1 is always true, or one that stops part-way.

**Alternatives considered.** Refusing one spelling, so that a float parameter must be written `7.0`:
a rule users would trip over. Converting each parameter once in `initialize()`: a constant given
with `-d` cannot be redefined in another type.

**How to check.** `tests/testthat/test-slim.R` refuses each case in a direct run with its parameter
named, and runs a five-day configuration with latency, an MDA round in the per-hypnozoite mode and a
pulse of imports in both spellings, with identical `events.tsv` and `state.tsv`.

### D44 — One parameter dictionary holds each parameter's facts, and an edited set is checked when it is used

**Decision.** `R/param-dictionary.R` holds one record per parameter: its typed default, unit,
meaning, single-value bounds, a joint or conditional rule in words, and its source. The defaults are
taken from it; roxygen generates the Parameters section of `?vivax_params` from it, and
`dev/update-doc-tables.R` the tables of section 10 of `MODEL_DESCRIPTION.md`; `param_info()` gives
it as a table, and `as.data.frame()` gives a parameter set's. `validate_params()` is the one check
(D42). A set edited with `$<-`, `[[<-` or `modifyList()` keeps its class and is checked in full when
it is used: `run_sim()` and `run_grid()` build every set again through `vivax_params()`, and
`p <- do.call(vivax_params, p)` checks one at once, as `?vivax_params` says. A default taken from
another parameter is set when a set is built, so editing `N_FOUNDER_HOSTS` leaves
`N_FOUNDER_LINEAGES` unchanged.

**Why.** A parameter's facts should be written once, where the code reads them, so that the help,
the model description and the checks cannot disagree; tests hold the documented bounds to the check,
the records' types to the model's, and the defaults to the model's own.

**Alternatives considered.** A `$<-` method that checks each edit: every caller that adds a value to
a set, `run_sim()` adding `SEED` and `OUTDIR` among them, and every test that builds a refused set
on purpose would have to strip the class first, for a check that use already makes. Generating
`validate_params()` from the dictionary: the joint rules are not single-value bounds, and each
refusal gives its rule's reason. Pointing section 10 at the help page: the model description is read
on its own.

**How to check.** `tests/testthat/test-param-dictionary.R`: a record and one group for each
parameter, the records' types against the model's, every documented bound against the one enforced,
each refusal naming its parameter, and `?vivax_params` and section 10 holding the current tables.
`tests/testthat/test-params.R`: the defaults against the model's `defineDefault()` values, each
evaluated as Eidos would, and an edited set refused by `run_sim()`.

### D45 — The loaders give every run the same columns and classes, and the host rows by default

**Decision.** `load_events()` and `load_state()` read their file against one schema per file
(`R/outputs.R`), the columns the model writes, in its order and classes, so every run gives the same
columns in the same order and classes whatever it holds: a run with no event gives a table with no
rows, a column the file lacks is added as `NA` of its class, a column the schema lacks follows in
the file's class, and a value that does not fit stops with the file and the column named, and its
line when it is text in a column of numbers. A missing or empty file, or an argument that is not a
run, is named in the error. The -1 that marks a field as not applicable is read as -1, not as `NA`.
`load_state()` gives the host rows unless `entity = "mosquito"` or `"all"` asks for the mosquitoes'
or both, in the file's order, and the columns and their classes do not depend on which rows are
asked for. The same schemas, with a meaning for each column, give `?events.tsv`, `?state.tsv` and
the tables of sections 16.1 and 16.2 of the model description.

**Why.** A table's columns are part of its contract, so a reader should not have to test classes
before using them; left to guess, `fread()` types a column by its values, and an empty file or a
column holding only `NA` would read as logical. Most rows of `state.tsv` are mosquitoes, 100 of
every 140 at the defaults, and a mosquito's `moi_active` counts the lineages it carries, so a
prevalence taken over every row mixes the two populations; the default is the rows nearly every
analysis wants, and no function of the package reads the mosquito rows.

**Alternatives considered.** Converting the -1 sentinels to `NA`: analyses test `!= -1`, and the
sentinel means not applicable, not missing. Every row of `state.tsv` by default: it keeps the
mistake easy. Separate loaders for hosts and mosquitoes: two functions for one file with one schema.

**How to check.** `tests/testthat/test-load.R`: an empty log, classes that do not depend on what a
run holds, missing and extra columns, a value that does not fit, a missing or empty file, normal,
extinct and never-started runs with the same columns and classes, each `entity` with its rows and
the same classes, and the help, the loaders and the model description's tables in agreement.

### D46 — The event dictionary is the authority, and the case, episode and recurrence tables are built from it

**Decision.** `event_types()` names every event type the model writes, with its `stage`, whether it
ends a patent episode (`ends_episode`) or begins one (`is_onset`), and a description. The package's
onset and episode-end filters are derived from those two flags, and its metrics take their cases
from `classify_cases()`, rather than restating event names, so a new way for an episode to begin or
end is picked up with no further edit. A type that cannot occur has no row: the founders are hosts
and every mosquito starts uninfected, so no event seeds a mosquito. `mosq_death` is described as
what it is, one row per genome an infected mosquito carried when it died; the deaths of uninfected
mosquitoes are not logged. Built on the dictionary, `classify_cases()` gives one row per clone
onset, founder seeding excluded, with its true origin and its `episode_id`, `"<human_id>:<tick>"`,
an episode being the clones that became patent in one host on one day; `episodes()` gives every
blood-stage clone's span from its onset to its first end, founders included; and `recurrences()`
pairs each clone of every episode after a host's first with each clone of the earlier episodes, or
of the latest one, never with a clone of its own episode, and reads each pair's IBD from the run's
IBD object when given it. Founder seeding is not an episode.

**Why.** A dictionary is read as a promise about the log, so a type that cannot occur, or a
description in the wrong unit, breaks it: counted as deaths, the `mosq_death` rows would miss every
uninfected mosquito's death and count an infected one once per genome. The unit a recurrence study
compares is the episode, and pairing row by row pairs the genotypes one bite delivered with each
other. One implementation of the spans and the pairings, checked against `state.tsv`, spares each
analysis a copy of its own, which could miss a way an episode begins or ends.

**Alternatives considered.** Seeding infected mosquitoes, so that a mosquito could start with a
founder genotype: no analysis asks for it, and the founders are hosts by design. Logging every
mosquito death: the rows of uninfected deaths would be most of a run's log, for a number the state
table already gives, since a mosquito's age starts again at 0 when it dies. One row per episode: the
genotypes of one episode differ in origin and in IBD, which a classifier needs. Founder seeding as
an earlier episode: every founder host's first acquired episode would become a recurrence, paired
with an unrelated founder genome.

**How to check.** `tests/testthat/test-load.R` (the dictionary, and what `mosq_death` counts),
`tests/testthat/test-run.R` (every type a run logs is in the dictionary) and
`tests/testthat/test-episodes.R` (the fixture's episodes, spans and pairings; with co-transmission,
several clones to an episode, never paired with each other; and on a run, every clone's span
covering exactly the days it counts in `moi_active`); `inst/validation/audit_events_and_ibd.R`
checks on each configuration that every type logged is documented and every clearance matches an
onset.

### D47 — state.tsv counts each host's infectious bites, which gives the EIR

**Decision.** Each host row of `state.tsv` carries `infectious_bites`, the bites the host received
that day from a mosquito carrying a product past the EIP, whatever the bite then did (`NA` on
mosquito rows). The model counts a bite as infectious exactly where it already draws the mosquito's
transmissible products, which come back empty exactly when no product has completed the EIP, so the
count draws no random number. The count is kept on the host and reset after each day's snapshot,
since `reproduction()` does not run on a day with no parasite. `summarise_run()` reports the EIR
from it as `eir_per_py`, infectious bites per person-year of every host over its window.

**Why.** The EIR is the usual measure of transmission intensity and, with prevalence, the quantity
to calibrate `BITING_RATE` against (section 13 of the model description). Counting at the bite is
exact, and per host it also gives the heterogeneity of exposure. The Ross–Macdonald formula, with
the sporozoite rate read from `n_transmissible`, gives about 1 − 1/`MOSQ_LIFESPAN` of the counted
EIR, about 0.9 at the default 10 days, because the snapshot is taken after the day's mosquito
deaths, and with seasonality, vector control or a boost its error can go either way (section 4.1 of
the model description).

**Alternatives considered.** A column on the mosquito rows: the EIR is per person, and the host row
is where a reader looks for exposure. A file of its own: one more output to load and document, for
one number per host and day.

**How to check.** `tests/testthat/test-model-invariants.R`: where every infectious bite infects and
no host is full, with case management on, each host's infectious bites a day are its inoculations
that day. `tests/testthat/test-summary.R`: `eir_per_py` is the window's infectious bites over its
person-years.

### D48 — A table's shape follows its arguments, never the run's content

**Decision.** Every reader returns the same columns, in the same order and classes, with the same
factor levels, on every path, whatever the run holds. `force_of_infection()` splits relapses into
`relapse (imported)` and `relapse (local)` unless `split_imported = FALSE`, whether or not the run
holds imports, and returns `origin` as a factor with all its levels; `rc_by_origin()` gives every
origin a row; `classify_travel()` gives a run without imports the same columns; the IBD tables keep
one column order and class whether a pair shares everything or nothing; a run with no parasite gives
every reader its usual columns, with no rows where there is nothing to list, counts of 0 and `NA`
where a quantity needs a parasite; and a table of the wrong kind is refused by name.
`new_infection_diversity()`, which caps the genomes a window compares, draws them with its own
`seed`, the k-th window with `seed + k - 1`, always under the same RNG kinds, and leaves the
caller's random stream untouched.

**Why.** A grid pools per-run tables, so a shape that varies with a run's content either fails on
one cell or, worse, binds by position and gives numerically plausible nonsense. An imported relapse
is a distinct truth category, so relapses are split by provenance, in every run. A seeded draw lets
repeat calls on one run agree, as a metric must, while the draw is unbiased in time and lineage.

**Alternatives considered.** Splitting relapses only in runs that hold imports: the same quantity
would come back as `relapse` in one run and `relapse (local)` in another, and a pooled `dcast()`
would put one origin in three columns without an error. Taking a window's first genomes instead of a
seeded draw: reproducible, but biased toward the window's early genomes.

**How to check.** `tests/testthat/test-windows.R` (the same `origin` levels for a run with no case),
`tests/testthat/test-ibd.R` (segment tables keep one schema), `tests/testthat/test-no-parasite.R`
(every reader takes a run with no parasite), `tests/testthat/test-inputs.R` (a reader given the
other kind of table says which it got) and `tests/testthat/test-rng.R` (a seeded draw is the same
under any RNG kind, and leaves the caller's kinds and stream alone).

### D49 — Windowed series have a row for every window, divided by the days it covers, to SIM_DAYS

**Decision.** `force_of_infection()` gives one row per window of the run and per origin, with
`n = 0` where there is no case, and `days`, the days the window covers within the run: the first
window holds days 1 to `bin - 1`, the last ends at `SIM_DAYS`, and `rate` divides by `days`. `bin`
is a whole number of days, at least 1, and `new_infection_diversity(bin =)` gives every window of
the same grid a row. Every function that summarises a window of a run takes the run's end from
`SIM_DAYS`, whether or not the run went extinct (D5), so the days after an extinction count as days
with no parasite: prevalence and the reservoir are means over every host-day of the window, while
COI and polyclonality are means over the infected host-days, `NA` when there are none. The per-day
tables, `coi_distribution()` and `host_state_prevalence()`, report the days `state.tsv` holds.

**Why.** A window without a case is a measurement of zero, not a gap that a line or a stacked area
draws straight across, and a short window divided by `bin` would start or end every series in a dip
that is pure binning artefact; full windows are the same either way, since there `days` is `bin`. A
run's window is the period the study asks about, not the period the parasite lasted: a run that died
out had no infection on its remaining days.

**Alternatives considered.** Dropping the partial windows: it throws away a real estimate, unbiased
though noisier. Zero-filling only in the plots: every other reader of the table would meet the gaps
again. Ending a run's series at its last state row or its last case: an extinct run's windows would
be shorter, and a trailing window whose true value is zero, the most informative reading there is,
would be lost. Padding `state.tsv` to `SIM_DAYS`: the file would hold rows for days the model never
ran.

**How to check.** `tests/testthat/test-windows.R`: every window and origin, the days of the first
and last windows, a run with no case and its plot, an extinct run's zeros after its extinction, and
the window statistics of runs that went extinct before day 200 of 400.

### D50 — summarise_run() is the one summary of a run, over one window of days

**Decision.** `summarise_run()` gives one row of typed, unrounded columns named by their units
(`_frac`, `_per_py`, `_tick`, `n_`), every one over the same window: days
`round(SIM_DAYS * burn_frac) + 1` to `SIM_DAYS`, the second half by default, or `from` and `to`
directly. The row holds prevalence, COI and polyclonality over the window's host-days, the
reservoir, the EIR (D47), the FOI and both molFOB readings, by clone and by episode (Koepfli et al.
2013; D33), the relapse rate per person-year at risk and both relapse fractions (D28), the imported
share, the counts, and the person-years and person-years at risk; the EIR, FOI and molFOB divide by
the person-years of every host. It records its window, its reading (`exclude_masked`, which leaves
the masked relapses of D26 out of every relapse column, while `n_relapse_masked` counts them either
way) and how the run ended. A quantity with nothing to divide is `NA`, and a run that stopped early
with no extinction or start-failure record is refused unless `to` is at most its last day.
`plot_diagnostics()` takes its subtitle's numbers and its burn-in line from the same row, and
`new_infection_diversity()` starts on the same first day for a given `burn_frac`.

**Why.** A summary pooled over a grid must mean the same thing in every column and every cell, and
the counts and denominators on the row let runs be pooled by summing them rather than by averaging
rates. One definition serves the user, the dashboard and the calibration:
`inst/validation/calibrate_nu.R` fits `NU_HYPNOZOITE_BATCH` to `summarise_run()`'s own columns.

**Alternatives considered.** Summing `force_of_infection()`'s rows for the FOI and molFOB: the
window need not line up with its bins. A class with a print method: `rbindlist()` drops it. One
column per masked reading: the reading is an argument, and the row records it.

**How to check.** `tests/testthat/test-summary.R`: names and classes on every kind of run; the
fixture's values over the whole run, the default window and two others; `exclude_masked`; extinct,
never-started and unfinished runs; and on a SLiM run, every column against its definition over the
window, FOI and molFOB against `force_of_infection()`, and the dashboard's subtitle against the row.

### D51 — host_state_prevalence() has every state on every day and counts subpatent carriers

**Decision.** `host_state_prevalence()` gives one row per day of the state table and per state, with
`n = 0` where no host is in it and `frac` its share of the day's hosts, in five levels:
`uninfected`, `active only`, `active + reservoir`, `subpatent` (a subpatent clone and no patent one,
with or without hypnozoites) and `reservoir only (invisible)`, the hosts with hypnozoites alone,
whom no blood survey sees. Each row also carries the day's mean hypnozoites per host, so the
reservoir's trajectory comes off the same table. At the default `P_RECRUD_GIVEN_FAIL = 0` no clone
is subpatent, and that state is 0 throughout.

**Why.** A subpatent carrier is infected and can infect a mosquito, so counting it as uninfected
would understate carriage exactly where treatment fails. A missing row is not a zero: without the
zero rows a stacked area joins the rows on either side of an empty state and draws hosts on days
that had none.

**Alternatives considered.** Putting subpatent carriers with hypnozoites in the reservoir-only
state: they carry a blood-stage clone, which a sensitive assay would find.

**How to check.** `tests/testthat/test-host-states.R` (every state on every day, the subpatent
carriers, an empty table, the reservoir on each row, the plot) and `tests/testthat/test-treatment.R`
(on a run with recrudescence, the subpatent count is each day's hosts with a subpatent clone and no
patent one).

### D52 — Transmission links follow both routes, with half weight to each parent of an outcrossed oocyst, and form the chains

**Decision.** `transmission_links()` joins each parent genome of an oocyst to each genome the
oocyst's products infected: `direct`, the bite's blood stage, or `liver`, each relapse of the
hypnozoite batch the same bite left, which carries the batch's oocyst. Only local infections are
targets: imports, founders, the relapses of an import's batch and recrudescences can only be
sources. A host's own relapse is never linked to its earlier infection, though a mosquito that feeds
on a host and later bites it makes a link like any other. Each row names both genomes and hosts, the
mosquito, the oocyst, four days (the source's onset, the feed, the bite and the target) and a
weight: 1, or 1/2 for each parent of an outcrossed oocyst, so a target's weights sum to 1. `kind`
keeps either route, and a log cut short of the rows a link needs gives one warning. A chain is a
connected component of the links taken in either direction, each source taken as the clone it came
from through repeated recrudescences; each link carries its `chain_id`, numbered in the order of the
chains' first links, and `chain_size`, the chain's genomes, founders included. A genome with no link
is in no chain. The components come from a union-find with path halving.

**Why.** A relapse in another host is a transmission delayed by the liver, and relapses are about
three in five new blood-stage clones in the model description's representative run, so leaving liver
links out would undercount onward transmission most for the infections whose batches relapse. A
recrudescence is the same infection as the clone it came from (D33), so its onward transmission
belongs to that clone's chain. On the links themselves a chain can be followed and joined, and the
union-find stays linear in practice where liver links make chains of thousands of genomes.

**Alternatives considered.** One row per target: it loses the second parent of an outcrossed oocyst.
Carrying a target's `detected` or `into_patent_lineage` on its link: `classify_cases()` holds a
case's attributes, and the help gives the join. Label propagation with pointer jumping for the
chains: it hooks nodes, not roots, so it needs about as many passes as a chain is long. A chain of
its own for each unlinked genome: an edge table cannot hold it, and `classify_chain()` labels every
case with its introduction (D41).

**How to check.** `tests/testthat/test-transmission.R`: the fixture's links, targets and weights, a
self link, `kind`, empty and cut logs, the chains of each kind with a recrudescence kept with its
clone, the components against a reference union-find on long, disjoint, star and random graphs, and
on a run, every local case the target of one transmission and the chains partitioning the linked
genomes. `inst/validation/audit_events_and_ibd.R`: its transmission-link checks on each
configuration.

### D53 — Rc is the case reproduction number by true origin, each secondary case counted once

**Decision.** `rc_by_origin()` gives the mean number of secondary cases each case caused (Fraser
2007), by the case's true origin. A case is a clone onset, a row of `classify_cases()`; its
secondary cases are the targets of its links through both routes (D52), each counted once, so one
from an outcrossed oocyst credits 1/2 to each parent. `kind` chooses the routes, and `direct` and
`liver` add up to `all`. The function takes the run: `from` and `followup` choose the onsets
counted, from day `from` to `SIM_DAYS - followup`, both 0 by default. With `unit = "acquisition"` a
recrudescence is not a case, and its onward transmission is credited to the clone it re-expanded
from, through repeated recrudescences. Every origin keeps a row, with `n_cases`, `n_secondary`, `Rc`
and `onward_frac` unrounded, and `Rc` and `onward_frac` `NA` where there is no case.

**Why.** A case reproduction number counts what each case actually caused, which the links record
exactly. Half credit counts each secondary case once, so that Rc sits near 1 at a self-sustaining
steady state, above 1 while an epidemic grows from its founders, and lower where imports sustain
transmission. A case near the end of a run has had less time to transmit, above all through
relapses: a hypnozoite stays dormant 1/(`ALPHA_ACTIVATE` + `MU_HYPNOZOITE`) days on average, about
four months at the defaults, and more than a year for 1 in 20.

**Alternatives considered.** Full credit to both parents of an outcrossed oocyst: it counts such a
secondary case twice. A default follow-up: the right one depends on the run's phase, which the help
explains. An events table and a separate run length as arguments: the run gives both.

**How to check.** `tests/testthat/test-transmission.R`: the fixture's values for both units, half
credit and the recrudescence walked back; the routes adding up; a run with no case; `from` and
`followup` and their refusals; and on a run, the counts as the weights of the links from cases.

### D54 — One function per question, every reader takes a run, and no plumbing is exported

**Decision.** Each export answers one question, and a function whose own help would tell the user
not to use it is not exported. Plumbing is internal: SLiM's path and installed version, which
`check_setup()` shows; the builder of SLiM's arguments, which `run_sim()` uses; and the plot theme,
which every plot applies. Every reader takes a run, a `vivax_run` or its directory, and a plot takes
a run or, where it draws one function's table, that table; a table of the other kind is refused by
name. One naming rule holds throughout: `from` and `to` for a window, `n` for a count, `_frac` for a
share from 0 to 1, and every table a reader returns a `data.table`. `example_run()` returns a small
run that ships with the package, 60 hosts for 400 days with imports, treatment, recrudescence and
co-transmission, so that every kind of case, recurrence and link occurs in it;
`data-raw/example_run.R` makes it from the settings in `R/example.R`.

**Why.** A user looking for how to answer a question should find one function, and plumbing is only
more to learn: a methods section needs the SLiM version a run ran on, which every run records
(`load_run()$slim_version`), and a user restyles a plot with ggplot2's own themes. One convention
lets any reader be pointed at a run. The bundled run lets the examples, the vignettes and a first
look run without SLiM, and the tests that compare it with its settings and with the model fail when
a default or the model changes, until it is made again.

**How to check.** `tests/testthat/test-exports.R` lists the exports exactly;
`tests/testthat/test-inputs.R` gives every reader a run and every plot a run or its table, and each
the wrong kind of table; `tests/testthat/test-docs.R` holds every help page to one description of a
run; `tests/testthat/test-example-run.R` checks that the shipped run is this model's run of its
settings and holds every kind of case, recurrence and link.
