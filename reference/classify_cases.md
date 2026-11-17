# Classify each blood-stage clone onset ("case") by its true origin

One row per clone that became patent in a host, in the order of the
event log: a bite or an import that gave a blood-stage primary, a
relapse, or a recrudescence. Clones that became patent in one host on
one day are one episode, sharing an `episode_id` (`"<human_id>:<tick>"`,
the same in every function and subset); a bite that delivered two
genotypes gives two rows of one episode. Founder seeding
(`vh_inoculation_seed`) is not a case: the founders are where the run
starts. Detection is drawn per bite, activation, recrudescence or
import, so the clones of one episode can differ in `detected`.

## Usage

``` r
classify_cases(events, exclude_masked = FALSE)
```

## Arguments

- events:

  an events table from
  [`load_events()`](https://shaziaruybal.github.io/malariafwd/reference/load_events.md),
  or a run (a `vivax_run` or its directory).

- exclude_masked:

  if `TRUE`, leave out the masked relapses.

## Value

a data.table, one row per clone onset: `tick` (the day it became
patent), `human_id`, `episode_id`, `lineage_id`, `oocyst_id`,
`pedigree_id`, `origin`, `imported`, `detected` (1 or 0) and
`into_patent_lineage` (`TRUE` for a masked relapse, `FALSE` for every
other row).

## Details

`origin` gives the mechanism by which the clone arose, and takes
**four** values: `relapse` = a hypnozoite activated; `imported` = a new
infection that arrived via importation; `local reinfection` = a new
local mosquito-borne infection; and `recrudescence` = a clone that
survived a treatment course, was suppressed below patency, and
re-expanded. Note the fourth is **not** an acquisition – no new parasite
entered the host – so it is excluded from
[`force_of_infection()`](https://shaziaruybal.github.io/malariafwd/reference/force_of_infection.md)'s
rates and from
[`new_infection_diversity()`](https://shaziaruybal.github.io/malariafwd/reference/new_infection_diversity.md).
`imported` is a separate genotype-origin flag that is `TRUE` whenever
the parasite genotype is import-derived – including **imported
relapses** (a relapse of an imported genotype), a key vivax truth
category that `origin` alone ("relapse") does not expose. So the four
truth classes are: imported new infection (`origin == "imported"`),
local reinfection, imported relapse (`origin == "relapse"` &
`imported`), and local relapse (`origin == "relapse"` & `!imported`).

A relapse into a patent clone of its own lineage, born on an earlier
day, is *masked*: a study that genotypes the host sees no new clone
(section 5 of the [model
description](https://shaziaruybal.github.io/malariafwd/MODEL_DESCRIPTION.html)).
`into_patent_lineage` marks it, and `exclude_masked = TRUE` leaves it
out.

## See also

[`episodes()`](https://shaziaruybal.github.io/malariafwd/reference/episodes.md)
for every clone's span,
[`recurrences()`](https://shaziaruybal.github.io/malariafwd/reference/recurrences.md)
for each episode against the earlier ones.

## Examples

``` r
cs <- classify_cases(example_run())
table(cs$origin, cs$imported)
#>                    
#>                     FALSE TRUE
#>   imported              0   23
#>   local reinfection   113    0
#>   recrudescence         4    1
#>   relapse             132   26
```
