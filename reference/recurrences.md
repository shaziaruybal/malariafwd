# Each recurrence, against the earlier episodes of its host

A recurrence is every episode of a host after its first: a relapse, a
reinfection, an import or a recrudescence. Each of its clones is paired
with each clone of the earlier episodes (`prior = "all"`) or of the
latest one (`"previous"`); the clones of one episode are never paired,
and founder seeding is not an episode. Whether the earlier clone was
still patent is not a condition: join
[`episodes()`](https://shaziaruybal.github.io/malariafwd/reference/episodes.md)
on `prior_pedigree_id` for that. Give the run's IBD object to read each
pair's IBD.

## Usage

``` r
recurrences(
  events,
  prior = c("all", "previous"),
  from = 0,
  ibd = NULL,
  exclude_masked = FALSE
)
```

## Arguments

- events:

  an events table from
  [`load_events()`](https://shaziaruybal.github.io/malariafwd/reference/load_events.md),
  or a run (a `vivax_run` or its directory).

- prior:

  `"all"` (every earlier episode) or `"previous"` (the latest one).

- from:

  keep the recurrences that begin on day `from` or later, paired with
  earlier episodes whenever they were.

- ibd:

  `NULL`, or the run's IBD object from
  [`ibd_pairs()`](https://shaziaruybal.github.io/malariafwd/reference/ibd_pairs.md),
  for an `ibd` column.

- exclude_masked:

  if `TRUE`, masked relapses (into a patent clone of their own lineage,
  see
  [`classify_cases()`](https://shaziaruybal.github.io/malariafwd/reference/classify_cases.md))
  are neither recurrences nor earlier episodes.

## Value

a data.table, one row per (recurrence clone, earlier clone), by host,
day and pedigree id: `human_id`, `episode_id`, `tick`, `pedigree_id`,
`lineage_id`, `origin`, `imported`, `detected`, `into_patent_lineage`,
the same of the earlier clone prefixed `prior_` (`prior_episode_id`,
`prior_tick`, `prior_pedigree_id`, `prior_lineage_id`, `prior_origin`,
`prior_imported`, `prior_detected`), `days_since` (1 or more), and with
`ibd` given, `ibd` (`NA` where either genome is not in the object).

## Details

With `prior = "all"` the table grows with the square of a host's
episodes. To pair each recurrence with the latest earlier episode on the
observed line list, filter `recurrences(ev)` to `prior_detected == 1`
and keep each clone's row with the largest `prior_tick`.

## See also

[`episodes()`](https://shaziaruybal.github.io/malariafwd/reference/episodes.md),
[`classify_cases()`](https://shaziaruybal.github.io/malariafwd/reference/classify_cases.md),
[`ibd_pairs()`](https://shaziaruybal.github.io/malariafwd/reference/ibd_pairs.md)

## Examples

``` r
head(recurrences(example_run()))
#>    human_id episode_id  tick pedigree_id lineage_id            origin imported
#>       <int>     <char> <int>       <int>     <char>            <char>   <lgcl>
#> 1:        2      2:237   237         596          2           relapse    FALSE
#> 2:        2      2:284   284         766          7 local reinfection    FALSE
#> 3:        2      2:284   284         766          7 local reinfection    FALSE
#> 4:        2      2:298   298         828          2           relapse    FALSE
#> 5:        2      2:298   298         828          2           relapse    FALSE
#> 6:        2      2:298   298         828          2           relapse    FALSE
#>    detected into_patent_lineage prior_episode_id prior_tick prior_pedigree_id
#>       <int>              <lgcl>           <char>      <int>             <int>
#> 1:        1                TRUE            2:230        230               577
#> 2:        1               FALSE            2:230        230               577
#> 3:        1               FALSE            2:237        237               596
#> 4:        0               FALSE            2:230        230               577
#> 5:        0               FALSE            2:237        237               596
#> 6:        0               FALSE            2:284        284               766
#>    prior_lineage_id      prior_origin prior_imported prior_detected days_since
#>              <char>            <char>         <lgcl>          <int>      <int>
#> 1:                2 local reinfection          FALSE              1          7
#> 2:                2 local reinfection          FALSE              1         54
#> 3:                2           relapse          FALSE              1         47
#> 4:                2 local reinfection          FALSE              1         68
#> 5:                2           relapse          FALSE              1         61
#> 6:                7 local reinfection          FALSE              1         14
run <- example_run()
x <- ibd_pairs(run, ids = classify_cases(run)$pedigree_id)
head(recurrences(run, ibd = x))   # with each pair's IBD
#>    human_id episode_id  tick pedigree_id lineage_id            origin imported
#>       <int>     <char> <int>       <int>     <char>            <char>   <lgcl>
#> 1:        2      2:237   237         596          2           relapse    FALSE
#> 2:        2      2:284   284         766          7 local reinfection    FALSE
#> 3:        2      2:284   284         766          7 local reinfection    FALSE
#> 4:        2      2:298   298         828          2           relapse    FALSE
#> 5:        2      2:298   298         828          2           relapse    FALSE
#> 6:        2      2:298   298         828          2           relapse    FALSE
#>    detected into_patent_lineage prior_episode_id prior_tick prior_pedigree_id
#>       <int>              <lgcl>           <char>      <int>             <int>
#> 1:        1                TRUE            2:230        230               577
#> 2:        1               FALSE            2:230        230               577
#> 3:        1               FALSE            2:237        237               596
#> 4:        0               FALSE            2:230        230               577
#> 5:        0               FALSE            2:237        237               596
#> 6:        0               FALSE            2:284        284               766
#>    prior_lineage_id      prior_origin prior_imported prior_detected days_since
#>              <char>            <char>         <lgcl>          <int>      <int>
#> 1:                2 local reinfection          FALSE              1          7
#> 2:                2 local reinfection          FALSE              1         54
#> 3:                2           relapse          FALSE              1         47
#> 4:                2 local reinfection          FALSE              1         68
#> 5:                2           relapse          FALSE              1         61
#> 6:                7 local reinfection          FALSE              1         14
#>      ibd
#>    <num>
#> 1:     1
#> 2:     0
#> 3:     0
#> 4:     1
#> 5:     1
#> 6:     0
```
