# Every blood-stage clone of a run, from onset to end

One row per blood-stage clone, founder seeding included: when it became
patent, how, and when and how it ended. A clone counts in `state.tsv`'s
`moi_active` from `onset_tick` to `end_tick - 1`, so a clone that ends
on its first day (an MDA round can clear one) is never in `state.tsv`. A
clone that ends below patency (`bloodcure_subpatent`) returns, if it
does, as a recrudescence: a row of its own.

## Usage

``` r
episodes(events)
```

## Arguments

- events:

  an events table from
  [`load_events()`](https://shaziaruybal.github.io/malariafwd/reference/load_events.md),
  or a run (a `vivax_run` or its directory).

## Value

a data.table, by host, onset and pedigree id: `episode_id`, `human_id`,
`pedigree_id`, `lineage_id`, `onset_tick`, `onset_event`
(`vh_inoculation`, `vh_inoculation_import`, `vh_inoculation_seed`,
`hypnozoite_activation` or `recrudescence`), `end_tick` (`NA` for a
clone still patent when the run ended) and `end_reason`
(`human_clearance`, `case_bloodstage_clear`, `mda_bloodstage_clear`,
`bloodcure_subpatent`, or `"ongoing"`). Join
[`classify_cases()`](https://shaziaruybal.github.io/malariafwd/reference/classify_cases.md)
on `pedigree_id` for each clone's origin; founders have none.

## See also

[`recurrences()`](https://shaziaruybal.github.io/malariafwd/reference/recurrences.md),
[`classify_cases()`](https://shaziaruybal.github.io/malariafwd/reference/classify_cases.md)

## Examples

``` r
ep <- episodes(example_run())
head(ep)
#>    episode_id human_id pedigree_id lineage_id onset_tick           onset_event
#>        <char>    <int>       <int>     <char>      <int>                <char>
#> 1:        0:1        0           0          0          1   vh_inoculation_seed
#> 2:        1:1        1           1          1          1   vh_inoculation_seed
#> 3:      1:340        1        1215          7        340        vh_inoculation
#> 4:        2:1        2           2          2          1   vh_inoculation_seed
#> 5:      2:230        2         577          2        230        vh_inoculation
#> 6:      2:237        2         596          2        237 hypnozoite_activation
#>    end_tick      end_reason
#>       <int>          <char>
#> 1:       73 human_clearance
#> 2:      107 human_clearance
#> 3:      353 human_clearance
#> 4:      213 human_clearance
#> 5:      279 human_clearance
#> 6:      261 human_clearance
```
