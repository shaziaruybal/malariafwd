# Each case as a travel rule classes it: recent travel or not

What a surveillance system that knows only travel history concludes: a
case is `"imported (recent travel)"` when its host brought an infection
in during the `window` days up to the case, and
`"local (no recent travel)"` otherwise. Score it against the truth, the
`imported` column of
[`classify_cases()`](https://shaziaruybal.github.io/malariafwd/reference/classify_cases.md)
or the `who_class` of
[`classify_chain()`](https://shaziaruybal.github.io/malariafwd/reference/classify_chain.md).
The two are expected to disagree: the disagreement is what a travel rule
costs.

## Usage

``` r
classify_travel(events, window = 60)
```

## Arguments

- events:

  an events table from
  [`load_events()`](https://shaziaruybal.github.io/malariafwd/reference/load_events.md),
  or a run (a `vivax_run` or its directory).

- window:

  the look-back, a whole number of days, 0 or more: a case on day `tick`
  is flagged when an import arrived in its host on a day from
  `tick - window + 1` to `tick`. An import's own case is flagged for any
  `window` of 1 or more; at 0 no case is.

## Value

a data.table, one row per row of
[`classify_cases()`](https://shaziaruybal.github.io/malariafwd/reference/classify_cases.md)
and in its order: its columns, then `recent_travel` (logical) and
`surveillance`, `"imported (recent travel)"` or
`"local (no recent travel)"`. A run without imports gives the same
columns, every case local.

## Details

The trip is recorded once, on the day an import arrives, whether it
established a blood-stage infection, hypnozoites or both. An import that
established nothing (`vh_inoculation_import_no_take`) is not a trip: the
model's only travel is the import of a parasite, and travellers who
bring none back are not modelled, so counting the imports that failed
would add a few of them and leave out the rest.

The window decides most of the result, and for vivax no window is right:
a relapse of an imported batch can come long after the trip, so a short
window misses the late relapses, while a long one flags more of the
cases the host acquired locally after the trip.

## See also

[`classify_chain()`](https://shaziaruybal.github.io/malariafwd/reference/classify_chain.md),
[`classify_cases()`](https://shaziaruybal.github.io/malariafwd/reference/classify_cases.md)

## Examples

``` r
tr <- classify_travel(example_run(), window = 90)
table(tr$surveillance, tr$imported)   # the travel rule against the truth
#>                           
#>                            FALSE TRUE
#>   imported (recent travel)    28   43
#>   local (no recent travel)   221    7
```
