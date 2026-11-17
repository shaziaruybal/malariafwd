# Classify each case by its transmission chain from an import (the WHO classes)

Each blood-stage case, a row of
[`classify_cases()`](https://shaziaruybal.github.io/malariafwd/reference/classify_cases.md),
with the class the WHO malaria terminology gives it by where it was
acquired and its chain of transmission from an imported infection:
`"imported"` (acquired abroad: an import, a relapse of an import's
batch, or a recrudescence of either), `"introduced"` (infected locally
by an imported case), `"indigenous (from an import)"` (two or more links
from one) and `"indigenous (local lineage)"` (no import in its chain).

## Usage

``` r
classify_chain(events)
```

## Arguments

- events:

  an events table from
  [`load_events()`](https://shaziaruybal.github.io/malariafwd/reference/load_events.md),
  or a run (a `vivax_run` or its directory).

## Value

a data.table, one row per row of
[`classify_cases()`](https://shaziaruybal.github.io/malariafwd/reference/classify_cases.md)
in its order, with its columns and `from_import` (logical), `generation`
(the fewest links from a clone acquired by importation, 0 for an
imported case, `NA` with no import in the chain), `import_pedigree_id`
(the genome of the import row the shortest chain starts from: an
imported infection's blood-stage genome, or the batch of an imported
relapse; the smaller id on a tie; `NA` with no import) and `who_class`.

## Details

The chain follows
[`transmission_links()`](https://shaziaruybal.github.io/malariafwd/reference/transmission_links.md)
through both routes, so a relapse is linked to the infector of the bite
that laid its batch. A recrudescence is the same infection as the clone
it came from. An outcrossed oocyst has two infectors: a case is from an
import when either is, and the shorter chain counts. Masked relapses are
classed like any case; `into_patent_lineage` leaves them out.

## References

World Health Organization (2021) WHO malaria terminology, 2021 update.
Geneva: World Health Organization.

## See also

[`classify_travel()`](https://shaziaruybal.github.io/malariafwd/reference/classify_travel.md),
[`transmission_links()`](https://shaziaruybal.github.io/malariafwd/reference/transmission_links.md),
[`classify_cases()`](https://shaziaruybal.github.io/malariafwd/reference/classify_cases.md)

## Examples

``` r
ch <- classify_chain(example_run())
table(ch$who_class)
#> 
#>                    imported indigenous (from an import) 
#>                          50                          38 
#>  indigenous (local lineage)                  introduced 
#>                         188                          23 
```
