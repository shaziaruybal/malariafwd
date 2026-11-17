# Pairwise IBD segments along the genome

The stretches of genome over which two parasites are identical by
descent – the detail behind the single IBD number from
[`ibd_pairs()`](https://shaziaruybal.github.io/malariafwd/reference/ibd_pairs.md).
A segment is a maximal stretch over which the two genomes descend from
one common ancestor genome, cut at each chromosome junction, the same
rows as `ibd_pairs(segments = TRUE)` gives for the pair; each is in both
concatenated and within-chromosome coordinates.

## Usage

``` r
ibd_segments(run, ids_a, ids_b, min_span = 1000)
```

## Arguments

- run:

  a run: a `vivax_run` (from
  [`run_sim()`](https://shaziaruybal.github.io/malariafwd/reference/run_sim.md)
  or
  [`load_run()`](https://shaziaruybal.github.io/malariafwd/reference/load_run.md))
  or the path of a run directory.

- ids_a, ids_b:

  equal-length vectors of pedigree ids, compared position by position.

- min_span:

  the shortest segment returned, in bp, after the cut at the junctions
  (a segment of exactly `min_span` is kept): below it a shared stretch
  is treated as noise, mirroring the detection floor of IBD inference.

## Value

a data.table, by pair and position: `pair` (the position in `ids_a` and
`ids_b`, from 1), `a` and `b` (as given), `chrom` (the PvP01
chromosome), `start` and `end` (bp within the chromosome, `end`
exclusive), `gstart` and `gend` (bp along the concatenated genome),
`length` (bp), `mrca_id` (the pedigree id of the common ancestor),
`mrca_time` and `depth` (its age; see
[`ibd_pairs()`](https://shaziaruybal.github.io/malariafwd/reference/ibd_pairs.md)).
A pair sharing nothing has no rows, and nor does a pair with a genome
the tree sequence does not hold, which warns; the table's `"missing"`
attribute lists those ids. A genome with itself is shared end to end:
one row per chromosome, its own ancestor at depth 0. A pair given twice
has its rows twice.

## Python

The IBD and ARG functions read a run's tree sequence with Python's
tskit, through the reticulate package. Either of two routes works.
*Automatic*: with no Python chosen, reticulate sets up its own
environment with tskit (\>= 1.0) and numpy the first time one of these
functions runs; this needs internet once. *Your own Python*: set
`RETICULATE_PYTHON` to a Python that has those modules (in `~/.Renviron`
to keep it) before the first such call of a session. A function that
cannot reach them stops with an error of class `malariafwd_no_python`
that says what to do, and
[`check_setup()`](https://shaziaruybal.github.io/malariafwd/reference/check_setup.md)
reports what is found.

## See also

[`plot_ibd_pair()`](https://shaziaruybal.github.io/malariafwd/reference/plot_ibd_pair.md),
[`ibd_pairs()`](https://shaziaruybal.github.io/malariafwd/reference/ibd_pairs.md)

## Examples

``` r
run <- example_run()
rel <- load_events(run)[event_type == "hypnozoite_activation" & parent1_pedigree_id != -1][1]
ibd_segments(run, rel$pedigree_id, rel$parent1_pedigree_id)   # all 14 chromosomes shared
#>      pair     a     b chrom start     end   gstart     gend  length mrca_id
#>     <int> <int> <int> <int> <num>   <num>    <num>    <num>   <num>   <int>
#>  1:     1   121    52     1     0 1021664        0  1021664 1021664      52
#>  2:     1   121    52     2     0  956327  1021664  1977991  956327      52
#>  3:     1   121    52     3     0  896704  1977991  2874695  896704      52
#>  4:     1   121    52     4     0 1012024  2874695  3886719 1012024      52
#>  5:     1   121    52     5     0 1524814  3886719  5411533 1524814      52
#>  6:     1   121    52     6     0 1042791  5411533  6454324 1042791      52
#>  7:     1   121    52     7     0 1652210  6454324  8106534 1652210      52
#>  8:     1   121    52     8     0 1761288  8106534  9867822 1761288      52
#>  9:     1   121    52     9     0 2237066  9867822 12104888 2237066      52
#> 10:     1   121    52    10     0 1548844 12104888 13653732 1548844      52
#> 11:     1   121    52    11     0 2131221 13653732 15784953 2131221      52
#> 12:     1   121    52    12     0 3182763 15784953 18967716 3182763      52
#> 13:     1   121    52    13     0 2093556 18967716 21061272 2093556      52
#> 14:     1   121    52    14     0 3153402 21061272 24214674 3153402      52
#>     mrca_time depth
#>         <num> <num>
#>  1:       387     0
#>  2:       387     0
#>  3:       387     0
#>  4:       387     0
#>  5:       387     0
#>  6:       387     0
#>  7:       387     0
#>  8:       387     0
#>  9:       387     0
#> 10:       387     0
#> 11:       387     0
#> 12:       387     0
#> 13:       387     0
#> 14:       387     0
```
