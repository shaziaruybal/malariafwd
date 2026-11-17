# Ground-truth IBD from a run's tree sequence

Computes genome-wide identity-by-descent (shared-segment length /
sequence length) among a set of genomes, keyed by pedigree id, from the
tskit tree sequence recorded by the model.

## Usage

``` r
ibd_pairs(run, ids = NULL, pairs = NULL, min_span = 1000, segments = FALSE)
```

## Arguments

- run:

  a run: a `vivax_run` (from
  [`run_sim()`](https://shaziaruybal.github.io/malariafwd/reference/run_sim.md)
  or
  [`load_run()`](https://shaziaruybal.github.io/malariafwd/reference/load_run.md))
  or the path of a run directory.

- ids:

  the pedigree ids of the genomes to compare, every pair among them;
  `NULL` for every genome the tree sequence records, which warns when
  that is more than 1e6 pairs.

- pairs:

  instead of `ids`, a data frame of pairs to compare, with the pedigree
  ids in columns `a` and `b`: only these pairs are computed, about 0.2
  ms each, so use it for a few pairs among many genomes, such as each
  recurrence against its host's earlier episodes.

- min_span:

  the shortest segment returned, in bp, after the cut at the junctions:
  a segment of exactly `min_span` is kept. It does not change `ibd`.

- segments:

  if `TRUE`, also return every shared segment with its common ancestor
  and that ancestor's age.

## Value

a `vivax_ibd` object, a list of:

- `pairs`, a data.table keyed by `a` and `b`: one row for every pair of
  distinct genomes found, `a < b` (integer pedigree ids), and `ibd`
  (double, in \[0, 1\]), 0 for a pair that shares nothing;

- `samples`, one row per genome found, by pedigree id: `pedigree_id`,
  `node_time` (ticks before the end of the run) and `tick` (the day the
  genome arose, as in `events.tsv`);

- `missing`, the requested ids the tree sequence does not hold;

- `segments`, `NULL` or a data.table, by pair and position: `a`, `b`,
  `chrom` (the PvP01 chromosome), `start` and `end` (bp within the
  chromosome, `end` exclusive), `gstart` and `gend` (bp along the
  concatenated genome), `length` (bp), `mrca_id` (the pedigree id of the
  common ancestor), `mrca_time` and `depth` (its age, see *Segment
  ages*);

- `sequence_length`, `min_span`, `all_pairs` (`TRUE`: every pair of the
  genomes is in `pairs`; `FALSE` for listed `pairs`, which hold each
  listed pair of two distinct genomes found, once, with `a < b`),
  `end_tick` (the run's last tick) and `trees` (the file read).

Read a pair's value with
[`ibd_lookup()`](https://shaziaruybal.github.io/malariafwd/reference/ibd_lookup.md),
or the table with [`as.matrix()`](https://rdrr.io/r/base/matrix.html) or
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html).

## What is counted

A pair's `ibd` counts **all** shared ancestry: every stretch over which
the two genomes descend from a common ancestor genome, of any age and
any length, over the sequence length. `min_span` does not change it. To
impose a cutoff, ask for `segments` and filter them on their length or
their age (see *Segment ages*).

## Segments

A segment is a maximal stretch over which two genomes descend from one
common ancestor genome, cut at each chromosome junction: a change of
ancestor starts a new segment even with no gap between them. A pair gets
the same segments whatever other genomes are passed, the same as
[`ibd_segments()`](https://shaziaruybal.github.io/malariafwd/reference/ibd_segments.md)
gives, and at `min_span = 0` they add up to its `ibd`.

## Ids that are not in the tree sequence

In a run written with `RETAIN_FULL_HISTORY = FALSE`, most genomes not
alive at its end leave the tree sequence. A requested id it does not
hold is listed in `missing`, gives one warning, and looks up as `NA`,
not as 0; asking only for such ids is an error.

## Segment ages

`mrca_time` is a tskit node time, ticks before the end of the run, not
an age. A segment's `depth` is how long before the older of the two
genomes its common ancestor lived,
`mrca_time - max(node_time[a], node_time[b])`, in ticks. A recency
filter is written against it, for example `x$segments[depth <= 180]` for
ancestry within about six months: filtering on `mrca_time` alone scores
a pair by when its ancestry last branched, not by how related the two
genomes are.

## Cost

`segments = TRUE` takes a few times as long and gives one row per
segment per pair, so restrict it with `ids` and `min_span`.

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

[`ibd_lookup()`](https://shaziaruybal.github.io/malariafwd/reference/ibd_lookup.md),
[`ibd_segments()`](https://shaziaruybal.github.io/malariafwd/reference/ibd_segments.md)

## Examples

``` r
run <- example_run()
x <- ibd_pairs(run, ids = classify_cases(run)$pedigree_id[1:30])
x
#> <vivax_ibd> 30 genomes, 435 pairs, min_span 1000 bp
#>   segments: not computed
#>   from /tmp/RtmpsLWJEz/malariafwd-example-run-60b2e786/plasmo.trees, ending on tick 400
as.matrix(x)[1:4, 1:4]
#>     51 82 98 121
#> 51   1  0  0   1
#> 82   0  1  0   0
#> 98   0  0  1   0
#> 121  1  0  0   1
```
