# Pairwise time to the most recent common ancestor, the relatedness measure of the ancestral recombination graph

For each pair, the span-weighted mean TMRCA **over the part of the
genome where a common ancestor exists**, plus `coalesced`, the fraction
of the genome where it does. This is the time-domain complement to IBD:
IBD asks *how much* of the genome is shared, `arg_tmrca()` asks *how far
back* the sharing goes. Unlike tskit's branch-mode divergence it is
meaningful on a forest, because non-coalescing spans are excluded rather
than scored as 0.

## Usage

``` r
arg_tmrca(run, ids_a, ids_b)
```

## Arguments

- run:

  a run: a `vivax_run` (from
  [`run_sim()`](https://shaziaruybal.github.io/malariafwd/reference/run_sim.md)
  or
  [`load_run()`](https://shaziaruybal.github.io/malariafwd/reference/load_run.md))
  or the path of a run directory.

- ids_a, ids_b:

  equal-length vectors of pedigree ids, compared pairwise.

## Value

a data.table, one row per pair: `a`, `b`, `coalesced` (0-1), `tmrca` and
`tmrca_node_time` (both ticks; `NA` if the pair never coalesces or a
genome is not in the tree sequence). A run that recorded no genome gives
an `NA` row per pair, and none of the ids being a genome of a run that
recorded some is an error.

## Details

Clonal pairs give `coalesced` = 1 with a very recent TMRCA; unrelated
genomes give `coalesced` = 0 and `tmrca` = `NA`.

## Which TMRCA you get

`tmrca` is a **depth**: ticks from the older of the two sampled genomes
back to their common ancestor, so an identical relationship returns an
identical number wherever in the run it occurred – a clonal parent-child
pair gives 0 whether it happened on tick 10 or tick 390.
`tmrca_node_time` is the raw tskit node time of that ancestor, i.e.
ticks before the **end of the run**; that is what `tree.tmrca()` returns
natively and it is dominated by when the pair was sampled rather than by
how related it is. Use `tmrca` for relatedness.

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

## Examples

``` r
run <- example_run()
rel <- load_events(run)[event_type == "hypnozoite_activation" & parent1_pedigree_id != -1][1:3]
arg_tmrca(run, rel$pedigree_id, rel$parent1_pedigree_id)   # relapses vs their batches
#>        a     b coalesced tmrca tmrca_node_time
#>    <int> <int>     <num> <num>           <num>
#> 1:   121    52         1     0             387
#> 2:   134    52         1     0             387
#> 3:   152    52         1     0             387
```
