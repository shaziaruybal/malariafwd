# Simulate marker genotypes from a run's genealogy

Biallelic genotypes at chosen loci for chosen genomes, painted on the
genealogy the run recorded. The model carries no mutation, so a run's
genomes differ only by descent: at every position each genome descends
from exactly one **root genome**, a day-1 founder template or an
imported genome, and root genomes are unrelated. Each root genome draws
its allele at each locus independently, ALT with probability `plaf`, and
every genome carries the allele of the root genome it descends from
there. `roots` records that root genome for every genome and locus, the
truth behind the calls: two genomes share a root genome at a locus
exactly where
[`ibd_segments()`](https://shaziaruybal.github.io/malariafwd/reference/ibd_segments.md)
finds them sharing a common ancestor.

## Usage

``` r
simulate_genotypes(
  run,
  ids,
  loci,
  sample_id = NULL,
  error = 0,
  missing = 0,
  seed = NULL
)

# S3 method for class 'vivax_genotypes'
print(x, ...)
```

## Arguments

- run:

  a run: a `vivax_run` (from
  [`run_sim()`](https://shaziaruybal.github.io/malariafwd/reference/run_sim.md)
  or
  [`load_run()`](https://shaziaruybal.github.io/malariafwd/reference/load_run.md))
  or the path of a run directory.

- ids:

  the pedigree ids of the genomes to genotype.

- loci:

  a data frame, one row per biallelic locus: `chrom` (the PvP01
  chromosome, as a number or a name such as `PvP01_03_v2`), `pos`
  (1-based, within the chromosome as this run simulated it) and `plaf`
  (the frequency of the ALT allele among the root genomes); optionally
  `locus` (an id, by default `"<chrom>:<pos>"`), and `ref` and `alt` for
  [`write_vcf()`](https://shaziaruybal.github.io/malariafwd/reference/write_vcf.md).
  Other columns are kept.

- sample_id:

  `NULL` for clone samples, one per genome, or a vector as long as `ids`
  grouping the genomes into host samples; a genome can be in several
  samples.

- error:

  the probability that a call is wrong: a clone call flips, and a host
  call becomes one of the other two at random.

- missing:

  the probability that a call is missing, after the errors.

- seed:

  `NULL`, or a whole number that seeds this call alone.

- x:

  a `vivax_genotypes` object.

- ...:

  ignored.

## Value

a `vivax_genotypes` object, a list of:

- `calls`: a sample x locus matrix of 0, 1, 0.5 (host samples only) and
  `NA`;

- `true_calls`: the calls before error and missingness;

- `alleles`: an integer genome x locus matrix of 0 (REF) and 1 (ALT), by
  pedigree id;

- `roots`: an integer genome x locus matrix, the pedigree id of the root
  genome behind each allele;

- `samples` (`sample_id`, `n_clones`) and `clones` (`sample_id`,
  `pedigree_id`, in input order);

- `loci`: `locus`, `chrom`, `pos`, `position` (0-based bp along the
  concatenated genome, as in the IBD segments), `plaf`, then the other
  columns given;

- `root_genomes` (`root_id`, `tick`: day 1 for a founder template, the
  arrival day for an import, and `imported`) and `root_alleles`, their
  alleles;

- `chromosomes`, the run's
  [`chrom_map()`](https://shaziaruybal.github.io/malariafwd/reference/chrom_map.md);
  `type` (`"clone"` or `"host"`), `error`, `missing`, `seed` (`NA` when
  not set) and `run`.

## Details

Each genome can be its own sample (a clone), or the genomes can be
grouped into host samples by `sample_id`, as a blood sample holds every
clone in the host: a host call is 0 when every clone carries REF, 1 when
every one carries ALT, and 0.5 when both are present.

## Reproducibility

The root genomes' alleles are drawn first, for every root genome of the
run and every locus, so with one seed, run and loci a genome's alleles
do not depend on which other genomes are asked for, how they are
grouped, or the error and missingness rates; then the errors, then the
missing calls. Large sets of genomes can therefore be genotyped in
batches with one seed, and several error rates laid over one truth. A
number in `seed` seeds this call alone and leaves the session's random
stream as it was; `NULL` draws from the session's stream.

## What painting leaves out

Root genomes are unrelated, so there is no background relatedness or
linkage disequilibrium among them; imports draw from the same
frequencies as the founders; an allele never changes along the
genealogy; a host sample records which alleles are present, not their
proportions; loci are biallelic; and `plaf` is the frequency among the
root genomes, which drifts among the genomes sampled later.

## Which genomes can be genotyped

Any genome in the tree sequence: case clones, hypnozoite batches and
mosquito products alike. A run with `RETAIN_FULL_HISTORY = FALSE` keeps
few genomes beyond those alive at its end; a run that never started
holds none.

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

[`write_vcf()`](https://shaziaruybal.github.io/malariafwd/reference/write_vcf.md),
[`ibd_segments()`](https://shaziaruybal.github.io/malariafwd/reference/ibd_segments.md),
[`ibd_pairs()`](https://shaziaruybal.github.io/malariafwd/reference/ibd_pairs.md),
[`chrom_map()`](https://shaziaruybal.github.io/malariafwd/reference/chrom_map.md),
[`episodes()`](https://shaziaruybal.github.io/malariafwd/reference/episodes.md)

## Examples

``` r
run <- example_run()
cm <- chrom_map(run)
loci <- data.frame(chrom = rep(cm$chrom, each = 5), pos = 1000 * rep(1:5, nrow(cm)), plaf = 0.3)
g <- simulate_genotypes(run, classify_cases(run)$pedigree_id[1:20], loci, error = 0.01,
                        missing = 0.05, seed = 1)
g
#> <vivax_genotypes> 20 clone samples x 70 loci on 14 chromosomes
#>   painted on 33 root genomes (23 imported); error 0.01, missing 0.05, seed 1
#>   run: /tmp/RtmpjMBTAb/malariafwd-example-run-60b2e786
g$calls[1:4, 1:5]   # what a method sees
#>     1:1000 1:2000 1:3000 1:4000 1:5000
#> 51       0      1      0      0      0
#> 82       0      0      1      0      0
#> 98       0      0      0      1      0
#> 121      0      1      0      0     NA
g$roots[1:4, 1:5]   # the truth behind it
#>     1:1000 1:2000 1:3000 1:4000 1:5000
#> 51       7      7      7      7      7
#> 82      82     82     82     82     82
#> 98       4      4      4      4      4
#> 121      7      7      7      7      7
```
