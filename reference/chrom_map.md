# Chromosome layout of a run's concatenated genome

The model simulates the 14 PvP01 chromosomes as ONE concatenated linkage
group. This recovers the chromosome boundaries in those concatenated
coordinates, from the genome parameters the run itself recorded
(`GENOME_MODE`, `CHROMOSOMES`), and stops if they do not add up to the
genome length the run recorded.

## Usage

``` r
chrom_map(run)
```

## Arguments

- run:

  a run: a `vivax_run` (from
  [`run_sim()`](https://shaziaruybal.github.io/malariafwd/reference/run_sim.md)
  or
  [`load_run()`](https://shaziaruybal.github.io/malariafwd/reference/load_run.md))
  or the path of a run directory.

## Value

a data.table of `chrom` (PvP01 number), `start`, `end`, `length` in
concatenated bp.

## Examples

``` r
chrom_map(example_run())
#>     chrom  length    start      end
#>     <int>   <num>    <num>    <num>
#>  1:     1 1021664        0  1021664
#>  2:     2  956327  1021664  1977991
#>  3:     3  896704  1977991  2874695
#>  4:     4 1012024  2874695  3886719
#>  5:     5 1524814  3886719  5411533
#>  6:     6 1042791  5411533  6454324
#>  7:     7 1652210  6454324  8106534
#>  8:     8 1761288  8106534  9867822
#>  9:     9 2237066  9867822 12104888
#> 10:    10 1548844 12104888 13653732
#> 11:    11 2131221 13653732 15784953
#> 12:    12 3182763 15784953 18967716
#> 13:    13 2093556 18967716 21061272
#> 14:    14 3153402 21061272 24214674
```
