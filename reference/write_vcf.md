# Write simulated genotypes as a VCF file

Writes a
[`simulate_genotypes()`](https://shaziaruybal.github.io/malariafwd/reference/simulate_genotypes.md)
object as VCF 4.2: one contig per simulated PvP01 chromosome, with the
run's lengths; each locus at its 1-based position within its chromosome
(so in `"core"` mode, a position on the shortened chromosome), its id in
`ID`, and the ALT frequency it was painted with as `PLAF` in `INFO`,
rather than `AF`, which tools read as the frequency among the samples.
Records are sorted by chromosome and position. Clone samples are written
haploid (`0`, `1`, `.`) and host samples unphased diploid, as VCFs
called from mixed *Plasmodium* infections are: `0/0`, `1/1`, `0/1` for
both alleles present, `./.`. The model has no sequence, so without `ref`
and `alt` columns in the loci REF and ALT are the placeholders A and C,
and the header says so.

## Usage

``` r
write_vcf(x, file, calls = c("observed", "true"), chrom_names = NULL)
```

## Arguments

- x:

  a `vivax_genotypes` object, from
  [`simulate_genotypes()`](https://shaziaruybal.github.io/malariafwd/reference/simulate_genotypes.md).

- file:

  the path to write; a name ending in `.gz` is gzip-compressed (for
  tabix, write a plain file and compress it with bgzip).

- calls:

  `"observed"` (the calls with error and missingness) or `"true"`.

- chrom_names:

  `NULL`, or the 14 contig names by PvP01 number; by default
  `PvP01_01_v2` to `PvP01_14_v2`, and `as.character(1:14)` gives
  numbered contigs.

## Value

`file`, invisibly.

## See also

[`simulate_genotypes()`](https://shaziaruybal.github.io/malariafwd/reference/simulate_genotypes.md)

## Examples

``` r
run <- example_run()
loci <- data.frame(chrom = 1:14, pos = 1e5, plaf = 0.3)
g <- simulate_genotypes(run, classify_cases(run)$pedigree_id[1:10], loci, seed = 1)
f <- write_vcf(g, tempfile(fileext = ".vcf"))
readLines(f, n = 12)
#>  [1] "##fileformat=VCFv4.2"                                                        
#>  [2] "##source=malariafwd 0.1.0 simulate_genotypes()"                              
#>  [3] "##malariafwd_run=/tmp/RtmpUQ7ECQ/malariafwd-example-run-60b2e786"            
#>  [4] "##malariafwd_samples=clone"                                                  
#>  [5] "##malariafwd_calls=observed (error=0, missing=0)"                            
#>  [6] "##malariafwd_seed=1"                                                         
#>  [7] "##malariafwd_alleles=REF and ALT are placeholders: the model has no sequence"
#>  [8] "##contig=<ID=PvP01_01_v2,length=1021664>"                                    
#>  [9] "##contig=<ID=PvP01_02_v2,length=956327>"                                     
#> [10] "##contig=<ID=PvP01_03_v2,length=896704>"                                     
#> [11] "##contig=<ID=PvP01_04_v2,length=1012024>"                                    
#> [12] "##contig=<ID=PvP01_05_v2,length=1524814>"                                    
```
