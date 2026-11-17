# The IBD of pairs of genomes, from ibd_pairs()

`ibd_lookup()` gives the IBD of each pair (`a[i]`, `b[i]`): 1 for a
genome with itself, the pair's value otherwise, and `NA` when either
genome is not in the object or, for an object of listed pairs
(`ibd_pairs(pairs =)`), when the pair was not listed. A length-1
argument is recycled against the other.
[`as.matrix()`](https://rdrr.io/r/base/matrix.html) gives the pairs as a
symmetric matrix named by pedigree id, with a unit diagonal, and
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) the
`pairs` table.

## Usage

``` r
ibd_lookup(x, a, b)

# S3 method for class 'vivax_ibd'
as.matrix(x, ids = NULL, ...)

# S3 method for class 'vivax_ibd'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  a `vivax_ibd` object, from
  [`ibd_pairs()`](https://shaziaruybal.github.io/malariafwd/reference/ibd_pairs.md).

- a, b:

  pedigree ids.

- ids:

  for [`as.matrix()`](https://rdrr.io/r/base/matrix.html), the genomes
  to give, in order; `NULL` for every one.

- ...:

  ignored.

- row.names, optional:

  ignored.

## Value

`ibd_lookup()`: a numeric vector;
[`as.matrix()`](https://rdrr.io/r/base/matrix.html): a matrix;
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html): a
data.frame of `a`, `b` and `ibd`.

## See also

[`ibd_pairs()`](https://shaziaruybal.github.io/malariafwd/reference/ibd_pairs.md)

## Examples

``` r
run <- example_run()
ids <- classify_cases(run)$pedigree_id[1:30]
x <- ibd_pairs(run, ids = ids)
ibd_lookup(x, ids[1], ids[2:5])
#> [1] 0 0 1 0
```
