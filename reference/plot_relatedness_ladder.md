# Plot the relatedness ladder

IBD by known true relationship, against the textbook rungs (0 unrelated,
1/4 half-sib, 1/3 the mean over a meiotic tetrad, 1/2 full-sib, 1
clonal). The model-validation view: clonal relationships should sit at
1, distinct founders at 0, and same-oocyst siblings (any two products of
one meiosis) in two groups: meiosis II sisters around 0.6 and the other
pairs around 0.2, a third on average for unrelated parents and more when
the parents are related. On the full map none shares nothing.

## Usage

``` r
plot_relatedness_ladder(x, style = c("points", "bars"))
```

## Arguments

- x:

  a run, or the table
  [`ibd_truth_pairs()`](https://shaziaruybal.github.io/malariafwd/reference/ibd_truth_pairs.md)
  returns.

- style:

  `"points"` plots every pair against the theoretical rungs; `"bars"`
  shows the composition of each class across relatedness bands, which
  reads better when classes are point masses (relapses are always
  clonal, distinct founders always unrelated). Both label the number of
  pairs per class, so a class resting on two observations is visible.

## Value

a ggplot object.

## Examples

``` r
plot_relatedness_ladder(ibd_truth_pairs(example_run()))
```
