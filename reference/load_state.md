# Load the per-tick state table of a run

Reads `state.tsv`, by default its host rows. A mosquito's row counts the
lineages it carries, so a prevalence or MOI taken over every row mixes
the two populations. Every run gives the same 12 columns, in the order
and classes the model writes them
([state.tsv](https://shaziaruybal.github.io/malariafwd/reference/malariafwd-outputs.md)
describes them), whatever the run holds; a column the file lacks is
`NA`.

## Usage

``` r
load_state(run, entity = c("human", "mosquito", "all"))
```

## Arguments

- run:

  a run: a `vivax_run` (from
  [`run_sim()`](https://shaziaruybal.github.io/malariafwd/reference/run_sim.md)
  or
  [`load_run()`](https://shaziaruybal.github.io/malariafwd/reference/load_run.md))
  or the path of a run directory.

- entity:

  the rows to give: `"human"` (the hosts), `"mosquito"`, or `"all"`,
  both in the order of the file.

## Value

a `data.table` with one row per host (or mosquito, or both) per day,
from day 1 to the run's last day.

## Examples

``` r
st <- load_state(example_run())
st[tick == 400]   # every host on the last day
#>      tick entity_type entity_id   tag moi_active moi_subpatent moi_dormant
#>     <int>      <char>     <int> <int>      <int>         <int>       <int>
#>  1:   400       human         0     0          0             0           0
#>  2:   400       human         1     0          0             0           1
#>  3:   400       human         2     0          0             0           3
#>  4:   400       human         3     0          0             0           0
#>  5:   400       human         4     1          1             0           2
#>  6:   400       human         5     1          1             0           1
#>  7:   400       human         6     1          3             0           3
#>  8:   400       human         7     0          0             0           0
#>  9:   400       human         8     1          2             0           0
#> 10:   400       human         9     0          0             0           1
#> 11:   400       human        10     1          1             0           1
#> 12:   400       human        11     0          0             0           2
#> 13:   400       human        12     0          0             0           2
#> 14:   400       human        13     0          0             0           1
#> 15:   400       human        14     0          0             0           0
#> 16:   400       human        15     1          1             0           3
#> 17:   400       human        16     1          2             0           1
#> 18:   400       human        17     0          0             0           0
#> 19:   400       human        18     0          0             0           0
#> 20:   400       human        19     0          0             0           1
#> 21:   400       human        20     0          0             0           0
#> 22:   400       human        21     0          0             0           1
#> 23:   400       human        22     0          0             0           0
#> 24:   400       human        23     0          0             0           0
#> 25:   400       human        24     0          0             0           0
#> 26:   400       human        25     0          0             0           0
#> 27:   400       human        26     0          0             0           1
#> 28:   400       human        27     0          0             0           1
#> 29:   400       human        28     0          0             0           0
#> 30:   400       human        29     0          0             0           0
#> 31:   400       human        30     1          1             0           1
#> 32:   400       human        31     1          1             0           1
#> 33:   400       human        32     0          0             0           0
#> 34:   400       human        33     1          1             0           1
#> 35:   400       human        34     1          2             0           1
#> 36:   400       human        35     0          0             0           1
#> 37:   400       human        36     0          0             0           2
#> 38:   400       human        37     1          1             0           1
#> 39:   400       human        38     0          0             0           0
#> 40:   400       human        39     0          0             0           0
#> 41:   400       human        40     0          0             0           4
#> 42:   400       human        41     0          0             0           1
#> 43:   400       human        42     0          0             0           0
#> 44:   400       human        43     1          1             0           2
#> 45:   400       human        44     0          0             0           0
#> 46:   400       human        45     0          0             0           0
#> 47:   400       human        46     1          1             0           3
#> 48:   400       human        47     0          0             0           0
#> 49:   400       human        48     0          0             0           0
#> 50:   400       human        49     0          0             0           1
#> 51:   400       human        50     0          0             0           2
#> 52:   400       human        51     0          0             0           2
#> 53:   400       human        52     1          2             0           2
#> 54:   400       human        53     0          0             0           0
#> 55:   400       human        54     0          0             0           1
#> 56:   400       human        55     0          0             0           0
#> 57:   400       human        56     0          0             0           0
#> 58:   400       human        57     0          0             0           2
#> 59:   400       human        58     1          4             0           4
#> 60:   400       human        59     1          1             0           1
#>      tick entity_type entity_id   tag moi_active moi_subpatent moi_dormant
#>     <int>      <char>     <int> <int>      <int>         <int>       <int>
#>     hyp_reservoir_size            lineage_ids n_transmissible mosq_age_days
#>                  <int>                 <char>           <int>         <int>
#>  1:                  0                   <NA>              NA            NA
#>  2:                  3                      7              NA            NA
#>  3:                  8               0|1034|7              NA            NA
#>  4:                  0                   <NA>              NA            NA
#>  5:                  3                275|697              NA            NA
#>  6:                  2                  7|923              NA            NA
#>  7:                 13               0|7|1005              NA            NA
#>  8:                  0                   <NA>              NA            NA
#>  9:                  0                225|697              NA            NA
#> 10:                  5                   1127              NA            NA
#> 11:                  1                      7              NA            NA
#> 12:                  2                  5|225              NA            NA
#> 13:                  4                    2|7              NA            NA
#> 14:                  1                    648              NA            NA
#> 15:                  0                   <NA>              NA            NA
#> 16:                  7           1240|225|275              NA            NA
#> 17:                  2                 0|1028              NA            NA
#> 18:                  0                   <NA>              NA            NA
#> 19:                  0                   <NA>              NA            NA
#> 20:                  1                    649              NA            NA
#> 21:                  0                   <NA>              NA            NA
#> 22:                  3                    225              NA            NA
#> 23:                  0                   <NA>              NA            NA
#> 24:                  0                   <NA>              NA            NA
#> 25:                  0                   <NA>              NA            NA
#> 26:                  0                   <NA>              NA            NA
#> 27:                  1                    800              NA            NA
#> 28:                  8                    649              NA            NA
#> 29:                  0                   <NA>              NA            NA
#> 30:                  0                   <NA>              NA            NA
#> 31:                  6                    697              NA            NA
#> 32:                  2                275|697              NA            NA
#> 33:                  0                   <NA>              NA            NA
#> 34:                  2                      2              NA            NA
#> 35:                  2               1017|411              NA            NA
#> 36:                  1                    415              NA            NA
#> 37:                  8               697|1041              NA            NA
#> 38:                  1                    697              NA            NA
#> 39:                  0                   <NA>              NA            NA
#> 40:                  0                   <NA>              NA            NA
#> 41:                  9            0|746|7|946              NA            NA
#> 42:                  3                    275              NA            NA
#> 43:                  0                   <NA>              NA            NA
#> 44:                  8                  973|0              NA            NA
#> 45:                  0                   <NA>              NA            NA
#> 46:                  0                   <NA>              NA            NA
#> 47:                  9            680|993|649              NA            NA
#> 48:                  0                   <NA>              NA            NA
#> 49:                  0                   <NA>              NA            NA
#> 50:                  1                   1566              NA            NA
#> 51:                  2                981|697              NA            NA
#> 52:                  3                649|648              NA            NA
#> 53:                  3              1114|1115              NA            NA
#> 54:                  0                   <NA>              NA            NA
#> 55:                  2                   1345              NA            NA
#> 56:                  0                   <NA>              NA            NA
#> 57:                  0                   <NA>              NA            NA
#> 58:                  3                 7|1005              NA            NA
#> 59:                  8 649|650|1281|1285|1279              NA            NA
#> 60:                  1                    5|2              NA            NA
#>     hyp_reservoir_size            lineage_ids n_transmissible mosq_age_days
#>                  <int>                 <char>           <int>         <int>
#>     infectious_bites
#>                <int>
#>  1:                0
#>  2:                0
#>  3:                0
#>  4:                0
#>  5:                0
#>  6:                0
#>  7:                0
#>  8:                0
#>  9:                0
#> 10:                0
#> 11:                0
#> 12:                0
#> 13:                0
#> 14:                0
#> 15:                0
#> 16:                0
#> 17:                0
#> 18:                0
#> 19:                0
#> 20:                0
#> 21:                0
#> 22:                0
#> 23:                0
#> 24:                0
#> 25:                0
#> 26:                0
#> 27:                0
#> 28:                0
#> 29:                0
#> 30:                0
#> 31:                0
#> 32:                0
#> 33:                0
#> 34:                0
#> 35:                0
#> 36:                0
#> 37:                0
#> 38:                0
#> 39:                0
#> 40:                0
#> 41:                0
#> 42:                0
#> 43:                0
#> 44:                0
#> 45:                0
#> 46:                0
#> 47:                0
#> 48:                0
#> 49:                0
#> 50:                0
#> 51:                0
#> 52:                0
#> 53:                0
#> 54:                0
#> 55:                0
#> 56:                0
#> 57:                0
#> 58:                0
#> 59:                0
#> 60:                0
#>     infectious_bites
#>                <int>
```
