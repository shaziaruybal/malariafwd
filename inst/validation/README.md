# Validation scripts

Checks that go beyond the unit tests in `tests/`: the unit tests protect the *code*, these
protect the *science*. Run them after any change to the SLiM model, and before a release.

| script | what it checks |
|---|---|
| `audit_events_and_ibd.R` | event-log integrity, hypnozoite & blood-stage conservation against `state.tsv`, tree retention, IBD invariants, ARG/TMRCA (62 checks) |
| `check_parameters.R` | 34 configurations: structural invariants hold everywhere, every parameter moves the model the right way, seed reproducibility (35 checks) |
| `audit_arg_recombination.R` | breakpoints follow the recombination map (uniform within chromosomes) |
| `audit_arg_junctions.R` | the 13 chromosome junctions assort independently; why same-oocyst sibs can be IBD = 0 |
| `explain_ibd.R` | renders the IBD segment and relatedness-ladder figures |

They need SLiM on `PATH` and Python with `tskit`. Run from anywhere:

```sh
Rscript -e 'source(system.file("validation", "check_parameters.R", package = "malariafwd"))'
```

Figures are written to `$MALARIAFWD_FIGS` (default: a temporary directory).
