# Validation scripts

Checks that go beyond the unit tests in `tests/`: the unit tests protect the *code*, these
protect the *science*. Run them after any change to the SLiM model, and before a release.

| script | what it checks |
|---|---|
| `audit_events_and_ibd.R` | event-log integrity, hypnozoite & blood-stage conservation against `state.tsv`, tree retention, transmission links, IBD invariants, ARG/TMRCA, genotype painting |
| `check_parameters.R` | 34 configurations: structural invariants hold everywhere, every parameter moves the model the right way |
| `calibrate_nu.R` | fits `NU_HYPNOZOITE_BATCH` against four calibration targets at moderate transmission: the table in section 13 of the model description, and that the shipped value sits inside every target (15 runs, about 6 minutes on one core, less with `MALARIAFWD_WORKERS` set) |
| `audit_arg_recombination.R` | breakpoints follow the recombination map (uniform within chromosomes); a product of meiosis changes parent at rate r, so the segments it shares with a parent are exponential with mean 1/r, cut short at the chromosome ends: (14 + rL)/2 of them, of mean L/(14 + rL) |
| `audit_arg_junctions.R` | chromosomes assort independently across the 13 junctions; how much the products of one meiosis share: about a third for unrelated parents, more for related ones, and on the full map never nothing; the segments two products share, (14 + 2rL)/3 a pair of mean L/(14 + 2rL) |

They need SLiM on `PATH` and Python with `tskit` (>= 1.0) and `numpy`;
`malariafwd::check_setup()` reports both. Run from anywhere:

```sh
Rscript -e 'source(system.file("validation", "check_parameters.R", package = "malariafwd"))'
```

Figures are written to `$MALARIAFWD_FIGS` (default: a temporary directory). `RESULTS.md` records
the last full run: its date, commit and versions, and what each script found.
