# Validation results

The validation scripts in this directory, run on malariafwd 0.1.0 at commit 611a301 on 2026-11-10 01:11 UTC.
Every script passed.

| | version |
|---|---|
| R | 4.5.2 |
| SLiM | 5.2 |
| Python | 3.13 |
| tskit | 1.0.4 |
| numpy | 2.5.4 |
| platform | aarch64-apple-darwin20 |

| script | result | minutes |
|---|---|---|
| `audit_events_and_ibd.R` | 123 passed, 0 failed | 0.1 |
| `check_parameters.R` | 40 passed, 0 failed | 1.5 |
| `audit_arg_recombination.R` | ran without error (a numeric audit: run it to see its numbers) | 0.1 |
| `audit_arg_junctions.R` | ran without error (a numeric audit: run it to see its numbers) | 0.1 |

`dev/validation-results.R`, in the package's source repository, writes this file.
