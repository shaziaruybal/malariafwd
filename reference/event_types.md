# The malariafwd event dictionary

Every row in a run's `events.tsv` carries an `event_type`. This returns
a tidy table of the event types the model emits and what each means, so
you can navigate an event log. An inoculation that establishes neither a
blood-stage infection nor a hypnozoite batch is logged with a `_no_take`
suffix (e.g. `vh_inoculation_no_take`), and so is an import whose host
already carries `MAX_MOI` blood-stage lineages
(`vh_inoculation_import_no_take`); a local bite refused at `MAX_MOI`
leaves no row.

## Usage

``` r
event_types()
```

## Value

a data.table with columns `event_type`, `stage`, `ends_episode` (does
this event end a patent blood-stage episode?), `is_onset` (does it begin
one?) and `description`.

## Details

`ends_episode` marks every way a patent blood-stage episode can end and
`is_onset` every way one can begin: filter on them rather than listing
event names by hand.

`mosq_death` is one row per parasite genome an infected mosquito carried
when it died; the deaths of uninfected mosquitoes are not logged, so its
rows do not count deaths.

## See also

[events.tsv](https://shaziaruybal.github.io/malariafwd/reference/malariafwd-outputs.md)

## Examples

``` r
event_types()
#>                        event_type  stage ends_episode is_onset
#>                            <char> <char>       <lgcl>   <lgcl>
#>  1:                vh_inoculation  blood        FALSE     TRUE
#>  2:         vh_inoculation_import  blood        FALSE     TRUE
#>  3:           vh_inoculation_seed  blood        FALSE     TRUE
#>  4:        vh_inoculation_no_take  blood        FALSE    FALSE
#>  5: vh_inoculation_import_no_take  blood        FALSE    FALSE
#>  6:         hypnozoite_activation  blood        FALSE     TRUE
#>  7:                 recrudescence  blood        FALSE     TRUE
#>  8:               human_clearance  blood         TRUE    FALSE
#>  9:         case_bloodstage_clear  blood         TRUE    FALSE
#> 10:          mda_bloodstage_clear  blood         TRUE    FALSE
#> 11:           mda_subpatent_clear  blood        FALSE    FALSE
#> 12:           bloodcure_subpatent  blood         TRUE    FALSE
#> 13:               bloodcure_given  blood        FALSE    FALSE
#> 14:                 radcure_given  liver        FALSE    FALSE
#> 15:               radcure_failure  liver        FALSE    FALSE
#> 16:             bloodcure_persist  blood        FALSE    FALSE
#> 17:           mda_bloodcure_given  blood        FALSE    FALSE
#> 18:             mda_radcure_given  liver        FALSE    FALSE
#> 19:           mda_radcure_failure  liver        FALSE    FALSE
#> 20:          hypnozoite_clearance  liver        FALSE    FALSE
#> 21:      radcure_hypnozoite_clear  liver        FALSE    FALSE
#> 22:          mda_hypnozoite_clear  liver        FALSE    FALSE
#> 23:              hv_uptake_clonal vector        FALSE    FALSE
#> 24:         hv_uptake_recombinant vector        FALSE    FALSE
#> 25:                    mosq_death vector        FALSE    FALSE
#>                        event_type  stage ends_episode is_onset
#>                            <char> <char>       <lgcl>   <lgcl>
#>                                                                                                                                                                                                    description
#>                                                                                                                                                                                                         <char>
#>  1:                                                                                                                                                             new local mosquito-borne blood-stage infection
#>  2:                                                                                                                                                                             imported blood-stage infection
#>  3:                                                                                                                                                            founder blood-stage infection (burn-in seeding)
#>  4:                                                                                                 a bite established neither a blood-stage infection nor hypnozoites (P_PRIMARY failed and the batch drew 0)
#>  5:                                                                                         an importation established neither compartment (P_PRIMARY failed and the batch drew 0), or its host was at MAX_MOI
#>  6:                                                                                                                                       hypnozoite activates -> a blood-stage relapse (clonal to its origin)
#>  7:                                                                              a suppressed clone re-expands to patency after a failed course; mints a new genome clonal to the holder (P_RECRUD_GIVEN_FAIL)
#>  8:                                                                                                                                                                     blood-stage infection clears naturally
#>  9:                                                                                                                                 blood-stage infection cleared by case-based treatment (BLOODCURE_COVERAGE)
#> 10:                                              blood-stage infection cleared by mass drug administration (MDA_TICK), also one that became patent earlier that day (the day's imports arrive after the round)
#> 11:                                                                                    a suppressed clone cleared by MDA, so it never recrudesces; its patent episode had already ended at bloodcure_subpatent
#> 12: a clone survived treatment but was suppressed BELOW PATENCY; the patent episode ends here and the clone recrudesces later (P_RECRUD_GIVEN_FAIL), unless an MDA round clears it first (mda_subpatent_clear)
#> 13:                                                                                                                                               a blood-stage course is given to a host (BLOODCURE_COVERAGE)
#> 14:                                                                                                                                        radical cure is given to a host, nested inside a blood-stage course
#> 15:                                                                                                                            radical cure was given but the whole reservoir survived (RADCURE_MODE = person)
#> 16:                                                                                                                a blood-stage clone survived treatment and continues at full density (P_BLOODCURE_SURVIVAL)
#> 17:                                                                                                                       a blood-stage course is given by MDA to a host not empty at the round (MDA_COVERAGE)
#> 18:                                                                                                                                                                               radical cure is given by MDA
#> 19:                                                                                                                                                    MDA radical cure given but the whole reservoir survived
#> 20:                                                                                                                                                                    a hypnozoite dies / clears in the liver
#> 21:                                                                                                                                 hypnozoite reservoir thinned by case-based radical cure (RADCURE_COVERAGE)
#> 22:                                                                                                                                                                        hypnozoite reservoir thinned by MDA
#> 23:                                                                                                                                                        mosquito takes up a clonal (selfed) oocyst genotype
#> 24:                                                                                                                                               mosquito takes up a recombinant (outcrossed) oocyst genotype
#> 25:                                                                                                                       a parasite genome lost when an infected mosquito dies: one row per genome it carried
#>                                                                                                                                                                                                    description
#>                                                                                                                                                                                                         <char>
subset(event_types(), is_onset)$event_type       # every way an episode can begin
#> [1] "vh_inoculation"        "vh_inoculation_import" "vh_inoculation_seed"  
#> [4] "hypnozoite_activation" "recrudescence"        
```
