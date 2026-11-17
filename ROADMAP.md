# malariafwd — ongoing and future work

Where sensitivity analysis is still needed, and the extensions planned for later versions. The
[model description](https://shaziaruybal.github.io/malariafwd/MODEL_DESCRIPTION.html) says what the
model does today, and [design decisions](https://shaziaruybal.github.io/malariafwd/DECISIONS.html)
says why.

## Where sensitivity analysis is needed

Ordered by how much the answer is likely to move.

| # | target | why | status |
|---|---|---|---|
| S1 | **`LAMBDA_SPOROZOITES`** (genotypes per bite) | no measurement of genotypes per *natural* bite exists; *Saimiri* monkeys given 13,000 sporozoites intravenously, from the glands of 20 mosquitoes, carried 2 and 4 genotypes (Hazzard et al. 2024), a far larger inoculum than a bite's | **priority** — sweep, don't trust a point value |
| S2 | **`NU_HYPNOZOITE_BATCH`** | fitted at 3 for one setting; drives relapse rate, reservoir size and COI | sweep per setting |
| S3 | **`LAMBDA_OOCYSTS` / `K_OOCYSTS`** | vivax membrane-feeding counts are overdispersed (Miura et al. 2020), with a median of 7.2 oocysts a midgut in wild *An. arabiensis* (Ayele et al. 2023), against the default's mean of 2.31 an uptake | sweep; NB option exists |
| S4 | **Per-oocyst sporozoite productivity** | we assume equal (weight 4); Childs & Prosper (2017) draw ~3,000 sporozoites an oocyst, Poisson, 20% of which reach the glands, so an early or productive oocyst would dominate | not implemented; test as sensitivity |
| S5 | **`GAMETOCYTE_DELAY` / `P_GAMETOCYTE`** | shipped at identity; vivax gametocytes appear early (Bousema & Drakeley 2011; Koepfli et al. 2021), so sweep a short delay and a high carriage (3 d and 0.9, say) | knobs exist for exactly this |
| S6 | **`P_DETECT_RELAPSE`** | relapses are often sub-patent; everything observational depends on it | sweep 0.15–1.0 |
| S7 | **`ALPHA_ACTIVATE` / `K_LATENCY`** | tropical vs temperate strongly changes relapse periodicity and identifiability | sweep; report tropical and temperate settings separately |
| S8 | **`MAX_MOI`** | a runaway guard that becomes the measured quantity if set low | confirm it never binds in published runs |
| S9 | **`MEAN_BLOOD_STAGE_DURATION`** | sets episode overlap and therefore polyclonality and outcrossing | sweep |
| S10 | **Scale (`N_HOSTS`, `N_MOSQ`)** | prevalence and COI are scale-dependent; a rate fitted on a small population does not transfer | always calibrate at publication scale |
| S11 | **`MEAN_RECRUD_DELAY`** | 10.5 d is a literature assumption: the one published fit fixes the rate with an informative prior, so no estimate from data exists; the delay sets how many recrudescences a 28-day follow-up sees | **pending parameterisation**; sweep (matters only with `P_RECRUD_GIVEN_FAIL > 0`) |

---

## Planned extensions

The extensions planned for later versions, each with the reason for it.

| extension | motivation |
|---|---|
| **Genotyping each hypnozoite independently of the inoculated primaries** | at present a bite's hypnozoites are clones of the genotypes it brings to the blood, so the reservoir is a subset of the blood-stage infection; this change makes it emergent. **The highest-priority change to the model.** |
| Structured import source pool | currently imports are unrelated founders, so import detection is an upper bound |
| Post-treatment prophylaxis | a course leaves no period of drug protection, so a treated host can be reinfected or relapse the next day, which overstates recurrence after treatment (section 11 of the model description); planned, with protection that wanes after each course |
| Calibration of the treatment and recurrence machinery | including a chloroquine-resistant setting, in which recrudescence is material rather than negligible |
| Founder genomes from recombuddy | root genomes painted with background relatedness and linkage disequilibrium, which independent painting leaves out (section 14.4 of the model description) |
| Burden-scaled inoculum | sporozoite inoculum tracks salivary-gland load (Kanatani et al. 2024); heavily-infected mosquitoes should seed proportionally more |
| Per-oocyst productivity | S4 |
| Origin-dependent detection | to separate active vs passive surveillance |
| Richer VCF output: `AF`, within-sample frequencies and phase | `write_vcf()` reports only `PLAF`, the frequency the root genomes were painted with, and collapses a host sample to unphased `0/0` / `0/1` / `1/1`, which says only that both alleles are present. Each clone's allele is recorded, so a host sample's phase is known exactly and could be written as phased polyploid calls, one allele per clone, beside `AF` among the samples written and the true within-sample ALT frequency that mixed-infection deconvolution sets out to recover. Simulated read depth would let `AD` and `DP` make the files read like called data. |
| Sample and genome naming | clone samples are named by their bare pedigree id (`51`), which VCF allows but which makes an awkward column name downstream, and clone, genome, lineage and genotype each mean something different across the outputs; the most informative convention needs thinking through |
| Fate determination (`P_HYPNOZOITE`) | 0.58–0.63 of liver-stage parasites became hypnozoites in culture, depending on the patient isolate (Vantaux et al. 2022); **deferred** — would make ν emergent but trades one well-identified parameter for two poorly-identified ones |
| Immunity | explicitly out of scope for now; would require a separate identifiability argument |
| Spatial structure / heterogeneous biting | currently homogeneous mixing |
| Vector fitness cost of burden | `MAX_LINEAGES_PER_MOSQ` is a hard cap with no graded effect; mosquito bottleneck simulations may be revisited |
