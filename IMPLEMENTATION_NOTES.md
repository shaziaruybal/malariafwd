# malariafwd v0.1.0 — implementation notes and traps

Traps for anyone who changes the model script
(`inst/slim/malariafwd_vivax.slim`) or reads its tree sequence. The
[model
description](https://shaziaruybal.github.io/malariafwd/MODEL_DESCRIPTION.html)
says what the model does; [design
decisions](https://shaziaruybal.github.io/malariafwd/DECISIONS.html)
says why.

**Eidos scoping.** A plain assignment to a global *inside a callback or
function creates a local*. Use `defineGlobal()`. The trap is silent: a
list filled this way inside a callback never reaches `early()`, which is
why `RECRUD_RETIRE`, the genomes a recrudescence replaces, is kept with
`defineGlobal()`. A hypnozoite batch is retired differently, by setting
its own count to 0, which the validation suite’s over-draw check
([validation](https://shaziaruybal.github.io/malariafwd/VALIDATION.html))
watches: a batch left unretired would keep activating after it had run
out.

**Vacant tskit nodes.** SLiM haploid individuals get **two** nodes — the
real haplosome (`is_vacant = [0]`) and an isolated null placeholder
(`is_vacant = [1]`, IBD 0 with everything). ~50% of `ts.samples()` are
vacant and **must be filtered**. Any new tskit code must do this.

**The genealogy is a forest.** Founders have independent roots; ~93% of
pairs never coalesce. `tskit` branch-mode `divergence()` is therefore
**not** a relatedness measure here, and a time to the most recent common
ancestor (TMRCA) derived from it can come out older than the run itself.
Use
[`ibd_pairs()`](https://shaziaruybal.github.io/malariafwd/reference/ibd_pairs.md)
/
[`arg_tmrca()`](https://shaziaruybal.github.io/malariafwd/reference/arg_tmrca.md).

**Integer and float spellings.** `-d X=0` makes an Eidos integer and
`-d X=0.0` a float, and Eidos will not mix the two in some places
([`rbinom()`](https://rdrr.io/r/stats/Binomial.html), `seqLen()`,
subsetting, a tick range). The model converts a parameter to the type it
needs where it uses it there (`asInteger()`, `asFloat()`), and
`initialize()` checks each parameter’s type and range (D43), so a run
started with `slim -d` takes either spelling. A new parameter must be
converted at each such use (a constant given with `-d` cannot be
redefined in `initialize()`). R writes every float with a decimal point.

**Iterating `NULL`.** A `for` loop over a zero-length vector runs no
iteration, and the model relies on that, but a loop over `NULL` is an
error: guard anything that can be `NULL`.

**Chromosome-aware joins.** When splitting IBD segments by chromosome,
split each segment on its own — a whole-group comparison against the
14-row chromosome table silently recycles for pairs with more than one
segment.

**Painting alleles with tskit.** To paint the root genomes with tskit
itself, simplify with `keep_input_roots=True` (otherwise the roots are
the MRCAs of the genomes kept, and change with them), give each mutation
`metadata={"mutation_list": []}` (SLiM’s schema requires it), and read
`genotype_matrix(isolated_as_missing=False)`: by default a genome alone
under its root genome, such as an import that infected none of the
genomes kept, reads as missing wherever its root is REF.

**tskit’s segments are path pieces.** Read straight from tskit, a pair’s
shared stretches come in pieces that depend on the samples the tree was
simplified to, and a stretch can run across a chromosome junction. Merge
abutting pieces with the same ancestor and cut at the junctions before
counting or drawing them (sections 6.3 and 7 of the model description).
