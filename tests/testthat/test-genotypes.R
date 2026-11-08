## simulate_genotypes() and write_vcf(): genotypes painted on the root genomes of the genealogy.

## A run with co-transmission, imports and recrudescence, made once for this file.
genotypes_run <- local({
  r <- NULL
  function() {
    if (is.null(r))
      r <<- suppressMessages(run_sim(vivax_params(
        N_HOSTS = 80, N_MOSQ = 200, BITING_RATE = 0.11, SIM_DAYS = 365, N_FOUNDER_HOSTS = 20, N_FOUNDER_LINEAGES = 10,
        LAMBDA_OOCYSTS = 2, LAMBDA_SPOROZOITES = 2, IMPORT_START = 1L, IMPORT_END = 365L, IMPORT_RATE_PER_TICK = 0.05,
        BLOODCURE_COVERAGE = 0.6, P_BLOODCURE_SURVIVAL = 0.3, P_RECRUD_GIVEN_FAIL = 0.7, RADCURE_COVERAGE = 0.3,
        CHROMOSOMES = "1-14"), seed = 7, outdir = tempfile("mfwd_geno_")))
    r
  }
})
## loci spread over every chromosome of a run
some_loci <- function(r, per = 10, plaf = 0.35) {
  cm <- chrom_map(r)
  data.table::data.table(chrom = rep(cm$chrom, each = per),
                         pos = as.integer(unlist(lapply(cm$length, function(l) round(seq(1, l, length.out = per + 2)[-c(1, per + 2)])))),
                         plaf = plaf)
}
cm_full <- function() data.table::data.table(chrom = 1:14, length = malariafwd:::.PVP01_FULL,
                                             start = cumsum(c(0, utils::head(malariafwd:::.PVP01_FULL, -1))),
                                             end = cumsum(malariafwd:::.PVP01_FULL))

## ---- pure R ----

test_that("loci are read as PvP01 chromosome, 1-based position and ALT frequency, and checked", {
  gl <- malariafwd:::.genotype_loci; cm <- cm_full()
  lo <- gl(cm, data.frame(chrom = c("1", "PvP01_02_v2", "PvP01_02_v1"), pos = c(1, 1, 5), plaf = c(0, 0.5, 1), note = "x"))
  expect_identical(lo$chrom, c(1L, 2L, 2L))
  expect_identical(lo$position, c(0, 1021664, 1021668))
  expect_identical(lo$locus, c("1:1", "2:1", "2:5"))
  expect_identical(names(lo), c("locus", "chrom", "pos", "position", "plaf", "note"))
  core <- data.table::copy(cm)[, length := round(length * 0.9)]
  expect_error(gl(core, data.frame(chrom = 3, pos = 896704, plaf = 0.1)), "loci row 1: pos 896704 is beyond chromosome 3, which is 807034 bp in this run")
  sub <- cm[c(3, 5, 6)]
  expect_error(gl(sub, data.frame(chrom = 1, pos = 1, plaf = 0.1)), "chromosome 1 is not simulated in this run \\(CHROMOSOMES gives 3, 5, 6\\)")
  expect_error(gl(cm, data.frame(chrom = 15, pos = 1, plaf = 0.1)), "loci row 1: chrom 15 is not a PvP01 chromosome")
  expect_error(gl(cm, data.frame(chrom = 1, pos = 0, plaf = 0.1)), "loci row 1: pos must be a whole number from 1")
  expect_error(gl(cm, data.frame(chrom = c(2, 1, 2), pos = c(5, 5, 5), plaf = 0.1)), "loci rows 1 and 3 are both chromosome 2, pos 5")
  expect_error(gl(cm, data.frame(chrom = 1:2, pos = 1:2, plaf = 0.1, locus = "a")), "loci rows 1 and 2 have the same locus id a")
  for (v in list(NA, -0.1, 1.1)) expect_error(gl(cm, data.frame(chrom = 1, pos = 1, plaf = v)), "loci row 1: plaf must be a frequency")
  expect_error(gl(cm, data.frame(chrom = 1, pos = 1, PLAF = 0.2)), "use plaf = 1 - PLAF", fixed = TRUE)
  expect_error(gl(cm, data.frame(chrom = 1, plaf = 0.2)), "`loci` lacks the column\\(s\\) pos")
  expect_identical(nrow(gl(cm, data.frame(chrom = integer(), pos = integer(), plaf = numeric()))), 0L)
  ## other columns are kept, whatever their names
  ex <- gl(cm, data.frame(chrom = 1:3, pos = c(100, 200, 300), plaf = c(0.1, 0.2, 0.3), k = 3L, num = 7, pl = 2,
                          cm = c(0.5, 1.2, 3.4)))
  expect_identical(ex$position, c(99, 1021863, 1978290))
  expect_identical(c(ex$chrom, ex$plaf), c(1:3, 0.1, 0.2, 0.3))
  expect_identical(ex$locus, c("1:100", "2:200", "3:300"))
  expect_identical(c(ex$k, ex$num, ex$pl, ex$cm), c(3, 3, 3, 7, 7, 7, 2, 2, 2, 0.5, 1.2, 3.4))
})

test_that("host calls pool their clones: 0, 1, or 0.5 for both", {
  al <- matrix(c(0L, 0L, 1L, 1L,  0L, 1L, 1L, 1L), 4, 2, dimnames = list(NULL, c("l1", "l2")))
  p <- malariafwd:::.pool_calls(al, c("s2", "s2", "s1", "s3"))
  expect_identical(rownames(p), c("s2", "s1", "s3"))
  expect_identical(unname(p), matrix(c(0, 1, 1, 0.5, 1, 1), 3, 2))
  expect_identical(dim(malariafwd:::.pool_calls(al[0, , drop = FALSE], character())), c(0L, 2L))
})

test_that("errors replace a call and missingness blanks one, at their rates", {
  oc <- malariafwd:::.observe_calls
  tc <- matrix(rep(c(0, 1), 100), 20, 10); hc <- matrix(rep(c(0, 0.5, 1, 0), 50), 20, 10)
  expect_identical(oc(tc, "clone", 0, 0), tc)
  expect_true(all(is.na(oc(tc, "clone", 0, 1))))
  expect_identical(oc(tc, "clone", 1, 0), 1 - tc)
  h <- oc(hc, "host", 1, 0)
  expect_true(all(h != hc) && all(h %in% c(0, 0.5, 1)))
  set.seed(3); big <- matrix(0, 200, 200)
  o <- oc(big, "clone", 0.1, 0.2)
  m <- mean(is.na(o)); e <- mean(o[!is.na(o)] != 0)
  expect_lt(abs(m - 0.2), 4 * sqrt(0.2 * 0.8 / 40000))
  expect_lt(abs(e - 0.1), 4 * sqrt(0.1 * 0.9 / (40000 * 0.8)))
})

test_that("arguments are refused before Python, naming them", {
  d <- tempfile("mfwd_no_run_")
  expect_error(simulate_genotypes(file.path(d, "plasmo.trees"), 1, data.frame()), "not a .trees file", fixed = TRUE)
  r <- fake_run()
  lo <- data.frame(chrom = 1, pos = 1, plaf = 0.1)
  for (v in list(-0.1, 1.1, NA, c(0.1, 0.2))) {
    expect_error(simulate_genotypes(r, 1, lo, error = v), "`error` must be a single number from 0 to 1")
    expect_error(simulate_genotypes(r, 1, lo, missing = v), "`missing` must be a single number from 0 to 1")
  }
  for (v in list(1.5, NA, c(1, 2))) expect_error(simulate_genotypes(r, 1, lo, seed = v), "`seed` must be NULL or a single whole number")
  expect_error(simulate_genotypes(r, 1:2, lo, sample_id = "a"), "`sample_id` must be a vector as long as `ids`")
  expect_error(simulate_genotypes(r, 1:2, lo, sample_id = c("a", NA)), "`sample_id` must not be NA or empty")
  expect_error(simulate_genotypes(r, 1:2, lo, sample_id = c("a", "")), "`sample_id` must not be NA or empty")
  expect_error(simulate_genotypes(r, c(1, NA), lo), "`ids` must be pedigree ids")
  expect_error(simulate_genotypes(r, c(1, 1), lo), "each genome is genotyped once")
  expect_error(simulate_genotypes(r, c(1, 1), lo, sample_id = c("a", "a")), "a genome is in one host sample only once")
})

test_that("write_vcf() writes clone samples haploid and host samples diploid, sorted, with the run's contigs", {
  cm <- cm_full()[1:2]
  lo <- malariafwd:::.genotype_loci(cm, data.frame(chrom = c(2, 1, 2, 1), pos = c(50, 10, 20, 30), plaf = c(0.1, 0.2, 0.3, 1/3)))
  calls <- matrix(c(0, 1, NA,  1, 1, 0,  0, 0, 1,  1, NA, 0), 3, 4, dimnames = list(c("s1", "s2", "s3"), lo$locus))
  mk <- function(type, calls) malariafwd:::.new_vivax_genotypes(calls, calls, NULL, NULL, NULL, NULL, lo, NULL, NULL, cm,
                                                                type, 0.01, 0.05, 7L, "runs/x")
  f <- tempfile(fileext = ".vcf"); write_vcf(mk("clone", calls), f)
  v <- readLines(f)
  expect_identical(v[1], "##fileformat=VCFv4.2")
  expect_true(all(c("##contig=<ID=PvP01_01_v2,length=1021664>", "##contig=<ID=PvP01_02_v2,length=956327>",
                    "##malariafwd_seed=7", "##malariafwd_samples=clone") %in% v))
  expect_true(any(grepl("^##malariafwd_alleles=REF and ALT are placeholders", v)))
  hdr <- strsplit(v[startsWith(v, "#CHROM")], "\t")[[1]]
  expect_identical(utils::tail(hdr, 3), c("s1", "s2", "s3"))
  rec <- do.call(rbind, strsplit(v[!startsWith(v, "#")], "\t"))
  expect_identical(rec[, 1], c("PvP01_01_v2", "PvP01_01_v2", "PvP01_02_v2", "PvP01_02_v2"))
  expect_identical(rec[, 2], c("10", "30", "20", "50"))
  o <- order(lo$chrom, lo$pos)
  back <- matrix(c("0" = 0, "1" = 1, "." = NA)[t(rec[, 10:12])], 3)
  expect_identical(unname(back), unname(calls[, o]))
  expect_lt(max(abs(as.numeric(sub("PLAF=", "", rec[, 8])) - lo$plaf[o])), 1e-12)
  expect_identical(unique(rec[, 4:5]), matrix(c("A", "C"), 1))
  ## host samples, diploid
  hc <- calls; hc[1, 1] <- 0.5
  fh <- tempfile(fileext = ".vcf"); write_vcf(mk("host", hc), fh)
  rh <- do.call(rbind, strsplit(readLines(fh)[!startsWith(readLines(fh), "#")], "\t"))
  expect_true(all(c("0/1", "0/0", "1/1", "./.") %in% rh[, 10:12]))
  ## the true calls, given alleles, numbered contigs
  lo2 <- data.table::copy(lo)[, `:=`(ref = "G", alt = "T")]
  x2 <- malariafwd:::.new_vivax_genotypes(calls, matrix(0, 3, 4, dimnames = dimnames(calls)), NULL, NULL, NULL, NULL, lo2, NULL, NULL, cm, "clone", 0, 0, NA_integer_, "r")
  f2 <- tempfile(fileext = ".vcf"); write_vcf(x2, f2, calls = "true", chrom_names = as.character(1:14))
  v2 <- readLines(f2)
  expect_true("##malariafwd_calls=true" %in% v2); expect_true("##malariafwd_seed=NA" %in% v2)
  expect_false(any(grepl("placeholders", v2)))
  r2 <- do.call(rbind, strsplit(v2[!startsWith(v2, "#")], "\t"))
  expect_identical(unique(r2[, 1]), c("1", "2")); expect_identical(unique(r2[, 4:5]), matrix(c("G", "T"), 1))
  expect_true(all(r2[, 10:12] == "0"))
  ## refusals
  lo3 <- data.table::copy(lo)[, `:=`(ref = "AT", alt = "T")]
  expect_error(write_vcf(malariafwd:::.new_vivax_genotypes(calls, calls, NULL, NULL, NULL, NULL, lo3, NULL, NULL, cm, "clone", 0, 0, NA_integer_, "r"), f2),
               "`ref` and `alt` must each be one of A, C, G, T or N")
  for (cn in list(as.character(1:13), c(NA, 2:14), rep("a", 14), c("a b", 2:14)))
    expect_error(write_vcf(x2, f2, chrom_names = cn), "`chrom_names` must be 14 distinct contig names")
  expect_error(write_vcf(list(), f2), "`x` must be a vivax_genotypes object, from simulate_genotypes\\(\\), not a list")
  bad <- calls; rownames(bad)[1] <- "a b"
  expect_error(write_vcf(mk("clone", bad), f2), "a sample id contains whitespace")
  lo4 <- data.table::copy(lo)[1, locus := "x;y"]
  expect_error(write_vcf(malariafwd:::.new_vivax_genotypes(calls, calls, NULL, NULL, NULL, NULL, lo4, NULL, NULL, cm, "clone", 0, 0, NA_integer_, "r"), f2),
               "a locus id contains whitespace or ';'")
  ## no samples: sites only; no loci: the header alone; gzip reads back to the same lines
  f0 <- tempfile(fileext = ".vcf"); write_vcf(mk("clone", calls[0, , drop = FALSE]), f0)
  r0 <- strsplit(readLines(f0)[!startsWith(readLines(f0), "##")], "\t")
  expect_true(all(lengths(r0) == 8L))
  fl <- tempfile(fileext = ".vcf")
  x0 <- malariafwd:::.new_vivax_genotypes(calls[, 0, drop = FALSE], calls[, 0, drop = FALSE], NULL, NULL, NULL, NULL, lo[0], NULL, NULL, cm, "clone", 0, 0, NA_integer_, "r")
  write_vcf(x0, fl)
  expect_identical(sum(!startsWith(readLines(fl), "#")), 0L)
  fz <- tempfile(fileext = ".vcf.gz"); write_vcf(mk("clone", calls), fz)
  expect_identical(readLines(gzfile(fz)), v)
})

## ---- on a run ----

test_that("every genome carries its root genome's allele, by three routes", {
  skip_if_no_slim(); skip_if_no_tskit()
  r <- genotypes_run(); ev <- load_events(r); cs <- classify_cases(ev)
  rel <- ev[event_type == "hypnozoite_activation"][seq_len(min(.N, 15))]
  rec <- ev[event_type == "recrudescence"][seq_len(min(.N, 10))]
  ids <- unique(c(utils::head(cs$pedigree_id, 60), rel$pedigree_id, rel$parent1_pedigree_id, rec$pedigree_id,
                  rec$parent1_pedigree_id))
  lo <- some_loci(r)
  g <- simulate_genotypes(r, ids, lo, seed = 1)
  ## 1 the bookkeeping
  k <- match(g$roots, g$root_genomes$root_id)
  expect_false(anyNA(k))
  expect_identical(as.vector(g$alleles), g$root_alleles[cbind(k, as.vector(col(g$roots)))])
  ## 2 tskit, independently: mutations on the root nodes, then genotype_matrix()
  py <- reticulate::import_main(convert = FALSE)
  m <- malariafwd:::.py_module("_malariafwd_genotype_test", '
def paint(path, ids, rid, A, positions):
    ts = tskit.load(path); p2n, _ = _ped_nodes(ts)
    sts = ts.simplify(samples=[p2n[int(p)] for p in ids], keep_input_roots=True)
    _, sped = _ped_nodes(sts)
    row = {int(p): i for i, p in enumerate(rid)}
    A = np.asarray(A); pos = np.asarray(positions, float)
    t = sts.dump_tables(); tree = sts.first()
    order = np.argsort(pos, kind="stable")
    for j in order:
        tree.seek(pos[j]); s = t.sites.add_row(position=pos[j], ancestral_state="0")
        for r in tree.roots:
            if A[row[int(sped[r])], j] == 1:
                t.mutations.add_row(site=s, node=r, derived_state="1", metadata={"mutation_list": []})
    t.sort(); t.build_index(); t.compute_mutation_parents()
    G = t.tree_sequence().genotype_matrix(isolated_as_missing=False)
    out = np.empty((len(ids), len(pos)), np.int32); out[:, order] = G.T
    return out
def walk(path, ids, positions):
    ts = tskit.load(path); p2n, ped = _ped_nodes(ts)
    W = np.empty((len(ids), len(positions)), np.int64); t = ts.first()
    for j, x in enumerate(positions):
        t.seek(float(x))
        for i, p in enumerate(ids):
            u = p2n[int(p)]
            while t.parent(u) != tskit.NULL: u = t.parent(u)
            W[i, j] = ped[u]
    return W
')
  G <- reticulate::py_to_r(m$paint(trees_path(r), as.list(ids), as.list(g$root_genomes$root_id), g$root_alleles,
                                   as.list(g$loci$position)))
  expect_identical(unname(g$alleles), G)
  ## 3 walking up the unsimplified genealogy reaches the root genome named
  W <- reticulate::py_to_r(m$walk(trees_path(r), as.list(ids), as.list(g$loci$position)))
  expect_true(all(unname(g$roots) == W))
  ## 4 root genomes alone are their own roots
  rid <- g$root_genomes$root_id[c(1, which(g$root_genomes$imported)[1:2])]
  g0 <- simulate_genotypes(r, rid, lo, seed = 1)
  expect_true(all(g0$roots == matrix(rid, length(rid), nrow(lo))))
  expect_identical(unname(g0$alleles), unname(g0$root_alleles[as.character(rid), ]))
})

test_that("two genomes share a root genome exactly where they share an IBD segment, and clones are identical", {
  skip_if_no_slim(); skip_if_no_tskit()
  r <- genotypes_run(); ev <- load_events(r)
  rel <- ev[event_type == "hypnozoite_activation"][seq_len(min(.N, 5))]
  rec <- ev[event_type == "recrudescence"][seq_len(min(.N, 3))]
  oo <- ev[event_type == "vh_inoculation" & is_dormant == 0 & oocyst_id != -1, .(pedigree_id, oocyst_id, lineage_id)]
  sib <- oo[, if (data.table::uniqueN(lineage_id) > 1) .(a = pedigree_id[1], b = pedigree_id[which(lineage_id != lineage_id[1])[1]]),
            by = oocyst_id][seq_len(min(.N, 5))]
  lo <- some_loci(r, per = 20)
  g0 <- simulate_genotypes(r, integer(), lo)
  roots <- g0$root_genomes$root_id
  pr <- data.table::rbindlist(list(rel[, .(a = pedigree_id, b = parent1_pedigree_id)], rec[, .(a = pedigree_id, b = parent1_pedigree_id)],
                                   sib[, .(a, b)], data.table::data.table(a = roots[1:2], b = roots[3:4])))
  g <- simulate_genotypes(r, unique(c(pr$a, pr$b)), lo, seed = 2)
  sg <- ibd_segments(r, pr$a, pr$b, min_span = 0)
  for (i in seq_len(nrow(pr))) {
    same <- unname(g$roots[as.character(pr$a[i]), ] == g$roots[as.character(pr$b[i]), ])
    s <- sg[pair == i]
    inside <- vapply(g$loci$position, function(x) any(s$gstart <= x & x < s$gend), logical(1))
    expect_identical(same, inside, label = paste("pair", i))
  }
  ## a relapse and a recrudescence have their origin's alleles; distinct root genomes share none
  for (d in list(rel, rec))
    expect_identical(unname(g$alleles[as.character(d$pedigree_id), ]), unname(g$alleles[as.character(d$parent1_pedigree_id), ]))
  expect_false(any(g$roots[as.character(roots[1]), ] == g$roots[as.character(roots[3]), ]))
})

test_that("root genomes are the founder templates and the imports, with their days", {
  skip_if_no_slim(); skip_if_no_tskit()
  r <- genotypes_run(); ev <- load_events(r)
  rg <- simulate_genotypes(r, integer(), some_loci(r))$root_genomes
  seeds <- ev[event_type == "vh_inoculation_seed"]; imp <- ev[event_type == "vh_inoculation_import"]
  expect_true(all(rg$root_id %in% c(seeds$pedigree_id, imp$pedigree_id)))
  expect_false(anyNA(rg$imported))
  expect_identical(rg$imported, rg$root_id %in% imp$pedigree_id)
  expect_identical(rg$tick, c(seeds$tick, imp$tick)[match(rg$root_id, c(seeds$pedigree_id, imp$pedigree_id))])
  expect_identical(sum(!rg$imported), 10L)
})

test_that("a seed gives the same root alleles whatever else is asked, and leaves the session's stream alone", {
  skip_if_no_slim(); skip_if_no_tskit()
  r <- genotypes_run(); cs <- classify_cases(load_events(r)); lo <- some_loci(r)
  ids <- utils::head(cs$pedigree_id, 30); more <- unique(c(ids, cs$pedigree_id[31:60]))
  a <- simulate_genotypes(r, ids, lo, seed = 10)
  set.seed(99); before <- .Random.seed
  b <- simulate_genotypes(r, rev(more), lo, sample_id = paste0("s", seq_along(more)), error = 0.1, missing = 0.2, seed = 10)
  expect_identical(.Random.seed, before)
  expect_identical(b$root_alleles, a$root_alleles); expect_identical(b$root_genomes, a$root_genomes)
  expect_identical(b$alleles[as.character(ids), ], a$alleles); expect_identical(b$roots[as.character(ids), ], a$roots)
  set.seed(10); c0 <- simulate_genotypes(r, ids, lo)
  expect_identical(c0$alleles, a$alleles); expect_identical(c0$seed, NA_integer_); expect_identical(a$seed, 10L)
})

test_that("host samples pool the clones of a bite, a genome can be in two samples, and the VCF reads back", {
  skip_if_no_slim(); skip_if_no_tskit()
  r <- genotypes_run(); ev <- load_events(r); lo <- some_loci(r)
  inoc <- ev[event_type == "vh_inoculation" & is_dormant == 0]
  multi <- inoc[, if (.N > 1) .(pedigree_id), by = inoculation_id][seq_len(min(.N, 60))]
  sid <- c(paste0("bite", multi$inoculation_id), "both", "both")
  ids <- c(multi$pedigree_id, multi$pedigree_id[1:2])
  g <- simulate_genotypes(r, ids, lo, sample_id = sid, seed = 3)
  expect_identical(g$type, "host")
  expect_identical(g$samples$n_clones, as.integer(table(factor(sid, unique(sid)))))
  expect_true(all(g$true_calls %in% c(0, 0.5, 1)))
  expect_gt(sum(g$true_calls == 0.5), 0)
  al <- g$alleles[as.character(ids), , drop = FALSE]
  expect_identical(g$true_calls, malariafwd:::.pool_calls(al, sid))
  expect_identical(nrow(g$alleles), length(unique(ids)))
  f <- tempfile(fileext = ".vcf"); write_vcf(g, f)
  v <- data.table::fread(f, skip = "#CHROM", sep = "\t", colClasses = "character")
  gt <- t(as.matrix(v[, -(1:9)]))
  back <- matrix(c("0/0" = 0, "1/1" = 1, "0/1" = 0.5, "./." = NA)[gt], nrow(gt))
  o <- order(g$loci$chrom, g$loci$pos)
  expect_identical(unname(back), unname(g$calls[, o]))
  expect_output(print(g), paste0(nrow(g$calls), " host samples \\(", length(unique(ids)), " genomes\\)"))
})

test_that("genomes not in the tree sequence, and runs that never started, are refused with the reason", {
  skip_if_no_slim(); skip_if_no_tskit()
  r <- genotypes_run(); lo <- some_loci(r)
  expect_error(simulate_genotypes(r, c(99999991L, 99999992L), lo), "pedigree id\\(s\\) 99999991, 99999992 are not genomes in this run's tree sequence")
  nv <- suppressMessages(run_sim(vivax_params(N_HOSTS = 20, N_MOSQ = 40, SIM_DAYS = 30, N_FOUNDER_HOSTS = 0, IMPORT_TICK = 5,
                                              N_IMPORTS = 1, P_PRIMARY = 0, NU_HYPNOZOITE_BATCH = 0, CHROMOSOMES = "1"),
                                 seed = 1, outdir = tempfile("mfwd_geno_never_")))
  expect_true(nv$never_started)
  lo1 <- data.frame(chrom = 1, pos = 1:3, plaf = 0.3)
  expect_error(simulate_genotypes(nv, 1L, lo1), "this run never started")
  e <- simulate_genotypes(nv, integer(), lo1)
  expect_identical(c(dim(e$calls), nrow(e$root_genomes)), c(0L, 3L, 0L))
})

test_that("the object keeps one schema, also with no genome or no locus", {
  skip_if_no_slim(); skip_if_no_tskit()
  r <- genotypes_run(); cs <- classify_cases(load_events(r)); lo <- some_loci(r)
  g <- simulate_genotypes(r, utils::head(cs$pedigree_id, 5), lo, seed = 1)
  nm <- c("calls", "true_calls", "alleles", "roots", "samples", "clones", "loci", "root_genomes", "root_alleles",
          "chromosomes", "type", "error", "missing", "seed", "run")
  for (x in list(g, simulate_genotypes(r, integer(), lo), simulate_genotypes(r, utils::head(cs$pedigree_id, 5), lo[0]))) {
    expect_s3_class(x, "vivax_genotypes")
    expect_identical(names(x), nm)
    expect_identical(c(storage.mode(x$calls), storage.mode(x$true_calls), storage.mode(x$alleles), storage.mode(x$roots),
                       storage.mode(x$root_alleles)), c("double", "double", "integer", "integer", "integer"))
    expect_identical(names(x$samples), c("sample_id", "n_clones"))
    expect_identical(names(x$clones), c("sample_id", "pedigree_id"))
    expect_identical(names(x$root_genomes), c("root_id", "tick", "imported"))
    expect_identical(names(x$loci)[1:5], c("locus", "chrom", "pos", "position", "plaf"))
    expect_identical(ncol(x$calls), nrow(x$loci))
  }
})
