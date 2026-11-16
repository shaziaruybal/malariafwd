## Invariants the model must keep. These need SLiM, so they skip when it is
## unavailable (e.g. on CRAN / CI without the engine) via skip_if_no_slim() in helper-slim.R.

## A full-genome run from 8 founder lineages, made once for this file.
invariants_run <- local({
  r <- NULL
  function() {
    if (is.null(r))
      r <<- suppressMessages(run_sim(vivax_params(N_HOSTS = 60, N_MOSQ = 120, BITING_RATE = 0.12, SIM_DAYS = 300,
                                                  N_FOUNDER_LINEAGES = 8),
                                     seed = 1, outdir = tempfile("mfw-invariants")))
    r
  }
})

test_that("hypnozoite batches never draw more than they hold, and ids are unique", {
  skip_if_no_slim()
  ev <- load_events(invariants_run())

  ## a batch's activations + deaths must never exceed the size it was established with:
  ## a retired batch draws nothing more
  est <- ev[grepl("^vh_inoculation", event_type) & is_dormant == 1,
            .(batch = pedigree_id, held = n_hypnozoites)]
  act <- ev[event_type == "hypnozoite_activation", .(batch = parent1_pedigree_id, n = n_hypnozoites)]
  die <- ev[grepl("hypnozoite_clear", event_type), .(batch = pedigree_id, n = n_hypnozoites)]
  drawn <- rbind(act, die)[, .(drew = sum(n)), by = batch]
  chk <- merge(est, drawn, by = "batch", all.x = TRUE)[is.na(drew), drew := 0]
  expect_equal(sum(chk$drew > chk$held), 0)

  ## event ids are unique, and inoculation ids do not restart each tick
  expect_equal(data.table::uniqueN(ev$event_id), nrow(ev))
  inoc <- ev[inoculation_id > 0, inoculation_id]
  expect_gt(max(inoc), 0.5 * data.table::uniqueN(inoc))
})

test_that("zero-truncated count helpers round-trip", {
  ## lambda is the PRE-truncation rate: the realised mean is higher
  expect_gt(uptake_moments(2)[["mean"]], 2)
  expect_equal(unname(uptake_moments(2)[["mean"]]), 2 / (1 - exp(-2)), tolerance = 1e-8)
  ## NB adds dispersion at the same lambda
  expect_gt(uptake_moments(2, k = 1)[["var"]], uptake_moments(2, k = 0)[["var"]])
  ## inverse recovers the targets
  p <- uptake_calibrate(2.5)
  expect_equal(unname(p$achieved[["mean"]]), 2.5, tolerance = 1e-6)
  expect_equal(p$k, 0)
  q <- uptake_calibrate(2.5, var = 6)
  expect_equal(unname(q$achieved[["mean"]]), 2.5, tolerance = 1e-3)
  expect_equal(unname(q$achieved[["var"]]),  6.0, tolerance = 1e-3)
  expect_gt(q$k, 0)
  expect_error(uptake_calibrate(0.8), "mean > 1")
})

test_that("ARG summaries behave on the forest genealogy", {
  skip_if_no_slim(); skip_if_no_tskit()
  run <- invariants_run()
  ev <- load_events(run)
  s <- arg_summary(run)
  expect_gt(s$num_trees, 1)          # recombination happened
  expect_gte(s$roots_min, 1)
  ## the genealogy is a FOREST: independent founder lineages never coalesce
  expect_gt(s$roots_median, 1)

  ## arg_tmrca is the relatedness measure that works on a forest
  rel <- ev[event_type == "hypnozoite_activation" & parent1_pedigree_id != -1][1:5]
  if (nrow(rel) >= 1 && !any(is.na(rel$pedigree_id))) {
    tm <- arg_tmrca(run, rel$pedigree_id, rel$parent1_pedigree_id)
    expect_true(all(tm$coalesced > 0.99, na.rm = TRUE))   # a relapse is clonal to its origin
  }
  fnd <- unique(ev[event_type == "vh_inoculation_seed"], by = "lineage_id")$pedigree_id
  if (length(fnd) >= 4) {
    tf <- arg_tmrca(run, fnd[1:2], fnd[3:4])
    expect_true(all(tf$coalesced < 0.5, na.rm = TRUE))    # distinct founders are independent roots
  }

  ## D17: tmrca is a DEPTH, so one fixed relationship gives one number wherever it was
  ## sampled. The ancestor's raw node time would instead change with the day the pair was sampled.
  spread <- ev[event_type == "hypnozoite_activation" & parent1_pedigree_id != -1]
  if (nrow(spread) >= 8) {
    spread <- spread[order(tick)][unique(round(seq(1, .N, length.out = 8)))]
    td <- arg_tmrca(run, spread$pedigree_id, spread$parent1_pedigree_id)
    td <- td[!is.na(td$tmrca), ]
    expect_true(all(td$tmrca >= 0))                       # an ancestor predates its descendants
    expect_lt(diff(range(td$tmrca)), 1e-6)                # identical relationship, identical depth
    expect_true("tmrca_node_time" %in% names(td))         # the raw node time is still available
  }
})

test_that("a new genome's coordinates are 0, also where SLiM reuses a dead genome's memory (D7)", {
  skip_if_no_slim()
  src <- readLines(slim_model(), warn = FALSE)
  fun <- function(name) {
    i <- grep(paste0("^function \\(.*\\) ", name, "\\("), src)
    src[i:(i + which(src[(i + 1):length(src)] == "}")[1])]
  }
  ## SLiM gives a new individual the memory of one that died, coordinates included: the script gives
  ## 40 genomes coordinates, kills them, makes 40 more and passes them through the model's two
  ## functions that make a genome, 20 each
  script <- tempfile(fileext = ".slim")
  writeLines(c(fun("setHostGenome"), fun("tagMosqProduct"), "initialize() {",
    "  initializeSLiMModelType('nonWF'); initializeSLiMOptions(keepPedigrees=T);",
    "  defineConstant('K_LATENCY', 0); defineConstant('GAMETOCYTE_DELAY', 0); defineConstant('EIP_DAYS', 10);",
    "  initializeChromosome(1, 100, 'H'); initializeMutationType('m1', 0.5, 'f', 0.0);",
    "  initializeGenomicElementType('g1', m1, 1.0); initializeGenomicElement(g1, 0, 99);",
    "  initializeMutationRate(0.0); initializeRecombinationRate(0.0); }",
    "1 early() { old = sim.addSubpop('p1', 40).individuals; old.x = 1.5; old.y = 2.5; old.z = 3.5;",
    "  sim.killIndividuals(old); }",
    "2 early() { inds = sim.addSubpop('p2', 40).individuals;",
    "  catn('BEFORE ' + paste(c(inds.x, inds.y, inds.z)));",
    "  for (g in inds[0:19]) setHostGenome(g, birth=2, lineage=g.pedigreeID, oocyst=-1, brood=-1, founder=T,",
    "    inoc_tick=2, inoc_id=-1, imported=F, activation=NULL, batch=0, color='grey');",
    "  for (g in inds[20:39]) tagMosqProduct(g, g.pedigreeID, -1, -1, 'grey');",
    "  catn('AFTER ' + paste(c(inds.x, inds.y, inds.z))); sim.simulationFinished(); }"), script)
  res <- processx::run(slim_path(), script, error_on_status = FALSE)
  expect_equal(res$status, 0L, info = res$stderr)
  read <- function(w) as.numeric(strsplit(sub(paste0("^", w, " "), "",
    grep(paste0("^", w, " "), strsplit(res$stdout, "\n")[[1]], value = TRUE)), " ")[[1]])
  ## without memory reuse the check below would prove nothing
  skip_if(all(read("BEFORE") == 0), "this SLiM did not reuse the dead genomes' memory")
  expect_length(read("AFTER"), 120L)
  expect_true(all(read("AFTER") == 0))
})

test_that("two runs of one seed give the same tree-sequence tables, provenance aside (D7)", {
  skip_if_no_slim(); skip_if_no_tskit()
  ## founders and their clonal copies, bites, mosquito products, relapses, imports and recrudescences
  p <- vivax_params(N_HOSTS = 40, N_MOSQ = 120, BITING_RATE = 0.2, SIM_DAYS = 120, CHROMOSOMES = "1",
                    N_FOUNDER_LINEAGES = 4, IMPORT_START = 1, IMPORT_END = 120, IMPORT_RATE_PER_TICK = 0.3,
                    BLOODCURE_COVERAGE = 0.5, P_BLOODCURE_SURVIVAL = 0.5, P_RECRUD_GIVEN_FAIL = 1)
  a <- run_sim(p, seed = 2, outdir = tempfile("mfw-xyz-a"))
  b <- run_sim(p, seed = 2, outdir = tempfile("mfw-xyz-b"))
  m <- malariafwd:::.py_module("_malariafwd_repeat_test", '
def same_tables(x, y):
    return tskit.load(x).dump_tables().equals(tskit.load(y).dump_tables(), ignore_provenance=True)
def locations(path):
    return tskit.load(path).tables.individuals.location.tolist()
')
  expect_true(reticulate::py_to_r(m$same_tables(trees_path(a), trees_path(b))))
  loc <- unlist(reticulate::py_to_r(m$locations(trees_path(a))))
  expect_gt(length(loc), 0)
  expect_true(all(loc == 0))
})

test_that("an inoculation's primary and its hypnozoite batch are the same genome", {
  skip_if_no_slim(); skip_if_no_tskit()
  ## D38: an import is seeded in early() and its second compartment is cloned one tick later,
  ## so an import, like an ordinary bite, must give IBD 1 between its two compartments.
  out <- tempfile("mfw-import-clone")
  ## imports on, so both the founder path and the ordinary path are exercised
  p <- vivax_params(N_HOSTS = 80, N_MOSQ = 160, BITING_RATE = 0.10, SIM_DAYS = 500,
                    N_FOUNDER_LINEAGES = 8, NU_HYPNOZOITE_BATCH = 5,
                    IMPORT_START = 1L, IMPORT_END = 500L, IMPORT_RATE_PER_TICK = 0.08)
  run <- run_sim(p, seed = 3, outdir = out)
  ev <- load_events(run)

  ## One inoculation delivers one sporozoite population, so the blood-stage infection and the
  ## hypnozoite batch established alongside it must be clones. Two unrelated genomes under one
  ## lineage_id would not show in the event log alone, which is why this check goes through the
  ## tree sequence.
  ## grouped on inoculation_id, NOT (tick, human_id): under D38 an import's two compartments
  ## land on CONSECUTIVE ticks, because the batch cannot be cloned until the primary is
  ## visible. inoculation_id is the thing that identifies one bite in both paths.
  pairs <- ev[grepl("^vh_inoculation", event_type) & pedigree_id >= 0 & inoculation_id > 0,
              .(pedigree_id, is_dormant, event_type), by = inoculation_id]
  pairs <- pairs[, if (uniqueN(is_dormant) == 2L)
                     .(primary = pedigree_id[is_dormant == 0][1],
                       batch   = pedigree_id[is_dormant == 1][1],
                       kind    = event_type[1]), by = inoculation_id]
  skip_if(nrow(pairs) == 0, "no inoculation produced both a primary and a batch")

  take <- pairs[seq_len(min(.N, 25))]
  x <- ibd_pairs(run, pairs = take[, .(a = primary, b = batch)])
  take[, ibd := ibd_lookup(x, primary, batch)]
  expect_equal(take$ibd, rep(1, nrow(take)), tolerance = 1e-8)

  ## and the import path on its own, since it makes the two compartments on different ticks
  imp <- take[grepl("import", kind)]
  if (nrow(imp)) expect_true(all(imp$ibd == 1))
})

test_that("the oocyst count is drawn exactly where a positive count is rare", {
  skip_if_no_slim()
  ## The model's drawZTCount(), with its rejection loop skipped, so that every draw is the exact
  ## draw it makes once 10000 tries have found no positive count.
  src <- readLines(slim_model(), warn = FALSE)
  i <- grep("^function \\(integer\\$\\) drawZTCount", src)
  j <- i + which(src[(i + 1):length(src)] == "}")[1]
  fn <- sub("seqLen(10000)", "seqLen(0)", src[i:j], fixed = TRUE)
  expect_true(any(grepl("seqLen(0)", fn, fixed = TRUE)))
  for (case in list(c(1.0153e-4, 0.001), c(0.05, 0.05), c(0.3, 2))) {
    script <- tempfile(fileext = ".slim")
    writeLines(c(fn, "initialize() {", "  setSeed(4);",
                 sprintf("  catn('X ' + paste(sapply(seqLen(20000), 'drawZTCount(%.17g, %.17g);')));", case[1], case[2]),
                 "}", "1 early() { sim.simulationFinished(); }"), script)
    r <- processx::run(slim_path(), script, error_on_status = FALSE)
    expect_equal(r$status, 0L)
    x <- scan(text = sub("^X ", "", grep("^X ", strsplit(r$stdout, "\n")[[1]], value = TRUE)), quiet = TRUE)
    want <- uptake_moments(case[1], case[2])
    lab <- sprintf("lambda %g, k %g", case[1], case[2])
    expect_lt(abs(mean(x) - want[["mean"]]), 5 * sd(x) / sqrt(length(x)), label = paste("mean,", lab))
    expect_lt(abs(var(x) - want[["var"]]), 5 * sd((x - mean(x))^2) / sqrt(length(x)), label = paste("variance,", lab))
  }
})

test_that("an outcrossed oocyst makes four products by four-strand meiosis with meiosis I pairing", {
  skip_if_no_slim()
  ## The model's own seedTetrad() and tagMosqProduct(), run on the PvP01 map against two parents
  ## told apart by a marker at every probe locus: 20 probes per chromosome and one at each
  ## centromere, 400 tetrads.
  src <- readLines(slim_model(), warn = FALSE)
  fun <- function(name) {
    i <- grep(paste0("^function \\(.*\\) ", name, "\\("), src)
    src[i:(i + which(src[(i + 1):length(src)] == "}")[1])]
  }
  def <- paste(src, collapse = " ")
  constant <- function(name)
    as.numeric(strsplit(gsub("[^0-9.,]", "", sub(paste0(".*defineConstant\\(\"", name, "\", c\\(([^)]*)\\).*"), "\\1", def)), ",")[[1]])
  lens <- constant("PVP01_FULL"); frac <- constant("PVP01_CEN_FRAC")
  expect_length(lens, 14)
  expect_length(frac, 14)
  r <- 7.4e-7
  starts <- c(0, cumsum(lens)[-14])
  cens <- starts + round(frac * lens)
  chr <- rep(1:14, each = 21)
  probes <- unlist(lapply(1:14, function(i) sort(c(starts[i] + round(lens[i] * (1:20) / 21), cens[i]))))
  at_cen <- match(cens, probes)
  ## the map, CHROM_STARTS and CHROM_CENTROMERES come from the model's own GENOME code
  genome <- src[grep('defineConstant("PVP01_FULL"', src, fixed = TRUE):grep('defineConstant("RECOMB_WITHIN"', src, fixed = TRUE)]
  genome_init <- function(mode, chromosomes)
    c(sprintf("  defineConstant('GENOME_MODE', '%s'); defineConstant('CHROMOSOMES', '%s');", mode, chromosomes),
      "  defineConstant('RECOMB_RATE', -1.0);", genome,
      "  initializeMutationType('m1', 0.5, 'f', 0.0); m1.convertToSubstitution = F;",
      "  initializeChromosome(1, GENOME_LENGTH, 'H'); initializeGenomicElementType('g1', m1, 1.0);",
      "  initializeGenomicElement(g1, 0, GENOME_LENGTH - 1); initializeMutationRate(0.0);",
      "  initializeRecombinationRate(rec_rates, rec_ends);")
  ## core mode and a subset place each centromere at the same fraction of its chromosome
  chk <- tempfile(fileext = ".slim")
  writeLines(c("initialize() {", "  initializeSLiMModelType('nonWF');", genome_init("core", "3,5-6"),
               "  catn('CEN ' + paste(CHROM_CENTROMERES)); }", "1 early() { sim.simulationFinished(); }"), chk)
  res <- processx::run(slim_path(), chk, error_on_status = FALSE)
  expect_equal(res$status, 0L, info = res$stderr)
  core <- round(lens[c(3, 5, 6)] * 0.9)
  got <- as.numeric(strsplit(sub("^CEN ", "", grep("^CEN ", strsplit(res$stdout, "\n")[[1]], value = TRUE)), " ")[[1]])
  expect_true(all(abs(got - (c(0, cumsum(core)[-3]) + round(frac[c(3, 5, 6)] * core))) <= 1))   # rounding
  script <- tempfile(fileext = ".slim")
  writeLines(c(fun("tagMosqProduct"), fun("seedTetrad"), "initialize() {",
    "  initializeSLiMModelType('nonWF'); initializeSLiMOptions(keepPedigrees=T); setSeed(11);",
    genome_init("full", "1-14"),
    "  defineConstant('EIP_DAYS', 10);",
    sprintf("  defineConstant('PROBES', c(%s));", paste(probes, collapse = ", ")), "}",
    "1 first() { sim.addSubpop('p1', 2); sim.addSubpop('p2', 0);",
    "  p1.individuals[0].haplosomes.addNewMutation(m1, 0.0, PROBES);",
    "  p1.individuals[0].setValue('lineage_id', -1); p1.individuals[1].setValue('lineage_id', -2); }",
    "reproduction(p1) { if (individual.index == 0) for (t in seqLen(400)) seedTetrad(p2, p1.individuals[0], p1.individuals[1], 'x', t, t); }",
    "1 late() { for (ind in p2.individuals) catn('P ' + ind.getValue('oocyst_id') + ' ' + ind.getValue('lineage_id') + ' ' +",
    "  ind.pedigreeID + ' ' + paste(asInteger(match(PROBES, ind.haplosomes.mutations.position) >= 0), sep = '')); sim.simulationFinished(); }"), script)
  res <- processx::run(slim_path(), script, error_on_status = FALSE)
  expect_equal(res$status, 0L, info = res$stderr)
  L <- grep("^P ", strsplit(res$stdout, "\n")[[1]], value = TRUE)
  x <- do.call(rbind, strsplit(L, " "))
  x <- x[order(as.numeric(x[, 4])), ]                               # each oocyst's products in the order made
  oocyst <- as.integer(x[, 2]); lineage <- as.integer(x[, 3])
  g <- do.call(rbind, lapply(strsplit(x[, 5], ""), as.integer))     # 1 = parent1's allele
  expect_equal(nrow(g), 1600)
  expect_equal(ncol(g), length(probes))
  tetrads <- split(seq_len(nrow(g)), oocyst)
  ## at every locus two of the four products carry each parent's allele
  expect_true(all(sapply(tetrads, function(i) colSums(g[i, , drop = FALSE])) == 2))
  ## no product is a whole-genome copy of a parent, and only a copy would keep a parent's lineage
  whole <- rowSums(g) %in% c(0, ncol(g))
  expect_equal(sum(whole), 0)
  expect_true(all(lineage[!whole] > 0))
  ## meiosis I: in each oocyst the products form two sister pairs, each pair carrying the same
  ## parent at all 14 centromeres, the two pairs opposite
  key <- apply(g[, at_cen], 1, paste, collapse = "")
  expect_true(all(vapply(tetrads, function(i) identical(sort(as.vector(table(key[i]))), c(2L, 2L)), logical(1))))
  sisters <- lapply(tetrads, function(i) split(i, key[i]))
  ## any two products share a third of the genome on average and never nothing; sisters share
  ## more and the others less, as P(same parent) = 1/3 + 2/3 e^(-3rt) or 1/3 - 1/3 e^(-3rt) at a
  ## map distance rt from the centromere, averaged over the probes off the centromeres
  off <- setdiff(seq_along(probes), at_cen)
  share <- function(a, b) mean(g[a, off] == g[b, off])
  pairs <- do.call(rbind, lapply(tetrads, function(i) t(utils::combn(i, 2))))
  sib <- apply(pairs, 1, function(ab) share(ab[1], ab[2]))
  is_sister <- apply(pairs, 1, function(ab) any(vapply(sisters[[as.character(oocyst[ab[1]])]], setequal, logical(1), ab)))
  expect_gt(min(sib), 0)
  expect_equal(mean(sib), 1 / 3, tolerance = 1e-9)
  e <- exp(-3 * r * abs(probes[off] - cens[chr[off]]))
  expect_equal(sum(is_sister), 800)
  expect_equal(mean(sib[is_sister]), mean(1 / 3 + 2 / 3 * e), tolerance = 0.03)
  expect_equal(mean(sib[!is_sister]), mean(1 / 3 - 1 / 3 * e), tolerance = 0.03)
  ## each product crosses over at the gamete rate, a Poisson process: switches between
  ## neighbouring probes of one chromosome have the Haldane mean over the probe intervals, and a
  ## variance about the mean
  within <- which(diff(chr) == 0)
  sw <- rowSums(g[, within] != g[, within + 1])
  expect_equal(mean(sw), sum((1 - exp(-2 * r * diff(probes)[within])) / 2), tolerance = 0.04)
  expect_gt(var(sw) / mean(sw), 0.8)
  expect_lt(var(sw) / mean(sw), 1.1)
  ## chromosomes assort independently: a product changes parent at a junction with probability
  ## 1/2, while neighbouring probes within a chromosome are linked
  junction <- which(diff(chr) == 1)
  expect_equal(mean(g[, junction] != g[, junction + 1]), 1 / 2, tolerance = 0.06)
  ibd <- g[pairs[, 1], ] == g[pairs[, 2], ]
  within_phi <- sapply(1:14, function(ch) { k <- which(chr == ch)[10]; stats::cor(as.numeric(ibd[, k]), as.numeric(ibd[, k + 1])) })
  expect_gt(mean(within_phi), 0.5)
  ## a crossover's two products are one from each pole. On a chromosome where one crossover shows
  ## (two products switch once, between the same probes), they are never the first two products
  ## made or the last two (seedTetrad()'s 0 and 1, 2 and 3: the meiosis II sisters). The pairing
  ## they make is drawn for each chromosome, so two such chromosomes of one oocyst share it half
  ## the time; a random deal of chromatids would give a third, one pairing for the genome always.
  one <- do.call(rbind, lapply(tetrads, function(i) do.call(rbind, lapply(1:14, function(ch) {
    k <- which(chr == ch); s <- g[i, k[-1], drop = FALSE] != g[i, k[-length(k)], drop = FALSE]
    who <- which(rowSums(s) > 0)
    if (length(who) != 2 || any(rowSums(s)[who] != 1) || which(s[who[1], ]) != which(s[who[2], ])) return(NULL)
    data.frame(oocyst = oocyst[i[1]], pair = paste(who, collapse = ""))
  }))))
  expect_gt(nrow(one), 500)
  expect_false(any(one$pair %in% c("12", "34")))
  pairing <- ifelse(one$pair %in% c("13", "24"), "13|24", ifelse(one$pair %in% c("14", "23"), "14|23", "12|34"))
  same <- unlist(lapply(split(pairing, one$oocyst), function(p) if (length(p) < 2) NULL else {
    ij <- utils::combn(length(p), 2); p[ij[1, ]] == p[ij[2, ]] }))
  expect_equal(mean(same), 1 / 2, tolerance = 0.15)
})

test_that("infectious_bites counts each host's bites from a mosquito past the EIP", {
  skip_if_no_slim()
  ## every infectious bite infects (P_TRANS_VH = 1) and no host is ever full, so each one is a
  ## local inoculation, taken or not, also on a host that relapsed or was bitten earlier that day;
  ## case management makes such hosts common
  run <- run_sim(vivax_params(N_HOSTS = 40, N_MOSQ = 100, SIM_DAYS = 150, BITING_RATE = 0.2,
                              P_TRANS_VH = 1, MAX_MOI = 1000L, BLOODCURE_COVERAGE = 0.8,
                              CHROMOSOMES = "1"),
                 seed = 4, outdir = tempfile("mfw-eir"))
  st <- load_state(run, entity = "all")
  ev <- load_events(run)
  expect_gt(nrow(st[entity_type == "mosquito"]), 0)
  expect_true(all(is.na(st[entity_type == "mosquito", infectious_bites])))
  bites <- st[entity_type == "human", .(tick, human_id = entity_id, infectious_bites)]
  inoc <- ev[event_type %in% c("vh_inoculation", "vh_inoculation_no_take"),
             .(n = data.table::uniqueN(inoculation_id)), by = .(tick, human_id)]
  m <- merge(bites, inoc, by = c("tick", "human_id"), all.x = TRUE)
  m[is.na(n), n := 0L]
  expect_gt(sum(m$infectious_bites), 100)
  expect_equal(m$infectious_bites, m$n)
})

test_that("MAX_MOI holds within a day, counting the clones a host gained earlier that day", {
  skip_if_no_slim()
  ## no relapses, which bypass the check, and one genotype per bite, since the check admits a bite
  ## whole; bites that always infect
  run <- run_sim(vivax_params(N_HOSTS = 30, N_MOSQ = 150, SIM_DAYS = 200, BITING_RATE = 0.3,
                              P_TRANS_VH = 1, MAX_MOI = 2L, NU_HYPNOZOITE_BATCH = 0,
                              LAMBDA_SPOROZOITES = 0, CHROMOSOMES = "1"),
                 seed = 1, outdir = tempfile("mfw-maxmoi"))
  st <- load_state(run)[entity_type == "human"]
  expect_gt(sum(st$moi_active > 1), 100)
  expect_lte(max(st$moi_active + st$moi_subpatent), 2L)
})

test_that("ibd_truth_pairs() gives typed rows, or none, whatever pairs a run holds", {
  skip_if_no_slim(); skip_if_no_tskit()
  founders <- function(k) suppressMessages(run_sim(
    vivax_params(N_HOSTS = 20, N_MOSQ = 10, BITING_RATE = 1e-6, SIM_DAYS = 10, N_FOUNDER_HOSTS = k,
                 N_FOUNDER_LINEAGES = k, CHROMOSOMES = "1"), seed = 1, outdir = tempfile("mfwd_founders_")))
  types <- function(d) vapply(d, function(x) class(x)[1], "")
  want <- c(cls = "factor", a = "integer", b = "integer", ibd = "numeric")
  z <- ibd_truth_pairs(founders(3))
  expect_identical(nrow(z), 0L)
  expect_identical(types(z), want)
  six <- ibd_truth_pairs(founders(6))
  expect_identical(types(six), want)
  expect_identical(as.character(unique(six$cls)), "distinct founders")
  expect_identical(nrow(six), 15L)
  expect_true(all(six$ibd == 0))
  expect_no_error(ibd_truth_pairs(w_run(2)))
})
