## Model invariants that guard against regressions. These need SLiM, so they skip when it is
## unavailable (e.g. on CRAN / CI without the engine) via skip_if_no_slim() in helper-slim.R.

test_that("hypnozoite batches never draw more than they hold, and ids are unique", {
  skip_if_no_slim()
  out <- file.path(tempdir(), "mfw-invariants")
  unlink(out, recursive = TRUE)
  p <- vivax_params(N_HOSTS = 60, N_MOSQ = 120, BITING_RATE = 0.12, SIM_DAYS = 300,
                    N_FOUNDER_LINEAGES = 8)
  run <- run_sim(p, seed = 1, outdir = out)
  ev <- data.table::as.data.table(load_events(run))

  ## (WS-G bug 1) a batch's activations + deaths must never exceed the size it was
  ## established with. The old PENDING_HYP_DEATHS queue was scoped local, so retired
  ## batches stayed alive and kept activating (64% of batches over-drew).
  est <- ev[grepl("^vh_inoculation", event_type) & is_dormant == 1,
            .(batch = pedigree_id, held = n_hypnozoites)]
  act <- ev[event_type == "hypnozoite_activation", .(batch = parent1_pedigree_id, n = n_hypnozoites)]
  die <- ev[grepl("hypnozoite_clear", event_type), .(batch = pedigree_id, n = n_hypnozoites)]
  drawn <- rbind(act, die)[, .(drew = sum(n)), by = batch]
  chk <- merge(est, drawn, by = "batch", all.x = TRUE)[is.na(drew), drew := 0]
  expect_equal(sum(chk$drew > chk$held), 0)

  ## (WS-G bugs 2 & 3) event ids unique; inoculation ids do not restart each tick
  expect_equal(data.table::uniqueN(ev$event_id), nrow(ev))
  inoc <- ev[inoculation_id > 0, inoculation_id]
  expect_gt(max(inoc), 0.5 * data.table::uniqueN(inoc))
})

test_that("classify_cases labels imported relapses", {
  skip_if_no_slim()
  out <- file.path(tempdir(), "mfw-imported")
  unlink(out, recursive = TRUE)
  p <- vivax_params(N_HOSTS = 60, N_MOSQ = 120, BITING_RATE = 0.06, SIM_DAYS = 400,
                    IMPORT_START = 1L, IMPORT_END = 99999L, IMPORT_RATE_PER_TICK = 0.15)
  run <- run_sim(p, seed = 1, outdir = out)
  cs <- classify_cases(load_events(run))
  expect_true("imported" %in% names(cs))
  ## a relapse of an imported genotype must be labelled imported (was hard-coded FALSE)
  expect_gt(nrow(cs[origin == "relapse" & imported == TRUE]), 0)
})

test_that("deprecated N_OOCYSTS_PER_BITE is rejected", {
  expect_error(vivax_params(N_OOCYSTS_PER_BITE = 2), "deprecated")
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
  skip_if_no_slim()
  skip_if_not(requireNamespace("reticulate", quietly = TRUE), "reticulate not available")
  out <- file.path(tempdir(), "mfw-arg")
  unlink(out, recursive = TRUE)
  run <- run_sim(vivax_params(N_HOSTS = 60, N_MOSQ = 120, BITING_RATE = 0.12,
                              SIM_DAYS = 300, N_FOUNDER_LINEAGES = 8), seed = 1, outdir = out)
  ev <- data.table::as.data.table(load_events(run))
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

  ## D22: tmrca is a DEPTH, so one fixed relationship gives one number wherever it was
  ## sampled. Before the fix this returned the ancestor's raw node time and the same
  ## relapse-vs-origin relationship spanned 43-391 ticks on a 400-tick run.
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

test_that("ibd_pairs returns all ancestry, and segments reconstruct it", {
  skip_if_no_slim()
  skip_if_not(requireNamespace("reticulate", quietly = TRUE), "reticulate not available")
  out <- file.path(tempdir(), "mfw-ibdseg")
  unlink(out, recursive = TRUE)
  run <- run_sim(vivax_params(N_HOSTS = 50, N_MOSQ = 100, BITING_RATE = 0.12,
                              SIM_DAYS = 250, N_FOUNDER_LINEAGES = 8), seed = 3, outdir = out)
  ev  <- data.table::as.data.table(load_events(run))
  ids <- unique(data.table::as.data.table(classify_cases(ev))$pedigree_id)
  skip_if(length(ids) < 10, "too few cases in this run")

  ## D22: no time cutoff is offered, so a clonal relapse reads 1 however old its ancestor is
  idx <- ibd_index(ibd_pairs(run, ids = ids))
  ro  <- data.table::as.data.table(ibd_relapse_origin(ev, idx))
  if (nrow(ro)) expect_true(all(ro$ibd > 0.999))

  sg <- ibd_pairs(run, ids = ids, segments = TRUE)
  expect_named(sg$segments,
               c("node_a", "node_b", "left", "right", "mrca_node", "mrca_time"))
  expect_gt(sg$sequence_length, 0)
  ## the segments must add back up to the pair table they came from
  rec <- sg$segments[, .(f = sum(right - left) / sg$sequence_length), by = .(node_a, node_b)]
  cmp <- merge(data.table::as.data.table(sg$ibd), rec, by = c("node_a", "node_b"))
  expect_equal(nrow(cmp), nrow(sg$ibd))
  expect_lt(max(abs(cmp$ibd_fraction - cmp$f)), 1e-9)
  ## mrca_time is a node time, never more recent than either sample it joins
  expect_true(all(sg$segments$mrca_time >= 0))
})

test_that("an inoculation's primary and its hypnozoite batch are the same genome", {
  skip_if_no_slim()
  ## Regression for D17. The import path used to build the batch with addEmpty(NULL), giving
  ## two genetically unrelated founders (measured IBD 0) under one lineage_id, while an
  ## ordinary bite gave IBD 1. Imports are now seeded in early() and the second compartment is
  ## cloned one tick later, so both paths must give IBD 1.
  out <- file.path(tempdir(), "mfw-import-clone")
  unlink(out, recursive = TRUE)
  ## imports on, so both the founder path and the ordinary path are exercised
  p <- vivax_params(N_HOSTS = 80, N_MOSQ = 160, BITING_RATE = 0.10, SIM_DAYS = 500,
                    N_FOUNDER_LINEAGES = 8, NU_HYPNOZOITE_BATCH = 5,
                    IMPORT_START = 1L, IMPORT_END = 500L, IMPORT_RATE_PER_TICK = 0.08)
  run <- run_sim(p, seed = 3, outdir = out)
  ev <- data.table::as.data.table(load_events(run))

  ## One inoculation delivers one sporozoite population, so the blood-stage infection and the
  ## hypnozoite batch established alongside it must be clones. The import path used to build
  ## the batch with addEmpty(NULL), giving two unrelated founders (IBD 0) under one lineage_id
  ## while a normal bite gave IBD 1 -- undetectable from the event log alone, which is why
  ## this check goes through the tree sequence.
  ## grouped on inoculation_id, NOT (tick, human_id): under D17 an import's two compartments
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
  idx  <- ibd_index(ibd_pairs(run, ids = unique(c(take$primary, take$batch))))
  take[, ibd := idx$getibd(primary, batch)]
  expect_equal(take$ibd, rep(1, nrow(take)), tolerance = 1e-8)

  ## and the import path specifically, since that is where the defect was
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
    "  paste(asInteger(match(PROBES, ind.haplosomes.mutations.position) >= 0), sep = '')); sim.simulationFinished(); }"), script)
  res <- processx::run(slim_path(), script, error_on_status = FALSE)
  expect_equal(res$status, 0L, info = res$stderr)
  L <- grep("^P ", strsplit(res$stdout, "\n")[[1]], value = TRUE)
  x <- do.call(rbind, strsplit(L, " "))
  oocyst <- as.integer(x[, 2]); lineage <- as.integer(x[, 3])
  g <- do.call(rbind, lapply(strsplit(x[, 4], ""), as.integer))     # 1 = parent1's allele
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
})

test_that("infectious_bites counts each host's bites from a mosquito past the EIP", {
  skip_if_no_slim()
  ## every infectious bite infects (P_TRANS_VH = 1), no host is ever full and nothing relapses, so
  ## each one is a local inoculation, taken or not
  run <- run_sim(vivax_params(N_HOSTS = 40, N_MOSQ = 100, SIM_DAYS = 120, BITING_RATE = 0.2,
                              P_TRANS_VH = 1, MAX_MOI = 1000L, NU_HYPNOZOITE_BATCH = 0,
                              CHROMOSOMES = "1"),
                 seed = 4, outdir = file.path(tempdir(), "mfw-eir"))
  st <- data.table::as.data.table(load_state(run)); ev <- data.table::as.data.table(load_events(run))
  expect_true(all(is.na(st[entity_type == "mosquito", infectious_bites])))
  bites <- st[entity_type == "human", .(tick, human_id = entity_id, infectious_bites)]
  inoc <- ev[event_type %in% c("vh_inoculation", "vh_inoculation_no_take"),
             .(n = data.table::uniqueN(inoculation_id)), by = .(tick, human_id)]
  m <- merge(bites, inoc, by = c("tick", "human_id"), all.x = TRUE)
  m[is.na(n), n := 0L]
  expect_gt(sum(m$infectious_bites), 100)
  expect_true(all(m$infectious_bites >= m$n))
  expect_true(all(m[n == 0, infectious_bites] == 0))
  ## except a further bite on a host whose only clone arrived earlier that day, which is dropped:
  ## the model does not see that clone until the day's bites are over (D47)
  expect_gt(sum(m$n) / sum(m$infectious_bites), 0.98)
})
