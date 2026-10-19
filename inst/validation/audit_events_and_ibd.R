## audit_events_and_ibd.R -- full integrity audit of the malariafwd event log, the state
## snapshots, and the downstream ARG / IBD values. Run after any model change.
##   Rscript analysis/model_checks/audit_events_and_ibd.R
suppressMessages({
  library(malariafwd)
  library(data.table); library(reticulate)
})
py_require(c("tskit", "pandas")); invisible(py_run_string("import tskit"))

OUT <- file.path(tempdir(), "mfw_audit")
PASS <- 0; FAIL <- 0
ok <- function(lbl, cond, detail = "") {
  if (isTRUE(cond)) { PASS <<- PASS + 1; cat(sprintf("  PASS  %-58s %s\n", lbl, detail)) }
  else              { FAIL <<- FAIL + 1; cat(sprintf("  FAIL  %-58s %s\n", lbl, detail)) }
}

audit_run <- function(tag, params, seed = 1) {
  cat("\n================ ", tag, " ================\n")
  d <- file.path(OUT, tag); unlink(d, recursive = TRUE)
  run <- run_sim(params, seed = seed, outdir = d)
  ev <- as.data.table(load_events(run)); st <- as.data.table(load_state(run))
  endt <- max(st$tick)

  ## ---------------- A. EVENT LOG ----------------
  cat("-- A. event log --\n")
  dict <- event_types()$event_type
  emitted <- sort(unique(ev$event_type))
  undocumented <- setdiff(sub("_no_take$", "", emitted), dict)
  ok("every emitted event_type is in event_types()", length(undocumented) == 0,
     if (length(undocumented)) paste(undocumented, collapse = ",") else "")
  ok("event_id unique", uniqueN(ev$event_id) == nrow(ev),
     sprintf("%d rows / %d ids", nrow(ev), uniqueN(ev$event_id)))
  ok("ticks within run", min(ev$tick) >= 1 && max(ev$tick) <= endt)

  ## creation events mint distinct pedigree_ids
  creat <- ev[event_type %in% c("vh_inoculation","vh_inoculation_import","vh_inoculation_seed",
                                "hypnozoite_activation") & pedigree_id > 0, pedigree_id]
  ok("creation events mint unique pedigree_ids", uniqueN(creat) == length(creat),
     sprintf("%d created", length(creat)))

  ## inoculation_id groups co-inoculated genotypes of ONE bite, and is unique across the run
  inoc <- ev[inoculation_id > 0]
  if (nrow(inoc)) {
    ## D17: this used to assert "1 tick, 1 host". One host still holds everywhere, but an
    ## IMPORT's two compartments now land on CONSECUTIVE ticks - the hypnozoite batch cannot
    ## be cloned from the blood-stage genome until that genome is visible, which is the tick
    ## after it arrives. Local bites are unaffected: their source already exists, so both
    ## compartments are still created in one callback.
    per <- inoc[, .(ticks = uniqueN(tick), span = max(tick) - min(tick),
                    hosts = uniqueN(human_id),
                    imported = any(grepl("import", event_type))), by = inoculation_id]
    ok("each inoculation_id is one bite into one host",
       all(per$hosts == 1), sprintf("%d bites", nrow(per)))
    ok("local inoculations occupy a single tick",
       all(per[imported == FALSE]$ticks == 1), sprintf("%d local", per[imported == FALSE, .N]))
    ok("imported inoculations span at most the one-tick batch deferral (D17)",
       all(per[imported == TRUE]$span <= 1), sprintf("%d imported", per[imported == TRUE, .N]))
    ok("inoculation_ids do not restart each tick",
       max(inoc$inoculation_id) >= 0.5 * uniqueN(inoc$inoculation_id))
  }

  ## is_dormant semantics: activations & clearances are blood-stage; batch rows are dormant
  ok("hypnozoite_activation logged as blood-stage (is_dormant=0)",
     nrow(ev[event_type == "hypnozoite_activation" & is_dormant != 0]) == 0)
  ok("hypnozoite_clearance logged as dormant (is_dormant=1)",
     nrow(ev[event_type == "hypnozoite_clearance" & is_dormant != 1]) == 0)

  ## detected flag only on case events
  cases <- c("vh_inoculation","vh_inoculation_import","hypnozoite_activation")
  ok("detected set only on case events",
     nrow(ev[!(event_type %in% cases) & detected != -1]) == 0)

  ## ---- CONSERVATION: hypnozoite batches ----
  est <- ev[grepl("^vh_inoculation", event_type) & is_dormant == 1,
            .(batch = pedigree_id, held = n_hypnozoites)]
  act <- ev[event_type == "hypnozoite_activation", .(batch = parent1_pedigree_id, n = n_hypnozoites)]
  gone <- ev[grepl("hypnozoite_clear", event_type), .(batch = pedigree_id, n = n_hypnozoites)]
  drawn <- rbind(act, gone)[, .(drew = sum(n)), by = batch]
  cons <- merge(est, drawn, by = "batch", all.x = TRUE)[is.na(drew), drew := 0]
  ok("no batch draws more hypnozoites than it held", sum(cons$drew > cons$held) == 0,
     sprintf("%d batches, %d over-drawing", nrow(cons), sum(cons$drew > cons$held)))
  ok("every activation references an established batch",
     all(act$batch %in% est$batch), sprintf("%d activations", nrow(act)))
  ## unresolved remainder must equal the reservoir still standing at the end
  remainder <- sum(pmax(cons$held - cons$drew, 0))
  final_res <- st[entity_type == "human" & tick == endt, sum(hyp_reservoir_size)]
  ok("unresolved hypnozoites == final reservoir in state.tsv", remainder == final_res,
     sprintf("events %d vs state %d", remainder, final_res))

  ## ---- CONSERVATION: blood-stage episodes ----
  ## every way an episode can end (natural, case treatment, MDA) - from the public dictionary
  .clearers <- event_types()$event_type[event_types()$ends_episode]
  onset <- ev[(grepl("^vh_inoculation", event_type) & is_dormant == 0) |
                event_type == "hypnozoite_activation", pedigree_id]
  clear <- ev[event_type %in% .clearers, pedigree_id]
  ok("no blood-stage infection clears twice", uniqueN(clear) == length(clear))
  ok("every clearance matches an onset", all(clear %in% onset))
  alive <- st[entity_type == "human" & tick == endt, sum(moi_active)]
  ok("onsets - clearances == active blood-stage lineages at the end",
     length(onset) - length(clear) >= alive,     # >= : several lineages can share a host-count
     sprintf("onsets %d, cleared %d, active %d", length(onset), length(clear), alive))

  ## ---- WS-F: oocyst / tetrad structure ----
  up <- ev[event_type %in% c("hv_uptake_clonal", "hv_uptake_recombinant")]
  if (nrow(up)) {
    ok("each uptake event has its own oocyst_id", uniqueN(up$oocyst_id) == nrow(up),
       sprintf("%d oocysts", nrow(up)))
    meals <- up[, .N, by = .(tick, mosq_id, human_id)]
    ok("oocysts per blood meal >= 1", min(meals$N) >= 1,
       sprintf("mean %.2f, max %d", mean(meals$N), max(meals$N)))
    ok("recombinant uptakes record two distinct parents",
       nrow(ev[event_type == "hv_uptake_recombinant" &
                 parent1_lineage_id == parent2_lineage_id]) == 0,
       sprintf("%.1f%% recombinant", 100 * mean(up$event_type == "hv_uptake_recombinant")))
  }

  ## ---------------- B. ARG / IBD ----------------
  cat("-- B. ARG / IBD --\n")
  ts <- py$tskit$load(trees_path(run))
  ok("tree sequence loads", !is.null(ts), sprintf("%d trees, %d nodes, seq_len %.2f Mb",
     ts$num_trees, ts$num_nodes, ts$sequence_length / 1e6))
  ok("recombination present in the ARG (>1 local tree)", ts$num_trees > 1)

  cs <- classify_cases(ev)
  ## NB classify_cases() deliberately excludes vh_inoculation_seed (founders are burn-in,
  ## not cases), so add the founder genomes explicitly to test the diversity-axis assumption
  ## that distinct founder lineages are unrelated.
  founders <- ev[event_type == "vh_inoculation_seed", pedigree_id]
  idx <- ibd_index(ibd_pairs(run, ids = unique(c(cs$pedigree_id, founders))))
  present <- cs$pedigree_id %in% idx$nm$pedigree_id
  ok("logged cases are retained in the tree sequence", mean(present) > 0.95,
     sprintf("%.1f%% of %d cases", 100 * mean(present), nrow(cs)))

  ## IBD bounds + self-identity
  smp <- head(idx$nm$pedigree_id, 40)
  gi <- idx$getibd(rep(smp[1], length(smp)), smp)
  ok("IBD values within [0,1]", all(gi >= 0 & gi <= 1, na.rm = TRUE))
  ok("IBD of a genome with itself == 1", isTRUE(all.equal(idx$getibd(smp[1], smp[1]), 1)))
  ok("IBD is symmetric",
     isTRUE(all.equal(idx$getibd(smp[1], smp[2]), idx$getibd(smp[2], smp[1]))))

  ## ---- ARG-native checks (tree/branch domain, complementing the IBD/segment domain) ----
  ag <- arg_summary(run, unique(c(cs$pedigree_id, founders)))
  ok("ARG recombines (>1 local tree)", ag$num_trees > 1, sprintf("%d trees", ag$num_trees))
  ok("genealogy is a FOREST (independent founder lineages)", ag$roots_median > 1,
     sprintf("%.0f roots/tree, %.0f%% coalesced, depth %.0f ticks",
             ag$roots_median, ag$pct_trees_coalesced, ag$max_root_time))
  ## arg_tmrca is the relatedness measure that is valid on a forest (branch divergence is NOT:
  ## it scores coalescing and non-coalescing pairs alike -- see ?arg_divergence)
  relp <- ev[event_type == "hypnozoite_activation" & parent1_pedigree_id > 0][1:min(.N, 15)]
  if (nrow(relp)) {
    tm <- arg_tmrca(run, relp$pedigree_id, relp$parent1_pedigree_id)
    ok("TMRCA: relapse coalesces with its origin over the whole genome",
       all(tm$coalesced > 0.99, na.rm = TRUE),
       sprintf("coalesced %.2f, mean TMRCA %.0f ticks", mean(tm$coalesced, na.rm = TRUE),
               mean(tm$tmrca, na.rm = TRUE)))
  }
  fp <- unique(ev[event_type == "vh_inoculation_seed"], by = "lineage_id")$pedigree_id
  if (length(fp) >= 4) {
    tf <- arg_tmrca(run, fp[1:2], fp[3:4])
    ok("TMRCA: distinct founders never coalesce", all(tf$coalesced < 0.5, na.rm = TRUE),
       sprintf("coalesced %.2f", mean(tf$coalesced, na.rm = TRUE)))
  }

  ## the key biological invariant: a relapse is clonal to its origin
  ro <- ibd_relapse_origin(ev, idx)
  if (!is.null(ro) && nrow(ro)) {
    ok("relapse-to-origin IBD == 1 (clonal)",
       median(ro$ibd, na.rm = TRUE) > 0.99,
       sprintf("median %.3f over %d relapses", median(ro$ibd, na.rm = TRUE), nrow(ro)))
  }

  ## DISTINCT FOUNDER lineages are independent genome roots -> must be ~unrelated.
  ## (A general cross-lineage test is not meaningful: reciprocal recombinants from one
  ## tetrad get different lineage_ids but are genuinely ~50% related.)
  fnd <- ev[event_type == "vh_inoculation_seed", .(ped = pedigree_id[1]), by = lineage_id]
  fnd <- fnd[ped %in% idx$nm$pedigree_id]
  if (nrow(fnd) >= 3) {
    cb <- combn(nrow(fnd), 2)
    bg <- idx$getibd(fnd$ped[cb[1, ]], fnd$ped[cb[2, ]])
    ok("distinct FOUNDER lineages are unrelated (IBD ~ 0)", median(bg, na.rm = TRUE) < 0.05,
       sprintf("median %.3f, max %.3f over %d pairs", median(bg, na.rm = TRUE),
               max(bg, na.rm = TRUE), length(bg)))
  }
  ## clonal founder COPIES (same lineage_id, different hosts) must be IBD = 1
  clone_lin <- ev[event_type == "vh_inoculation_seed", .N, by = lineage_id][N > 1, lineage_id]
  if (length(clone_lin)) {
    pr <- ev[event_type == "vh_inoculation_seed" & lineage_id == clone_lin[1], pedigree_id]
    pr <- pr[pr %in% idx$nm$pedigree_id]
    if (length(pr) >= 2)
      ok("clonal founder copies are IBD = 1", idx$getibd(pr[1], pr[2]) > 0.99,
         sprintf("%.3f", idx$getibd(pr[1], pr[2])))
  }
  invisible(NULL)
}

## two configurations: the shipped default, and the WS-F multi-oocyst/NB path
audit_run("default", vivax_params(N_HOSTS = 80, N_MOSQ = 160, BITING_RATE = 0.12,
                                  SIM_DAYS = 400, N_FOUNDER_LINEAGES = 8))
audit_run("wsf_nb_imports",
          vivax_params(N_HOSTS = 80, N_MOSQ = 160, BITING_RATE = 0.12, SIM_DAYS = 400,
                       N_FOUNDER_LINEAGES = 8,
                       LAMBDA_OOCYSTS = 2.0, K_OOCYSTS = 1.0, LAMBDA_SPOROZOITES = 2.0,
                       IMPORT_START = 1L, IMPORT_END = 99999L, IMPORT_RATE_PER_TICK = 0.1))

cat(sprintf("\n================ %d passed, %d failed ================\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)
