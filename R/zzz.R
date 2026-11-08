# Declare the data.table column names used unquoted inside [ ] expressions, so
# R CMD check does not flag them as "no visible binding for global variable".
# This is the standard idiom for data.table-based packages.
utils::globalVariables(c(
  ".", "N",
  # state / entity columns
  "entity_type", "entity_id", "tick", "moi_active", "moi_dormant", "hyp_reservoir_size",
  # event-log columns
  "event_type", "is_dormant", "is_imported", "lineage_id", "oocyst_id", "pedigree_id", "inoculation_id",
  "parent1_pedigree_id", "parent2_pedigree_id", "human_id", "origin", "detected",
  # derived / aggregate columns
  "cls", "frac", "win", "chain_id", "chain_size", "chain", "hypnozoites_per_host",
  "mean_ibd", "n_pairs",
  # IBD
  "node_a", "node_b", "ibd_fraction", "node_id", "a", "b", "ibd",
  "left", "right", "mrca_node", "mrca_time", "mrca_id", "node_time", "depth", "tract",
  # time-to-relapse / dormancy
  "t0", "tr", "time", "event", "nrel", "py", "batch_pedigree_id", "batch_tick",
  "n_hypnozoites", "masked", "relapse", "t_relapse", "t_empty",
  # within-host trajectory
  "kind",
  # classify_travel
  "imp_tick", "recent", "recent_travel", "surveillance",
  # Kaplan-Meier
  "d_i", "c_i", "n_i", "surv",
  # simulate_genotypes / write_vcf
  "root_id", "locus", "pos", "position", "plaf", ".row", "INFO", "FORMAT",
  # classify_chain
  "from_import", "generation", "import_pedigree_id", "who_class", "src", "to",
  # rc_by_origin
  "n_cases", "n_secondary", "n_onward", "Rc", "onward_frac",
  # transmission_links
  "source_pedigree_id", "target_pedigree_id", "source_human_id", "target_human_id", "mosq_id",
  "uptake_tick", "inoculation_tick", "target_tick", "source_tick", "weight", "p1", "p2", "batch",
  # plot_diagnostics / force_of_infection
  "w", "y", "rate", "n", "days", "series", "is_acquisition", "recrud_tick", "inf", "hyp",
  # plot_relatedness_ladder
  "band",
  # classify_cases imported-origin and masked-relapse flags
  "imported", "into_patent_lineage",
  # episodes / recurrences
  "episode_id", "onset_tick", "onset_event", "end_tick", "end_reason", "prev_tick", "prior_episode_id",
  "prior_tick", "prior_pedigree_id", "prior_lineage_id", "prior_origin", "prior_imported", "prior_detected",
  "days_since", "x.prior_tick", "i.tick", "is_imported",
  # segments / ladder
  "gstart", "gend", "segid", "cs", "ce", "gs", "ge", "chrom", "x0", "x1", "lab", "cls", "ped",
  # data.table columns that shadow stats::start / stats::end
  "start", "end", "pair", "len_mb", "ypos", "lineage", "r", "V1"
))
