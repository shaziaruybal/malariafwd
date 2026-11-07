# Declare the data.table column names used unquoted inside [ ] expressions, so
# R CMD check does not flag them as "no visible binding for global variable".
# This is the standard idiom for data.table-based packages.
utils::globalVariables(c(
  ".", "N",
  # state / entity columns
  "entity_type", "entity_id", "tick", "moi_active", "moi_dormant", "hyp_reservoir_size",
  # event-log columns
  "event_type", "is_dormant", "is_imported", "lineage_id", "oocyst_id", "pedigree_id",
  "parent1_pedigree_id", "parent2_pedigree_id", "human_id", "origin", "detected",
  # derived / aggregate columns
  "pct", "cls", "frac", "win", "secondary", "source", "target", "chain_size", "node", "comp",
  "mean_ibd", "n_pairs",
  # IBD
  "node_a", "node_b", "ibd_fraction", "node_id", "a", "b", "ibd",
  "left", "right", "mrca_node", "mrca_time",
  # time-to-relapse / dormancy
  "t0", "tr", "time", "event", "act_tick", "est_tick", "nrel", "py",
  # moi / within-host trajectory
  "peak_moi", "n_relapse", "kind",
  # classify_who
  "imp_tick", "recent", "recent_travel", "surveillance",
  # ibd_relapse_origin
  "rp", "rt", "pp", "pt",
  # Kaplan-Meier
  "d_i", "c_i", "n_i", "surv",
  # transmission_links
  "A1", "A2", "C",
  # plot_diagnostics / force_of_infection
  "w", "y", "rate", "n", "days", "series", "is_acquisition", "recrud_tick", "inf", "hyp",
  # plot_moi_vs_relapse / plot_relatedness_ladder
  "grp", "band",
  # classify_cases imported-origin and masked-relapse flags
  "imported", "into_patent_lineage",
  # segments / ladder
  "gstart", "gend", "segid", "cs", "ce", "gs", "ge", "chrom", "x0", "x1", "lab", "cls", "ped",
  # data.table columns that shadow stats::start / stats::end
  "start", "end", "pair", "len_mb", "ypos", "lineage", "r", "V1"
))
