## Regenerates the figures in man/figures/ that README.md and MODEL_DESCRIPTION.md embed.
##
## These used to be made ad hoc, so when one went stale there was no record of the parameters
## behind it. Everything here is seeded, so running it twice gives the same figure.
##   Rscript dev/make-readme-figures.R
suppressMessages({library(data.table); library(ggplot2); devtools::load_all(".", quiet = TRUE)})

FIG  <- "man/figures"
SEED <- 7L
OUT  <- file.path(tempdir(), "mfw-figs")

## Two years, at the same transmission setting as the representative run in
## MODEL_DESCRIPTION.md, so the figures in the two documents describe one model and not two.
p <- vivax_params(
  N_HOSTS = 200L, N_MOSQ = 400L, SIM_DAYS = 730L, BITING_RATE = 0.08,
  N_FOUNDER_HOSTS = 40L, N_FOUNDER_LINEAGES = 10L, NU_HYPNOZOITE_BATCH = 3,
  P_DETECT_PRIMARY = 0.9, P_DETECT_RELAPSE = 0.15
)
unlink(OUT, recursive = TRUE)
run <- run_sim(p, seed = SEED, outdir = OUT)
ev  <- as.data.table(load_events(run))
st  <- as.data.table(load_state(run))

cat(sprintf("run: %d hosts, %d days, seed %d\n", p$N_HOSTS, p$SIM_DAYS, SEED))
print(summarise_run(run))

## Choose the host DELIBERATELY, against what the figure is meant to show: several distinct
## lineages (the panel stacks them one colour each, so a single-lineage host is just one block),
## several relapses, and an early first infection so the two years are actually used.
## pick_host() optimises for relapse count alone and on this run returns a host with ONE lineage.
onsets <- ev[event_type %in% .onset_events() & is_dormant == 0,
             .(human_id, lineage_id, tick, event_type)]
cand <- onsets[, .(n_lin = uniqueN(lineage_id),
                   n_rel = sum(event_type == "hypnozoite_activation"),
                   first = min(tick)), by = human_id]
cand <- cand[n_lin >= 3 & n_rel >= 3][order(first, -n_lin)]
stopifnot(nrow(cand) > 0)
h <- cand$human_id[1]
cat(sprintf("host %d: %d distinct lineages, %d relapses, first infection day %d\n",
            h, cand$n_lin[1], cand$n_rel[1], cand$first[1]))

g <- plot_host_trajectory(st, ev, host_id = h, layout = "lineages")
ggsave(file.path(FIG, "README-host-trajectory.png"), g,
       width = 9, height = 4.6, dpi = 150, bg = "white")
cat("wrote", file.path(FIG, "README-host-trajectory.png"), "\n")
