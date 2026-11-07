## Runs shared by the window and episode tests, made once per session.

## W goes extinct early for seeds 1 to 6 (days 43 to 196), before its burn-in of 200 days ends.
w_params <- function() vivax_params(N_HOSTS = 60, N_MOSQ = 240, BITING_RATE = 0.12, SIM_DAYS = 400,
                                    N_FOUNDER_HOSTS = 12, NU_HYPNOZOITE_BATCH = 0,
                                    MEAN_BLOOD_STAGE_DURATION = 10, CHROMOSOMES = "1")
.w_runs <- new.env()
w_run <- function(seed) {
  key <- as.character(seed)
  if (is.null(.w_runs[[key]]))
    .w_runs[[key]] <- suppressMessages(run_sim(w_params(), seed = seed, outdir = tempfile("mfwd_w_")))
  .w_runs[[key]]
}
