# Clade-focused nowcast plots: chosen clades (all others summed into
# "other"), selected horizons, target dates on the x-axis.
#
# Set everything in "Settings" and run top to bottom. Plot settings can be
# changed by re-running just "Plot one state" and "Save"; re-run "Load" and
# "Clades" when nowcast_date, locations, models, or clades change.

library(here)
source(here("results", "nowcast_sample_densities", "plot_nowcast_sample_densities.R"))

# ---- Settings ----
hub_path     <- "~/GitHub/variant-nowcast-hub"
nowcast_date <- "2024-12-25"
locations    <- c("MA")                          # or c("CA", "IL", "MA"), or "all"
models       <- c("UMass-HMLR", "Hub-baseline")  # 1 or 2
clades       <- NULL                             # e.g. c("25B", "25C"); NULL = top two for the round

# Horizons: days relative to the nowcast date (hub range -31 to 10).
# The plot shows the corresponding target dates.
horizons     <- -31:-28                             # e.g. -3:0 , -31:10, seq(-28, 7, by = 7), c(-7, 0, 7)

model_colors <- c("#881C1C", "#0ba7e0")          # same order as models
alpha        <- 0.4                              # density fill opacity
free_y       <- FALSE                            # TRUE = per-clade y scales
trim         <- TRUE                             # FALSE = densities extend past the samples
rug          <- FALSE                            # TRUE = tick per sample
date_axis    <- TRUE                             # FALSE = days relative to nowcast date
slab_scale   <- 0.9                              # share of the gap between dates each density may fill
page_width   <- 11                               # PDF page size, inches
page_height  <- 5

stopifnot(all(horizons >= -31 & horizons <= 10))
target_dates <- as.Date(nowcast_date) + horizons

# ---- Load ----
samples   <- load_hub_samples(hub_path = hub_path, nowcast_date = nowcast_date,
                              location = locations, models = models)
eval_data <- load_eval_data(hub_path = hub_path, nowcast_date = nowcast_date,
                            location = locations, final_days = 90)  # NULL if not final
if (identical(locations, "all")) locations <- sort(unique(samples$location))

# ---- Clades (skip this block to plot all clades) ----
if (is.null(clades)) {
  clades <- default_clades(hub_path = hub_path, nowcast_date = nowcast_date,
                           samples = samples, n = 2)
}
samples <- collapse_clades(df = samples, clades = clades, value_col = "value")
if (!is.null(eval_data)) {
  eval_data <- collapse_clades(df = eval_data, clades = clades, value_col = "oracle_value")
}

# Plot settings shared by the single-state preview and the saved PDF.
plot_args <- list(
  horizons     = horizons,
  model_colors = model_colors,
  alpha        = alpha,
  free_y       = free_y,
  trim         = trim,
  rug          = rug,
  date_axis    = date_axis,
  slab_scale   = slab_scale
)

# ---- Plot one state ----
loc <- locations[1]
p <- do.call(plot_sample_densities, c(
  list(samples   = filter(samples, location == loc),
       eval_data = if (is.null(eval_data)) NULL else filter(eval_data, location == loc)),
  plot_args
))
p

# ---- Save every state, one page each ----
loc_label <- if (length(locations) > 5) {
  paste0(length(locations), "locs")
} else {
  paste(locations, collapse = "-")
}
out_file <- here("results", "nowcast_sample_densities", "output", sprintf(
  "%s_%s_%s_%s_%s-to-%s.pdf", nowcast_date, loc_label,
  paste(models, collapse = "_vs_"), paste(clades, collapse = "-"),
  format(min(target_dates), "%m%d"), format(max(target_dates), "%m%d")
))
do.call(save_state_pages, c(
  list(samples   = samples,
       file      = out_file,
       locations = locations,
       eval_data = eval_data,
       width     = page_width,
       height    = page_height),
  plot_args
))
