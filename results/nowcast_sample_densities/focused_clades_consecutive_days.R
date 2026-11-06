# Chosen clades over a few consecutive target dates, using
# plot_nowcast_sample_densities.R: all other clades are summed into "other",
# and target dates are on the x-axis.
#
# Edit the settings below and run the script top to bottom. After the first
# run, plot settings (colours, alpha, dates, free_y, trim, rug) can be
# changed by re-running just the "Plot" and "Save" sections; re-run "Load
# data" when nowcast_date, locations, models, or clades change.

# ---- Packages and functions ----
library(arrow)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggdist)
library(here)  # resolves paths from the repo root (scoring-multiclass-props)

source(here("results", "nowcast_sample_densities",
            "plot_nowcast_sample_densities.R"))

# ---- Settings to toggle ----
hub_path     <- "~/GitHub/variant-nowcast-hub"
nowcast_date <- "2026-05-06"                     # eval data needs a round 90+ days old (pull the hub)
locations    <- "CA"                             # one page per state; "all" = every submitted state
models       <- c("UMass-HMLR", "Hub-baseline")  # 1 or 2 models
model_colors <- c("#881C1C", "#0ba7e0")          # same order as models; hex or R names
clades       <- NULL                             # e.g. c("25B", "25C"); NULL = top two for the round
first_date   <- as.Date(nowcast_date) - 3        # first target date shown
n_days       <- 4                                # consecutive target dates shown
alpha        <- 0.4                              # density fill opacity
free_y       <- FALSE                            # TRUE = per-clade y scales instead of [0, 1]
trim         <- TRUE                             # FALSE = densities extend past the samples
rug          <- FALSE                            # TRUE = tick for each of the 100 samples
show_eval    <- TRUE                             # FALSE = skip observed / possible outcomes

# Target dates as days relative to the nowcast date (hub range: -31 to 10).
target_dates <- seq(first_date, by = 1, length.out = n_days)
horizons     <- as.integer(target_dates - as.Date(nowcast_date))

# ---- Load data ----
samples <- load_hub_samples(hub_path, nowcast_date, locations, models)
eval_data <- if (show_eval) {
  load_eval_data(hub_path, nowcast_date, locations)
} else {
  NULL
}
if (identical(locations, "all")) locations <- sort(unique(samples$location))

# Keep the chosen clades; sum all others into "other".
if (is.null(clades)) clades <- default_clades(hub_path, nowcast_date, samples)
samples <- collapse_clades(samples, clades, "value")
if (!is.null(eval_data)) {
  eval_data <- collapse_clades(eval_data, clades, "oracle_value")
}

# ---- Plot ----
# A quick look at the first state in the plot pane.
plot_args <- list(
  horizons = horizons,
  model_colors = model_colors,
  alpha = alpha,
  free_y = free_y,
  trim = trim,
  rug = rug,
  date_axis = TRUE
)
loc1 <- locations[1]
p <- do.call(plot_sample_densities, c(
  list(filter(samples, location == loc1),
       eval_data = if (is.null(eval_data)) NULL else filter(eval_data, location == loc1)),
  plot_args
))
p

# ---- Save: one page per state in a single PDF ----
out_dir <- here("results", "nowcast_sample_densities", "output")
loc_label <- if (length(locations) > 5) {
  paste0(length(locations), "locs")
} else {
  paste(locations, collapse = "-")
}
out_file <- file.path(out_dir, sprintf(
  "%s_%s_%s_%s_%s-to-%s.pdf", nowcast_date, loc_label,
  paste(models, collapse = "_vs_"), paste(clades, collapse = "-"),
  format(min(target_dates), "%m%d"), format(max(target_dates), "%m%d")
))
plots <- do.call(save_state_pages, c(
  list(samples, out_file, locations = locations, eval_data = eval_data,
       height = 5),
  plot_args
))
