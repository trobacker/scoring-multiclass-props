#!/usr/bin/env Rscript
# Plot sideways densities of submitted sample clade proportions from the
# Variant Nowcast Hub, faceted by clade.
#
# For one nowcast date and each requested location, each panel (clade) shows, at each
# requested horizon (target_date - nowcast_date, in days), a density of the
# 100 sample proportions submitted by each model. All densities are drawn on
# the right of their time point (semi-transparent, so two models overlap), with
# each model's median marked at the time point itself.
#
# If final evaluation data exist for the round (hub oracle output with as_of at
# least 90 days after the nowcast date, as used by the hub's scoring), each
# time point also shows every proportion the observed counts could take,
# k / N for k = 0..N (hollow black circles), and the observed proportion
# (filled black circle). N is printed under each date on the top axis.
#
# With several locations, each location is one page of a single PDF.
#
# Optionally, only chosen clades are kept and all others are summed into
# "other" (per sample for forecasts, per date for observed counts); see
# collapse_clades() and default_clades().
#
# Usage (run from the scoring-multiclass-props repo root):
#   Rscript results/nowcast_sample_densities/plot_nowcast_sample_densities.R \
#     --nowcast-date=2026-10-14 --location=MA \
#     --models=UMass-HMLR,Hub-baseline --horizons=-28:10:7
#
#   # Several states, one page each in one PDF:
#   Rscript results/nowcast_sample_densities/plot_nowcast_sample_densities.R \
#     --nowcast-date=2026-07-22 --location=CA,NY,TX \
#     --models=UMass-HMLR,Hub-baseline
#
# Arguments:
#   --nowcast-date  Round ID (a Wednesday), YYYY-MM-DD.            [required]
#   --location      Two-letter state abbreviation(s), comma-separated, or
#                   "all" for every location the models submitted. More
#                   than one location gives a multi-page PDF.    [required]
#   --models        One or two model IDs, comma-separated.         [required]
#   --horizons      Days relative to nowcast date. Either comma-separated
#                   (e.g. -21,-14,-7,0,7) or from:to:by (e.g. -28:10:7).
#                   Hub horizons run from -31 to 10.     [default -28:7:7]
#   --hub-path      Path to the variant-nowcast-hub repo.
#                                              [default ../variant-nowcast-hub]
#   --colors        Model colours, comma-separated, in the same order as
#                   --models. Hex codes or R colour names
#                   (e.g. "#1b9e77,#d95f02" or steelblue,firebrick).
#                                               [default #2b6cb0,#c05621]
#   --alpha         Fill opacity of the densities, 0-1.          [default 0.4]
#   --free-y        If set, use free y scales per clade instead of [0, 1].
#   --no-trim       If set, extend densities to the [0, 1] bounds instead of
#                   trimming them to the range of the samples.
#   --rug           If set, draw each sample as a short tick to the left of
#                   its time point (one column of ticks per model).
#   --no-eval       If set, do not overlay evaluation (oracle) data.
#   --clades        Clades to keep, comma-separated (e.g. 25B,25C); all
#                   others are summed into "other". "auto" picks the top two
#                   (see default_clades()).        [default: all clades]
#   --date-axis     If set, label the x-axis with target dates instead of
#                   days relative to the nowcast date (N stays on top).
#   --out           Output file. A single location is saved as PNG unless
#                   the name ends in .pdf; several locations always give a
#                   PDF.            [default results/nowcast_sample_densities/
#                                    output/<date>_<locs>_<models>.png/.pdf]
#
# The functions can also be used interactively (see all_clades_by_horizon.R,
# clades_nowcast.R and focused_clades_consecutive_days.R):
#   source("results/nowcast_sample_densities/plot_nowcast_sample_densities.R")
#   s <- load_hub_samples("../variant-nowcast-hub", "2026-10-14", "MA",
#                         c("UMass-HMLR", "Hub-baseline"))
#   s <- load_hub_samples("~/GitHub/variant-nowcast-hub", "2026-10-14", "CA", c("UMass-HMLR", "Hub-baseline"));
#   plot_sample_densities(s, horizons = seq(-28, 7, by = 7))
#   s <- load_hub_samples("../variant-nowcast-hub", "2026-07-22",
#                         c("CA", "NY"), c("UMass-HMLR", "Hub-baseline"))
#   save_state_pages(s, "out.pdf", horizons = seq(-28, 7, by = 7))

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(ggplot2)
  library(ggdist)
})

# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

#' Load sample-type predictions for one round, given locations, and models.
#' `location` is a vector of two-letter abbreviations, or "all" for every
#' location in the submissions.
#' Returns columns: model_id, nowcast_date, target_date, horizon, clade,
#' location, sample_id, value.
load_hub_samples <- function(hub_path, nowcast_date, location, models) {
  all_locs <- identical(location, "all")
  nowcast_date <- as.Date(nowcast_date)
  out <- lapply(models, function(m) {
    f <- file.path(hub_path, "model-output", m,
                   sprintf("%s-%s.parquet", nowcast_date, m))
    if (!file.exists(f)) {
      stop("No submission found for model '", m, "' on ", nowcast_date,
           ":\n  ", f)
    }
    locs <- location
    read_parquet(f) |>
      filter(output_type == "sample", all_locs | location %in% locs) |>
      transmute(
        model_id = m,
        nowcast_date = as.Date(nowcast_date),
        target_date = as.Date(target_date),
        horizon = as.integer(target_date - nowcast_date),
        clade, location,
        sample_id = output_type_id,
        value
      )
  })
  samples <- bind_rows(out)
  missing <- setdiff(models, unique(samples$model_id))
  if (length(missing) > 0) {
    stop("No sample-type output for location(s) ", paste(location, collapse = ", "),
         " from: ", paste(missing, collapse = ", "))
  }
  if (!all_locs) {
    absent <- setdiff(location, unique(samples$location))
    if (length(absent) > 0) {
      warning("No samples from any model for: ", paste(absent, collapse = ", "))
    }
  }
  samples
}

#' Load final evaluation data (hub oracle output) for one round and the given
#' location(s), or "all".
#'
#' Follows the hub's scoring scheme: oracle output for a round is only final
#' once sequences are pulled at least `final_days` days after the nowcast
#' date. Earlier oracle files are interim and are rebuilt weekly.
#'
#' Returns columns location, target_date, horizon, clade, oracle_value, N, or
#' NULL (with a message) if no final evaluation data exist.
load_eval_data <- function(hub_path, nowcast_date, location,
                           final_days = 90) {
  nowcast_date <- as.Date(nowcast_date)
  f <- file.path(hub_path, "target-data", "oracle-output",
                 paste0("nowcast_date=", nowcast_date), "oracle.parquet")
  if (!file.exists(f)) {
    message("No eval data exists for nowcast date ", nowcast_date,
            " (no oracle output file).")
    return(NULL)
  }
  locs <- location
  oracle <- read_parquet(f) |>
    filter(identical(locs, "all") | location %in% locs)
  as_of <- max(as.Date(oracle$as_of))
  lag <- as.integer(as_of - nowcast_date)
  if (lag < final_days) {
    message("No eval data exists for nowcast date ", nowcast_date,
            ": oracle output is interim (as_of ", as_of, ", ", lag,
            " days after nowcast date). Final data are expected as of ",
            nowcast_date + final_days, ".")
    return(NULL)
  }
  oracle |>
    transmute(
      location,
      target_date = as.Date(target_date),
      horizon = as.integer(target_date - nowcast_date),
      clade, oracle_value
    ) |>
    group_by(location, target_date) |>
    mutate(N = sum(oracle_value)) |>
    ungroup()
}

#' Parse a horizons string: "a:b:by", "a:b", or "h1,h2,...".
parse_horizons <- function(x) {
  if (grepl(":", x)) {
    p <- as.integer(strsplit(x, ":", fixed = TRUE)[[1]])
    by <- if (length(p) == 3) p[3] else 1L
    seq(p[1], p[2], by = by)
  } else {
    as.integer(strsplit(x, ",", fixed = TRUE)[[1]])
  }
}

#' Pick default clades to plot for a round: the `n` named clades (excluding
#' "other" and "recombinant") with the most sequences nationally in the 3
#' weeks before the round, from the hub's modeled-clades file. Rounds before
#' 2025-04-16 have no counts in that file; for those, fall back to the clades
#' with the highest mean predicted proportion in `samples`.
default_clades <- function(hub_path, nowcast_date, samples = NULL, n = 2) {
  exclude <- c("other", "recombinant")
  f <- file.path(hub_path, "auxiliary-data", "modeled-clades",
                 paste0(as.Date(nowcast_date), ".json"))
  counts <- if (file.exists(f)) {
    meta <- jsonlite::fromJSON(f)
    by_clade <- unlist(meta$meta$sequence_counts$sequences_by_clade)
    # Counts can include clades that were not modeled this round.
    by_clade[names(by_clade) %in% setdiff(meta$clades, exclude)]
  }
  if (length(counts) > 0) {
    top <- names(sort(counts, decreasing = TRUE))[seq_len(min(n, length(counts)))]
    message("Default clades (hub national sequence counts, 3 weeks before ",
            "the round): ", paste(top, collapse = ", "))
    return(top)
  }
  if (is.null(samples)) {
    stop("No sequence counts in ", f, "; pass `samples` to rank clades by ",
         "predicted proportion instead.")
  }
  top <- samples |>
    filter(!clade %in% exclude) |>
    group_by(clade) |>
    summarise(mean = mean(value), .groups = "drop") |>
    slice_max(mean, n = n, with_ties = FALSE) |>
    pull(clade)
  message("Default clades (no hub sequence counts for this round; highest ",
          "mean predicted proportion): ", paste(top, collapse = ", "))
  top
}

#' Keep `clades` and sum every other clade into "other".
#'
#' Works for samples (value_col = "value", summed within each sample, so each
#' sample still sums to 1) and for eval data (value_col = "oracle_value").
#' Clade becomes a factor ordered as `clades`, then "other", which sets the
#' facet order.
collapse_clades <- function(df, clades, value_col = "value") {
  missing <- setdiff(clades, unique(df$clade))
  if (length(missing) > 0) {
    warning("Clades not in data: ", paste(missing, collapse = ", "))
  }
  # Summarising sorts rows by the grouping columns; keep the original model
  # order, which plot_sample_densities() uses for colours and the legend.
  model_order <- if ("model_id" %in% names(df)) unique(df$model_id)
  out <- df |>
    mutate(clade = ifelse(clade %in% clades, clade, "other")) |>
    group_by(across(-all_of(value_col))) |>
    summarise(!!value_col := sum(.data[[value_col]]), .groups = "drop") |>
    mutate(clade = factor(clade, levels = unique(c(clades, "other"))))
  if (!is.null(model_order)) {
    out <- arrange(out, match(model_id, model_order))
  }
  out
}

# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

#' Facet-by-clade plot of sideways sample densities at each horizon.
#'
#' @param samples Output of load_hub_samples().
#' @param horizons Integer days relative to nowcast date to show.
#' @param model_colors Colours for the models, in model order. May be named
#'   by model ID. Defaults to blue and orange.
#' @param alpha Fill opacity of the densities.
#' @param free_y Use free y scales per clade (default: fixed [0, 1]).
#' @param trim Trim each density to the range of its samples (ggdist default).
#'   If FALSE, densities extend to the [0, 1] bounds.
#' @param rug Draw each sample as a short tick left of its time point.
#' @param eval_data Output of load_eval_data(), or NULL to skip.
#' @param date_axis Label the x-axis with target dates instead of days
#'   relative to the nowcast date. N (with eval data) stays on the top axis.
#' @param slab_scale Fraction of the spacing between horizons that a density
#'   may occupy to the right of its time point.
plot_sample_densities <- function(samples, horizons,
                                  model_colors = c("#2b6cb0", "#c05621"),
                                  alpha = 0.4, free_y = FALSE,
                                  trim = TRUE, rug = FALSE,
                                  eval_data = NULL, date_axis = FALSE,
                                  slab_scale = 0.9) {
  models <- unique(samples$model_id)
  if (length(models) > 2) stop("At most 2 models are supported for now.")
  if (length(model_colors) < length(models)) {
    stop("Need ", length(models), " colours, got ", length(model_colors))
  }
  if (is.null(names(model_colors))) {
    model_colors <- setNames(model_colors[seq_along(models)], models)
  }

  available <- sort(unique(samples$horizon))
  dropped <- setdiff(horizons, available)
  if (length(dropped) > 0) {
    warning("Horizons not in submissions, dropped: ",
            paste(dropped, collapse = ", "),
            " (available: ", min(available), " to ", max(available), ")")
  }
  horizons <- intersect(horizons, available)
  if (length(horizons) == 0) stop("None of the requested horizons are available.")

  d <- samples |>
    filter(horizon %in% horizons) |>
    mutate(model_id = factor(model_id, levels = models))

  # Width available to the right of a time point, in x (day) units.
  spacing <- if (length(horizons) > 1) min(diff(sort(horizons))) else 1

  p <- ggplot(d, aes(x = horizon, y = value, fill = model_id,
                     colour = model_id)) +
    stat_slab(
      aes(group = interaction(model_id, horizon)),
      side = "right",
      scale = slab_scale,
      # Normalize within each density so every horizon fills its slot;
      # shape is comparable, height across horizons is not.
      normalize = "groups",
      density = density_bounded(bounds = c(0, 1), bandwidth = "nrd0"),
      trim = trim,
      alpha = alpha, linewidth = 0.3,
      width = spacing
    )

  if (rug) {
    # One column of ticks per model, stacked leftwards from the time point.
    tick_len <- 0.06 * spacing
    gap <- 0.03 * spacing
    rug_d <- d |>
      mutate(
        i = as.integer(model_id) - 1L,
        xend = horizon - gap - i * (tick_len + gap),
        x = xend - tick_len
      )
    p <- p + geom_segment(
      data = rug_d,
      aes(x = x, xend = xend, y = value, yend = value, colour = model_id),
      inherit.aes = FALSE, alpha = 0.5, linewidth = 0.25,
      show.legend = FALSE
    )
  }

  nowcast_date <- unique(samples$nowcast_date)
  location <- unique(samples$location)
  date_labels <- format(nowcast_date + horizons, "%b %d")
  # Top axis: dates (and N) in horizon mode; only N in date mode.
  top_labels <- if (date_axis) NULL else date_labels
  subtitle <- "Densities of 100 submitted samples per model and target date; points are medians"
  # Observed-proportion layer, added after the medians so it is drawn on top.
  observed_layer <- NULL

  if (!is.null(eval_data)) {
    ev <- eval_data |>
      filter(horizon %in% horizons, clade %in% unique(d$clade))
    # Every proportion the observed counts could take at each time point.
    support <- ev |>
      filter(N > 0) |>
      distinct(clade, horizon, N) |>
      tidyr::uncount(N + 1L, .id = "k") |>
      mutate(value = (k - 1L) / N)
    observed <- ev |>
      filter(N > 0) |>
      mutate(value = oracle_value / N)

    p <- p +
      geom_point(data = support, aes(x = horizon, y = value),
                 inherit.aes = FALSE, shape = 21, fill = NA,
                 colour = "black", size = 1.3, stroke = 0.35)
    observed_layer <- geom_point(data = observed, aes(x = horizon, y = value),
                                 inherit.aes = FALSE, shape = 16,
                                 colour = "black", size = 1.6)

    n_by_h <- distinct(ev, horizon, N)
    n_lab <- n_by_h$N[match(horizons, n_by_h$horizon)]
    n_text <- paste0("N=", ifelse(is.na(n_lab), "?", n_lab))
    top_labels <- if (date_axis) n_text else paste0(top_labels, "\n", n_text)
    subtitle <- paste0(
      subtitle,
      "\nHollow circles: possible observed proportions k/N; filled circle: observed proportion"
    )
  }

  # Only mark the nowcast date when it falls within the plotted dates;
  # otherwise the line would stretch the x-axis out to it.
  show_nowcast <- min(horizons) <= 0 && max(horizons) >= 0
  nowcast_line <- if (show_nowcast) {
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40",
               linewidth = 0.3)
  }
  if (date_axis && show_nowcast) {
    subtitle <- paste0(subtitle, "; dashed line: nowcast date")
  }

  p +
    stat_summary(aes(group = interaction(model_id, horizon)),
                 fun = median, geom = "point", shape = 21, size = 1.3,
                 stroke = 0.3, colour = "white", show.legend = FALSE) +
    observed_layer +
    nowcast_line +
    facet_wrap(~ clade, scales = if (free_y) "free_y" else "fixed") +
    scale_x_continuous(
      breaks = horizons,
      labels = if (date_axis) date_labels else waiver(),
      sec.axis = if (is.null(top_labels)) {
        waiver()
      } else {
        dup_axis(name = NULL, breaks = horizons, labels = top_labels)
      }
    ) +
    # Zoom rather than set scale limits, so densities touching 0 or 1 are not
    # dropped.
    coord_cartesian(ylim = if (free_y) NULL else c(0, 1)) +
    scale_fill_manual(values = model_colors, name = NULL) +
    scale_colour_manual(values = model_colors, name = NULL) +
    labs(
      x = if (date_axis) "Target date" else "Days relative to nowcast date",
      y = "Clade proportion",
      title = sprintf("Submitted sample distributions: %s, nowcast date %s",
                      location, nowcast_date),
      subtitle = subtitle
    ) +
    theme_bw(base_size = 10) +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank(),
      axis.text.x.top = if (date_axis) {
        element_text(size = 7)
      } else {
        element_text(size = 6, angle = 45, hjust = 0)
      },
      axis.text.x.bottom = element_text(size = 7)
    )
}

#' Save one page per location to a single PDF.
#'
#' @param samples Output of load_hub_samples() covering one or more locations.
#' @param file Output PDF path.
#' @param locations Locations to plot, in page order. Defaults to all
#'   locations in `samples`, alphabetically.
#' @param eval_data Output of load_eval_data(), or NULL.
#' @param width,height Page size in inches.
#' @param ... Passed to plot_sample_densities() (horizons, model_colors, ...).
#' @return Invisibly, a named list of the ggplot objects.
save_state_pages <- function(samples, file,
                             locations = sort(unique(samples$location)),
                             eval_data = NULL, width = 11, height = 8, ...) {
  models <- unique(samples$model_id)
  plots <- list()
  for (loc in locations) {
    s <- filter(samples, location == loc)
    have <- unique(s$model_id)
    if (length(have) < length(models)) {
      message("Skipping ", loc, ": no samples from ",
              paste(setdiff(models, have), collapse = ", "))
      next
    }
    ev <- if (is.null(eval_data)) NULL else filter(eval_data, location == loc)
    plots[[loc]] <- tryCatch(
      plot_sample_densities(s, eval_data = ev, ...),
      error = function(e) {
        message("Skipping ", loc, ": ", conditionMessage(e))
        NULL
      }
    )
  }
  if (length(plots) == 0) stop("No locations could be plotted.")

  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  pdf(file, width = width, height = height)
  on.exit(dev.off())
  for (p in plots) print(p)
  message("Saved ", length(plots), " page(s): ", file)
  invisible(plots)
}

# ---------------------------------------------------------------------------
# Command line entry point
# ---------------------------------------------------------------------------

parse_args <- function(args) {
  flags <- c("free-y", "no-trim", "rug", "no-eval", "date-axis")
  out <- list()
  for (a in args) {
    a <- sub("^--", "", a)
    if (a %in% flags) {
      out[[a]] <- TRUE
    } else {
      kv <- regmatches(a, regexpr("=", a), invert = TRUE)[[1]]
      out[[kv[1]]] <- kv[2]
    }
  }
  out
}

if (sys.nframe() == 0) {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  for (req in c("nowcast-date", "location", "models")) {
    if (is.null(args[[req]])) stop("Missing required argument --", req)
  }
  hub_path <- args[["hub-path"]] %||% "../variant-nowcast-hub"
  models <- strsplit(args[["models"]], ",", fixed = TRUE)[[1]]
  horizons <- parse_horizons(args[["horizons"]] %||% "-28:7:7")
  free_y <- isTRUE(args[["free-y"]])
  model_colors <- if (is.null(args[["colors"]])) {
    c("#2b6cb0", "#c05621")
  } else {
    strsplit(args[["colors"]], ",", fixed = TRUE)[[1]]
  }
  alpha <- as.numeric(args[["alpha"]] %||% 0.4)
  trim <- !isTRUE(args[["no-trim"]])
  rug <- isTRUE(args[["rug"]])
  locations <- strsplit(args[["location"]], ",", fixed = TRUE)[[1]]

  samples <- load_hub_samples(hub_path, args[["nowcast-date"]],
                              locations, models)
  eval_data <- if (isTRUE(args[["no-eval"]])) {
    NULL
  } else {
    load_eval_data(hub_path, args[["nowcast-date"]], locations)
  }
  if (identical(locations, "all")) locations <- sort(unique(samples$location))

  clades <- args[["clades"]]
  if (!is.null(clades)) {
    clades <- if (clades == "auto") {
      default_clades(hub_path, args[["nowcast-date"]], samples)
    } else {
      strsplit(clades, ",", fixed = TRUE)[[1]]
    }
    samples <- collapse_clades(samples, clades, "value")
    if (!is.null(eval_data)) {
      eval_data <- collapse_clades(eval_data, clades, "oracle_value")
    }
  }

  multi <- length(locations) > 1
  loc_label <- if (length(locations) > 5) {
    paste0(length(locations), "locs")
  } else {
    paste(locations, collapse = "-")
  }
  out <- args[["out"]] %||% file.path(
    "results/nowcast_sample_densities/output",
    sprintf("%s_%s_%s%s%s%s%s.%s", args[["nowcast-date"]], loc_label,
            paste(models, collapse = "_vs_"),
            if (is.null(clades)) "" else paste0("_", paste(clades, collapse = "-")),
            if (free_y) "_freey" else "",
            if (trim) "" else "_notrim", if (rug) "_rug" else "",
            if (multi) "pdf" else "png")
  )
  plot_args <- list(horizons = horizons, model_colors = model_colors,
                    alpha = alpha, free_y = free_y, trim = trim, rug = rug,
                    date_axis = isTRUE(args[["date-axis"]]))

  if (multi || grepl("\\.pdf$", out)) {
    do.call(save_state_pages, c(list(samples, out, locations = locations,
                                     eval_data = eval_data), plot_args))
  } else {
    p <- do.call(plot_sample_densities,
                 c(list(samples, eval_data = eval_data), plot_args))
    dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
    ggsave(out, p, width = 11, height = 8, dpi = 200)
    cat("Saved:", out, "\n")
  }
}
