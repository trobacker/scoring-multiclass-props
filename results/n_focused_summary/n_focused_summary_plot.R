#!/usr/bin/env Rscript

# N-Focused Summary Figure
#
# For a single N and the symmetric alpha=(10,10,10) setting, combines:
#   - exact Dirichlet-multinomial probability tripanel (Dispersed/Correct/Narrow),
#     with the Dirichlet density blanket in the background
#   - empirical tripanel (Dirichlet samples + exact multinomial outcome
#     probabilities as sized/colored circles)
#   - discrimination (win %) vs n_counts, filtered to the Dispersed/Narrow
#     wrong-alpha scenarios that match the tripanels (scenario3
#     extreme_dispersed, scenario5 extreme_narrow)
#
# Layout: tripanels stacked on the left, discrimination plot on the right,
# with a vertical line marking the focal N.
#
# See n_focused_summary_energy_decomposition_plot.R for the companion figure
# with mean energy score (LHS) and term-decomposed energy score (RHS).

library(ggplot2)
library(dplyr)
library(magick)

# ============================================================================
# Config
# ============================================================================

focal_n_values <- c(2, 5)

exact_dir     <- "results/plots/quantization/exact_probabilities"
empirical_dir <- "results/plots/quantization/empirical_probabilities"
data_dir      <- "results/data/alpha_10_10_10"
output_dir    <- "results/n_focused_summary"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

scenario_files <- c(
  dispersed = "scenario3_extreme_dispersed_alpha_10_10_10.rds",
  narrow    = "scenario5_extreme_narrow_alpha_10_10_10.rds"
)

# ============================================================================
# Load and prepare ES / discrimination data (shared across focal N values)
# ============================================================================

cat("Loading scenario data...\n")

scenario_data <- lapply(names(scenario_files), function(key) {
  df <- readRDS(file.path(data_dir, scenario_files[key]))
  df$scenario_key <- key
  df
})
combined_df <- do.call(rbind, scenario_data)

combined_df <- combined_df %>%
  mutate(
    scenario_label = case_when(
      scenario_key == "dispersed" ~ "Dispersed (α=0.5,0.5,0.5)",
      scenario_key == "narrow"    ~ "Narrow (α=200,200,200)"
    ),
    method_label = case_when(
      method == "mn" ~ "Multinomial",
      method == "non-mn" ~ "Non-Multinomial"
    )
  )

summary_stats <- combined_df %>%
  group_by(scenario_label, method_label, n_counts) %>%
  summarise(
    prop_correct_wins = mean(es_diff < 0),
    .groups = "drop"
  )

# ============================================================================
# Function to build the discrimination-vs-N plot with an N marker
# ============================================================================

make_discrimination_plot <- function(focal_n) {
  ggplot(summary_stats,
         aes(x = n_counts, y = prop_correct_wins,
             color = method_label,
             linetype = method_label,
             shape = method_label)) +
    geom_vline(xintercept = focal_n, linetype = "dashed", color = "#984EA3", linewidth = 1.1) +
    geom_hline(yintercept = 0.5, linetype = "dotted", color = "black", alpha = 0.5) +
    geom_line(linewidth = 1.2, alpha = 0.7) +
    geom_point(size = 3, alpha = 0.7) +
    facet_wrap(~ scenario_label, ncol = 1) +
    scale_x_log10(breaks = c(1, 5, 25, 100, 500)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = scales::percent) +
    scale_color_manual(
      values = c("Multinomial" = "#2166AC", "Non-Multinomial" = "#FF7F00"),
      name = "Method"
    ) +
    scale_linetype_manual(
      values = c("Multinomial" = "solid", "Non-Multinomial" = "solid"),
      name = "Method"
    ) +
    scale_shape_manual(
      values = c("Multinomial" = 16, "Non-Multinomial" = 17),
      name = "Method"
    ) +
    labs(
      title = "Discrimination Power vs n_counts",
      subtitle = "Proportion of simulations where correct forecast wins",
      x = "N (log scale)",
      y = "Proportion Correct Wins"
    ) +
    theme_bw(base_size = 24) +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 22),
      legend.text = element_text(size = 20),
      plot.title = element_text(face = "bold", size = 26),
      plot.subtitle = element_text(size = 20),
      axis.title = element_text(size = 22),
      axis.text = element_text(size = 19),
      strip.background = element_rect(fill = "gray90"),
      strip.text = element_text(face = "bold", size = 20)
    )
}

# ============================================================================
# Function to assemble the full composite figure for one focal N
# ============================================================================

assemble_composite <- function(focal_n) {
  cat("Assembling composite figure for N =", focal_n, "...\n")

  exact_path <- file.path(exact_dir, paste0("exact_prob_symmetric_tripanel_N", focal_n, ".png"))
  empirical_path <- file.path(empirical_dir, paste0("quantization_symmetric_tripanel_exact_N", focal_n, ".png"))

  if (!file.exists(exact_path)) stop("Missing exact tripanel PNG: ", exact_path)
  if (!file.exists(empirical_path)) stop("Missing empirical tripanel PNG: ", empirical_path)

  exact_img <- image_read(exact_path)
  empirical_img <- image_read(empirical_path)

  panel_width <- image_info(exact_img)$width
  panel_height <- image_info(exact_img)$height

  # Left column: exact tripanel over empirical tripanel
  left_col <- image_append(image_scale(c(exact_img, empirical_img), paste0(panel_width)), stack = TRUE)
  left_height <- image_info(left_col)$height

  # Right column: discrimination-vs-N, rendered tall enough to match the
  # left column's full stacked height
  disc_plot_path <- file.path(output_dir, paste0("_tmp_disc_N", focal_n, ".png"))
  ggsave(disc_plot_path, make_discrimination_plot(focal_n),
         width = panel_width / 300, height = 2 * panel_height / 300, dpi = 300, units = "in")

  disc_img <- image_read(disc_plot_path)
  right_col <- image_scale(disc_img, paste0("x", left_height))

  composite <- image_append(c(left_col, right_col), stack = FALSE)

  # Add an overall title banner
  title_img <- image_blank(image_info(composite)$width, 680, color = "white")
  title_img <- image_annotate(
    title_img,
    paste0("Symmetric α=(10,10,10): Discretization and Discrimination at N = ", focal_n),
    size = 120, gravity = "center", weight = 700, color = "black"
  )

  final <- image_append(c(title_img, composite), stack = TRUE)

  out_path <- file.path(output_dir, paste0("n_focused_summary_symmetric_N", focal_n, ".png"))
  image_write(final, out_path, density = "300x300")

  file.remove(disc_plot_path)

  cat("  Saved:", out_path, "\n")
}

# ============================================================================
# Run
# ============================================================================

for (n in focal_n_values) {
  assemble_composite(n)
}

cat("\nDone!\n")
