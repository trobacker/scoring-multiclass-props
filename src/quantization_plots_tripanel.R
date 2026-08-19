#!/usr/bin/env Rscript

# Quantization Effect Visualization - Three Panel Version
# Shows dispersed, correct, and narrow forecasts side-by-side
# to illustrate the full spectrum of forecast variance

library(mc2d)
library(brms)
library(Ternary)
library(viridisLite)

# Source helper functions
source("src/helper-functions.R")

# Set output directory (empirical probabilities subdirectory)
output_dir <- "results/plots/quantization/empirical_probabilities"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("=========================================\n")
cat("QUANTIZATION TRIPANEL VISUALIZATION\n")
cat("=========================================\n\n")

# Set random seed for reproducibility
set.seed(42)

# Parameters
n_samples <- 100  # Number of Dirichlet samples to draw
n_values <- c(1, 2, 3, 4, 5, 10, 25, 50, 100)  # Different sample sizes to visualize

# Define tripanel scenarios (dispersed - correct - narrow)
tripanel_scenarios <- list(
  # Symmetric tripanel
  list(
    name = "symmetric_tripanel",
    alpha_dispersed = c(0.5, 0.5, 0.5),
    alpha_correct = c(10, 10, 10),
    alpha_narrow = c(200, 200, 200),
    dispersed_label = "Dispersed: α=(0.5,0.5,0.5)",
    correct_label = "Correct: α=(10,10,10)",
    narrow_label = "Narrow: α=(200,200,200)"
  ),
  # Asymmetric tripanel
  list(
    name = "asymmetric_tripanel",
    alpha_dispersed = c(0.4, 0.6, 3),
    alpha_correct = c(2, 3, 15),
    alpha_narrow = c(200, 300, 1500),
    dispersed_label = "Dispersed: α=(0.4,0.6,3)",
    correct_label = "Correct: α=(2,3,15)",
    narrow_label = "Narrow: α=(200,300,1500)"
  )
)

# Function to create three-panel ternary comparison
create_tripanel_comparison <- function(n, alpha_dispersed, alpha_correct, alpha_narrow,
                                       dispersed_label, correct_label, narrow_label,
                                       scenario_name) {

  cat("  Creating tripanel for N =", n, "...\n")

  # Generate Dirichlet samples
  dirichlet_dispersed <- brms::rdirichlet(n_samples, alpha_dispersed)
  dirichlet_correct <- brms::rdirichlet(n_samples, alpha_correct)
  dirichlet_narrow <- brms::rdirichlet(n_samples, alpha_narrow)

  # Generate multinomial resamples
  mn_dispersed <- t(sapply(1:n_samples, function(i) {
    counts <- mc2d::rmultinomial(1, n, dirichlet_dispersed[i, ])
    counts / n  # Convert to proportions
  }))

  mn_correct <- t(sapply(1:n_samples, function(i) {
    counts <- mc2d::rmultinomial(1, n, dirichlet_correct[i, ])
    counts / n
  }))

  mn_narrow <- t(sapply(1:n_samples, function(i) {
    counts <- mc2d::rmultinomial(1, n, dirichlet_narrow[i, ])
    counts / n
  }))

  # Count frequencies of multinomial outcomes for color mapping
  mn_dispersed_counts <- table(apply(mn_dispersed, 1, paste, collapse = "_"))
  mn_correct_counts <- table(apply(mn_correct, 1, paste, collapse = "_"))
  mn_narrow_counts <- table(apply(mn_narrow, 1, paste, collapse = "_"))

  # Create color mapping for multinomial samples based on frequency
  n_colors <- 100
  color_palette <- viridis(n_colors, begin = 0.1, end = 0.95, option = "viridis")

  # Function to get colors based on frequency
  get_freq_colors <- function(data, counts_table) {
    keys <- apply(data, 1, paste, collapse = "_")
    freqs <- as.numeric(counts_table[keys])

    # Handle case where all frequencies are the same
    if (max(freqs) == min(freqs)) {
      return(rep(color_palette[n_colors %/% 2], length(freqs)))
    }

    breaks <- seq(min(freqs) - 0.1, max(freqs) + 0.1, length.out = n_colors + 1)
    colors <- sapply(freqs, function(f) {
      idx <- findInterval(f, breaks, all.inside = TRUE)
      idx <- min(idx, n_colors)
      color_palette[idx]
    })
    return(colors)
  }

  mn_dispersed_colors <- get_freq_colors(mn_dispersed, mn_dispersed_counts)
  mn_correct_colors <- get_freq_colors(mn_correct, mn_correct_counts)
  mn_narrow_colors <- get_freq_colors(mn_narrow, mn_narrow_counts)

  # Create PNG with layout for legend
  filename <- paste0("quantization_", scenario_name, "_N", n, ".png")
  png(file.path(output_dir, filename), width = 2200, height = 700, res = 150)

  # Create layout: 3 ternary plots + 1 legend panel
  layout(matrix(c(1, 2, 3, 4), nrow = 1), widths = c(3, 3, 3, 1))

  par(mar = c(1, 1, 2.5, 1))

  # LEFT PANEL: Dispersed forecast
  TernaryPlot(
    alab = "p1",
    blab = "p2",
    clab = "p3",
    grid.lines = 5,
    point = "right",
    lab.cex = 0.8,
    grid.minor.lines = 0,
    grid.lty = "solid",
    col = rgb(0.9, 0.9, 0.9),
    axis.cex = 0.7,
    padding = 0.08
  )

  # Plot original Dirichlet samples (dodgerblue circles with transparency)
  TernaryPoints(
    dirichlet_dispersed,
    pch = 16,
    cex = 1.8,
    col = rgb(30/255, 144/255, 255/255, 0.5)  # dodgerblue with alpha=0.5
  )

  # Plot multinomial resampled points (colored triangles, slightly larger)
  TernaryPoints(
    mn_dispersed,
    pch = 17,
    cex = 2.2,
    col = mn_dispersed_colors
  )

  title(main = "Dispersed", cex.main = 1.5, font.main = 2, line = -2)

  # MIDDLE PANEL: Correct forecast
  TernaryPlot(
    alab = "p1",
    blab = "p2",
    clab = "p3",
    grid.lines = 5,
    point = "right",
    lab.cex = 0.8,
    grid.minor.lines = 0,
    grid.lty = "solid",
    col = rgb(0.9, 0.9, 0.9),
    axis.cex = 0.7,
    padding = 0.08
  )

  # Plot original Dirichlet samples (dodgerblue circles with transparency)
  TernaryPoints(
    dirichlet_correct,
    pch = 16,
    cex = 1.8,
    col = rgb(30/255, 144/255, 255/255, 0.5)
  )

  # Plot multinomial resampled points
  TernaryPoints(
    mn_correct,
    pch = 17,
    cex = 2.2,
    col = mn_correct_colors
  )

  title(main = "Correct", cex.main = 1.5, font.main = 2, line = -2)

  # RIGHT PANEL: Narrow forecast
  TernaryPlot(
    alab = "p1",
    blab = "p2",
    clab = "p3",
    grid.lines = 5,
    point = "right",
    lab.cex = 0.8,
    grid.minor.lines = 0,
    grid.lty = "solid",
    col = rgb(0.9, 0.9, 0.9),
    axis.cex = 0.7,
    padding = 0.08
  )

  # Plot original Dirichlet samples (dodgerblue circles with transparency)
  TernaryPoints(
    dirichlet_narrow,
    pch = 16,
    cex = 1.8,
    col = rgb(30/255, 144/255, 255/255, 0.5)
  )

  # Plot multinomial resampled points
  TernaryPoints(
    mn_narrow,
    pch = 17,
    cex = 2.2,
    col = mn_narrow_colors
  )

  title(main = "Narrow", cex.main = 1.5, font.main = 2, line = -2)

  # LEGEND PANEL
  par(mar = c(2, 0, 3, 2))
  plot.new()

  # Add distribution legend
  legend("topleft",
         legend = c("Dirichlet Samples", "Multinomial Resampled"),
         pch = c(16, 17),
         pt.cex = c(1.2, 1.8),
         col = c(rgb(30/255, 144/255, 255/255, 0.5), color_palette[50]),
         bty = "n",
         cex = 0.95,
         title = "Distribution",
         title.adj = 0)

  # Add color scale for frequency
  # Determine the frequency range across all three forecasts
  all_freqs <- c(as.numeric(mn_dispersed_counts),
                 as.numeric(mn_correct_counts),
                 as.numeric(mn_narrow_counts))
  max_freq <- max(all_freqs)
  min_freq <- min(all_freqs)

  # Create a vertical color bar
  n_bar_steps <- 50
  bar_colors <- viridis(n_bar_steps, begin = 0.1, end = 0.95, option = "viridis")

  # Position for color bar
  bar_left <- 0.15
  bar_right <- 0.35
  bar_bottom <- 0.15
  bar_top <- 0.75

  # Draw color bar
  for (i in 1:(n_bar_steps - 1)) {
    y_bottom <- bar_bottom + (i - 1) * (bar_top - bar_bottom) / n_bar_steps
    y_top <- bar_bottom + i * (bar_top - bar_bottom) / n_bar_steps
    rect(bar_left, y_bottom, bar_right, y_top,
         col = bar_colors[i], border = NA)
  }
  rect(bar_left, bar_bottom, bar_right, bar_top, border = "black", lwd = 1)

  # Add labels to color bar
  text(bar_right + 0.05, bar_top, sprintf("%d", max_freq), adj = 0, cex = 0.8)
  text(bar_right + 0.05, (bar_top + bar_bottom) / 2,
       sprintf("%d", round((max_freq + min_freq) / 2)), adj = 0, cex = 0.8)
  text(bar_right + 0.05, bar_bottom, sprintf("%d", min_freq), adj = 0, cex = 0.8)

  # Add title for color scale
  text((bar_left + bar_right) / 2, bar_top + 0.08, "Frequency",
       cex = 0.95, font = 2, adj = 0.5)

  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 7, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot.new()
  mtext(paste0("N = ", n),
        side = 3, line = 4, cex = 1.4, font = 2)

  dev.off()

  cat("    Saved:", filename, "\n")
}

# Generate plots for each scenario and N value
for (scenario in tripanel_scenarios) {
  cat("\nScenario:", scenario$name, "\n")
  cat("  Three-panel:", scenario$dispersed_label, "|",
      scenario$correct_label, "|", scenario$narrow_label, "\n")

  for (n in n_values) {
    create_tripanel_comparison(
      n,
      scenario$alpha_dispersed,
      scenario$alpha_correct,
      scenario$alpha_narrow,
      scenario$dispersed_label,
      scenario$correct_label,
      scenario$narrow_label,
      scenario$name
    )
  }
}

cat("\n=========================================\n")
cat("SUMMARY\n")
cat("=========================================\n")
cat("Generated tripanel quantization plots for:\n")
cat("  Scenarios:", length(tripanel_scenarios), "\n")
cat("  N values:", paste(n_values, collapse = ", "), "\n")
cat("  Total plots:", length(tripanel_scenarios) * length(n_values), "\n")
cat("\nAll plots saved to:", output_dir, "\n\n")
cat("Done!\n")
