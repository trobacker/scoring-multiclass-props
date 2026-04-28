#!/usr/bin/env Rscript

# Quantization Effect: Exact Probabilities Visualization
# Shows the exact probabilities of all possible multinomial outcomes
# using the Dirichlet-multinomial distribution

library(extraDistr)
library(brms)
library(Ternary)
library(viridisLite)

# Source helper functions
source("src/helper-functions.R")

# ============================================================================
# Helper Function: Create Dirichlet Density Grid
# ============================================================================

create_dirichlet_grid <- function(alpha, grid_res = 150) {
  grid_coords <- expand.grid(
    a = seq(0, 1, length.out = grid_res),
    b = seq(0, 1, length.out = grid_res)
  )
  grid_coords$c <- 1 - grid_coords$a - grid_coords$b

  # Filter to simplex interior
  grid_coords <- grid_coords[
    grid_coords$c >= 0.001 & grid_coords$c <= 0.999 &
    grid_coords$a >= 0.001 & grid_coords$a <= 0.999 &
    grid_coords$b >= 0.001 & grid_coords$b <= 0.999,
  ]

  # Renormalize
  row_sums <- grid_coords$a + grid_coords$b + grid_coords$c
  grid_coords$a <- grid_coords$a / row_sums
  grid_coords$b <- grid_coords$b / row_sums
  grid_coords$c <- grid_coords$c / row_sums

  # Calculate Dirichlet density
  grid_coords$density <- dd_func(grid_coords$a, grid_coords$b, grid_coords$c, alpha)

  # Map to colors (viridis palette - brighter and colorblind-friendly)
  n_dens_colors <- 100
  dens_colors <- viridisLite::viridis(n_dens_colors, alpha = 0.7, begin = 0, end = 1, direction = 1)

  # Handle edge case where density values are too similar
  if (max(grid_coords$density) - min(grid_coords$density) < 1e-10) {
    grid_coords$color <- rep(dens_colors[50], nrow(grid_coords))
  } else {
    dens_breaks <- seq(min(grid_coords$density),
                       max(grid_coords$density),
                       length.out = n_dens_colors + 1)
    grid_coords$color <- dens_colors[cut(grid_coords$density,
                                         breaks = dens_breaks,
                                         labels = FALSE,
                                         include.lowest = TRUE)]
  }

  return(grid_coords)
}

# Set output directory (subdirectory of quantization plots)
output_dir <- "results/plots/quantization/exact_probabilities"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("=========================================\n")
cat("QUANTIZATION EXACT PROBABILITIES\n")
cat("=========================================\n\n")

# Parameters - limit to N <= 10 as requested
n_values <- c(1, 2, 3, 4, 5, 10)

# Define scenarios to visualize (same as original quantization plots)
scenarios <- list(
  # Symmetric scenarios
  list(
    name = "symmetric_narrow",
    alpha_correct = c(10, 10, 10),
    alpha_wrong = c(200, 200, 200),
    correct_label = "Correct: α=(10,10,10)",
    wrong_label = "Narrow: α=(200,200,200)"
  ),
  list(
    name = "symmetric_dispersed",
    alpha_correct = c(10, 10, 10),
    alpha_wrong = c(0.5, 0.5, 0.5),
    correct_label = "Correct: α=(10,10,10)",
    wrong_label = "Dispersed: α=(0.5,0.5,0.5)"
  ),
  # Asymmetric scenarios
  list(
    name = "asymmetric_misspecified",
    alpha_correct = c(2, 3, 15),
    alpha_wrong = c(15, 3, 2),
    correct_label = "Correct: α=(2,3,15)",
    wrong_label = "Misspecified: α=(15,3,2)"
  ),
  list(
    name = "asymmetric_narrow",
    alpha_correct = c(2, 3, 15),
    alpha_wrong = c(200, 300, 1500),
    correct_label = "Correct: α=(2,3,15)",
    wrong_label = "Narrow: α=(200,300,1500)"
  ),
  list(
    name = "asymmetric_dispersed",
    alpha_correct = c(2, 3, 15),
    alpha_wrong = c(0.4, 0.6, 3),
    correct_label = "Correct: α=(2,3,15)",
    wrong_label = "Dispersed: α=(0.4,0.6,3)"
  )
)

# Function to compute exact probabilities for all multinomial outcomes
compute_exact_probabilities <- function(n, alpha) {
  # Enumerate all possible multinomial outcomes (hardcoded for 3 classes)
  outcomes <- enumerate_multinomial_sample_space(n)

  # Compute Dirichlet-multinomial probability for each outcome
  probs <- apply(outcomes, 1, function(counts) {
    extraDistr::ddirmnom(x = counts, size = n, alpha = alpha)
  })

  # Convert counts to proportions for plotting
  proportions <- outcomes / n

  return(list(
    proportions = proportions,
    probabilities = probs,
    counts = outcomes
  ))
}

# Function to create side-by-side ternary comparison with exact probabilities
create_probability_comparison <- function(n, alpha_correct, alpha_wrong,
                                         correct_label, wrong_label,
                                         scenario_name) {

  cat("  Creating exact probability plot for N =", n, "...\n")

  # Create Dirichlet density grids for background
  grid_correct <- create_dirichlet_grid(alpha_correct)
  grid_wrong <- create_dirichlet_grid(alpha_wrong)

  # Compute exact probabilities
  exact_correct <- compute_exact_probabilities(n, alpha_correct)
  exact_wrong <- compute_exact_probabilities(n, alpha_wrong)

  # Create color mapping based on probability (use viridis to match density blanket)
  n_colors <- 100
  color_palette <- viridis(n_colors, alpha = 1.0, begin = 0, end = 1, direction = 1)

  # Function to map probabilities to colors
  get_prob_colors <- function(probs) {
    if (length(probs) == 1) {
      return(color_palette[n_colors %/% 2])
    }

    # Use log scale if probabilities span multiple orders of magnitude
    log_probs <- log10(probs + 1e-10)  # Add small value to avoid log(0)
    breaks <- seq(min(log_probs), max(log_probs), length.out = n_colors + 1)

    colors <- sapply(log_probs, function(lp) {
      idx <- findInterval(lp, breaks, all.inside = TRUE)
      idx <- min(idx, n_colors)
      color_palette[idx]
    })
    return(colors)
  }

  correct_colors <- get_prob_colors(exact_correct$probabilities)
  wrong_colors <- get_prob_colors(exact_wrong$probabilities)

  # Scale point sizes by probability (cube root for better visibility)
  scale_size <- function(probs) {
    base_size <- 1.5
    max_size <- 4.0
    scaled <- (probs / max(probs))^(1/3)  # Cube root for better visual distribution
    base_size + scaled * (max_size - base_size)
  }

  correct_sizes <- scale_size(exact_correct$probabilities)
  wrong_sizes <- scale_size(exact_wrong$probabilities)

  # Create PNG with layout for legend
  filename <- paste0("exact_prob_", scenario_name, "_N", n, ".png")
  png(file.path(output_dir, filename), width = 1600, height = 700, res = 150)

  # Create layout: 2 ternary plots + 1 legend panel
  layout(matrix(c(1, 2, 3), nrow = 1), widths = c(3, 3, 1))

  # LEFT PANEL: Correct forecast
  par(mar = c(2, 2, 2, 2))
  TernaryPlot(
    alab = "p1 \u2192",
    blab = "\u2190 p2",
    clab = "p3 \u2192",
    grid.lines = 5,
    point = "right",
    lab.cex = 0.8,
    grid.minor.lines = 0,
    grid.lty = "solid",
    col = rgb(0.9, 0.9, 0.9),
    grid.col = "white",
    axis.col = rgb(0.6, 0.6, 0.6),
    ticks.col = rgb(0.6, 0.6, 0.6),
    axis.cex = 0.7,
    axis.rotate = FALSE,
    padding = 0.08
  )

  # Add Dirichlet density blanket
  AddToTernary(graphics::points,
               cbind(grid_correct$a, grid_correct$b, grid_correct$c),
               pch = 15,
               col = grid_correct$color,
               cex = 0.5)

  # Plot exact probabilities (colored circles, size = probability)
  AddToTernary(graphics::points,
               exact_correct$proportions,
               pch = 21,  # Circle with border
               cex = correct_sizes,
               col = "black",  # Border color
               bg = correct_colors,  # Fill color
               lwd = 0.8)

  # Subplot title position (line = -0.5 means slightly below default)
  # Lower values move title closer to plot, higher values move it up
  title(main = correct_label, cex.main = 1.1, font.main = 2, line = -0.5)

  # RIGHT PANEL: Wrong forecast
  par(mar = c(2, 2, 2, 2))
  TernaryPlot(
    alab = "p1 \u2192",
    blab = "\u2190 p2",
    clab = "p3 \u2192",
    grid.lines = 5,
    point = "right",
    lab.cex = 0.8,
    grid.minor.lines = 0,
    grid.lty = "solid",
    col = rgb(0.9, 0.9, 0.9),
    grid.col = "white",
    axis.col = rgb(0.6, 0.6, 0.6),
    ticks.col = rgb(0.6, 0.6, 0.6),
    axis.cex = 0.7,
    axis.rotate = FALSE,
    padding = 0.08
  )

  # Add Dirichlet density blanket
  AddToTernary(graphics::points,
               cbind(grid_wrong$a, grid_wrong$b, grid_wrong$c),
               pch = 15,
               col = grid_wrong$color,
               cex = 0.5)

  # Plot exact probabilities
  AddToTernary(graphics::points,
               exact_wrong$proportions,
               pch = 21,
               cex = wrong_sizes,
               col = "black",
               bg = wrong_colors,
               lwd = 0.8)

  # Subplot title position (line = -0.5 means slightly below default)
  title(main = wrong_label, cex.main = 1.1, font.main = 2, line = -0.5)

  # LEGEND PANEL
  par(mar = c(2, 0, 3, 2))
  plot.new()

  # Add distribution legend
  legend("topleft",
         legend = c("Dirichlet Prior", "Multinomial Outcomes"),
         pch = c(15, 21),
         pt.cex = c(1, 2),
         col = c(viridis(1, alpha = 0.7, begin = 0.5, end = 0.5, direction = 1), "black"),
         pt.bg = c(NA, viridis(1, alpha = 1.0, begin = 0.8, end = 0.8, direction = 1)),
         bty = "n",
         cex = 0.95,
         title = "Distribution",
         title.adj = 0)

  # Add color scale for probabilities
  # Create a vertical color bar (use viridis to match density blanket)
  n_bar_steps <- 50
  bar_heights <- seq(0, 1, length.out = n_bar_steps)
  bar_colors <- viridis(n_bar_steps, alpha = 1.0, begin = 0, end = 1, direction = 1)

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
  max_prob <- max(c(exact_correct$probabilities, exact_wrong$probabilities))
  min_prob <- min(c(exact_correct$probabilities, exact_wrong$probabilities))

  # Use log scale for labels if range is large
  if (max_prob / min_prob > 100) {
    text(bar_right + 0.05, bar_top, sprintf("%.1e", max_prob), adj = 0, cex = 0.8)
    text(bar_right + 0.05, (bar_top + bar_bottom) / 2,
         sprintf("%.1e", 10^((log10(max_prob) + log10(min_prob)) / 2)), adj = 0, cex = 0.8)
    text(bar_right + 0.05, bar_bottom, sprintf("%.1e", min_prob), adj = 0, cex = 0.8)
  } else {
    text(bar_right + 0.05, bar_top, sprintf("%.4f", max_prob), adj = 0, cex = 0.8)
    text(bar_right + 0.05, (bar_top + bar_bottom) / 2,
         sprintf("%.4f", (max_prob + min_prob) / 2), adj = 0, cex = 0.8)
    text(bar_right + 0.05, bar_bottom, sprintf("%.4f", min_prob), adj = 0, cex = 0.8)
  }

  # Add title for color scale
  text((bar_left + bar_right) / 2, bar_top + 0.08, "Probability",
       cex = 0.95, font = 2, adj = 0.5)

  # Add note about size
  text(0.5, 0.05, "Size ∝ P^(1/3)", cex = 0.8, adj = 0.5)

  # ============================================================================
  # TITLE SPACING CONTROLS - Adjust these values to change title separation
  # ============================================================================
  # oma[3] controls overall top margin (higher = more space at top)
  # line controls how far down the main title sits (higher = lower on page)
  # To increase separation: increase oma[3] and/or increase line
  # To decrease separation: decrease oma[3] and/or decrease line
  # ============================================================================

  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 5.0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot.new()
  mtext(paste0("Exact Multinomial Outcome Probabilities at N = ", n, " (",
               nrow(exact_correct$proportions), " possible outcomes)"),
        side = 3, line = 3.5, cex = 1.0, font = 2)

  dev.off()

  cat("    Saved:", filename, "\n")
  cat("    Number of outcomes:", nrow(exact_correct$proportions), "\n")
  cat("    Probability range (correct):", sprintf("%.6f - %.6f",
                                                   min(exact_correct$probabilities),
                                                   max(exact_correct$probabilities)), "\n")
  cat("    Probability range (wrong):", sprintf("%.6f - %.6f",
                                                 min(exact_wrong$probabilities),
                                                 max(exact_wrong$probabilities)), "\n")
}

# Generate plots for each scenario and N value
for (scenario in scenarios) {
  cat("\nScenario:", scenario$name, "\n")
  cat("  ", scenario$correct_label, "vs", scenario$wrong_label, "\n")

  for (n in n_values) {
    create_probability_comparison(
      n,
      scenario$alpha_correct,
      scenario$alpha_wrong,
      scenario$correct_label,
      scenario$wrong_label,
      scenario$name
    )
  }
}

cat("\n=========================================\n")
cat("SUMMARY\n")
cat("=========================================\n")
cat("Generated exact probability plots for:\n")
cat("  Scenarios:", length(scenarios), "\n")
cat("  N values:", paste(n_values, collapse = ", "), "\n")
cat("  Total plots:", length(scenarios) * length(n_values), "\n")
cat("\nAll plots saved to:", output_dir, "\n\n")
cat("Done!\n")
