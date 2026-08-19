#!/usr/bin/env Rscript

# Quantization Effect: Exact Probabilities - Three Panel Version
# Shows dispersed, correct, and narrow forecasts side-by-side
# using exact Dirichlet-multinomial probabilities

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
cat("QUANTIZATION EXACT PROBABILITIES - TRIPANEL\n")
cat("=========================================\n\n")

# Parameters - extended to larger N values
n_values <- c(1, 2, 3, 4, 5, 10, 25, 50, 100)

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

# Function to create three-panel ternary comparison with exact probabilities
create_tripanel_probability_comparison <- function(n, alpha_dispersed, alpha_correct, alpha_narrow,
                                                   dispersed_label, correct_label, narrow_label,
                                                   scenario_name) {

  cat("  Creating exact probability tripanel for N =", n, "...\n")

  # Create Dirichlet density grids for background
  grid_dispersed <- create_dirichlet_grid(alpha_dispersed)
  grid_correct <- create_dirichlet_grid(alpha_correct)
  grid_narrow <- create_dirichlet_grid(alpha_narrow)

  # Compute exact probabilities
  exact_dispersed <- compute_exact_probabilities(n, alpha_dispersed)
  exact_correct <- compute_exact_probabilities(n, alpha_correct)
  exact_narrow <- compute_exact_probabilities(n, alpha_narrow)

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

  dispersed_colors <- get_prob_colors(exact_dispersed$probabilities)
  correct_colors <- get_prob_colors(exact_correct$probabilities)
  narrow_colors <- get_prob_colors(exact_narrow$probabilities)

  # Scale point sizes by probability (cube root for better visibility)
  # Adjust size range based on N to avoid overcrowding at large N
  scale_size <- function(probs, n) {
    if (n <= 10) {
      base_size <- 1.5
      max_size <- 4.0
    } else if (n <= 25) {
      base_size <- 1.0
      max_size <- 3.0
    } else if (n <= 50) {
      base_size <- 0.8
      max_size <- 2.5
    } else {  # n > 50
      base_size <- 0.5
      max_size <- 1.8
    }
    scaled <- (probs / max(probs))^(1/3)  # Cube root for better visual distribution
    base_size + scaled * (max_size - base_size)
  }

  dispersed_sizes <- scale_size(exact_dispersed$probabilities, n)
  correct_sizes <- scale_size(exact_correct$probabilities, n)
  narrow_sizes <- scale_size(exact_narrow$probabilities, n)

  # Adjust border width for large N
  border_width <- if (n > 50) 0.4 else if (n > 25) 0.5 else if (n > 10) 0.6 else 0.8

  # Create PNG with layout for legend
  filename <- paste0("exact_prob_", scenario_name, "_N", n, ".png")
  png(file.path(output_dir, filename), width = 2200, height = 700, res = 150)

  # Create layout: 3 ternary plots + 1 legend panel
  layout(matrix(c(1, 2, 3, 4), nrow = 1), widths = c(3, 3, 3, 1))

  # LEFT PANEL: Dispersed forecast
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
               cbind(grid_dispersed$a, grid_dispersed$b, grid_dispersed$c),
               pch = 15,
               col = grid_dispersed$color,
               cex = 0.5)

  # Plot exact probabilities (colored circles, size = probability)
  AddToTernary(graphics::points,
               exact_dispersed$proportions,
               pch = 21,  # Circle with border
               cex = dispersed_sizes,
               col = "black",  # Border color
               bg = dispersed_colors,  # Fill color
               lwd = border_width)

  title(main = "Dispersed", cex.main = 1.5, font.main = 2, line = -0.5)

  # MIDDLE PANEL: Correct forecast
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

  # Plot exact probabilities
  AddToTernary(graphics::points,
               exact_correct$proportions,
               pch = 21,
               cex = correct_sizes,
               col = "black",
               bg = correct_colors,
               lwd = border_width)

  title(main = "Correct", cex.main = 1.5, font.main = 2, line = -0.5)

  # RIGHT PANEL: Narrow forecast
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
               cbind(grid_narrow$a, grid_narrow$b, grid_narrow$c),
               pch = 15,
               col = grid_narrow$color,
               cex = 0.5)

  # Plot exact probabilities
  AddToTernary(graphics::points,
               exact_narrow$proportions,
               pch = 21,
               cex = narrow_sizes,
               col = "black",
               bg = narrow_colors,
               lwd = border_width)

  title(main = "Narrow", cex.main = 1.5, font.main = 2, line = -0.5)

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
  all_probs <- c(exact_dispersed$probabilities, exact_correct$probabilities, exact_narrow$probabilities)
  max_prob <- max(all_probs)
  min_prob <- min(all_probs)

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
  text(0.5, 0.05, "Size \u221d P^(1/3)", cex = 0.8, adj = 0.5)

  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 5.0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot.new()
  mtext(paste0("N = ", n),
        side = 3, line = 3.5, cex = 1.4, font = 2)

  dev.off()

  cat("    Saved:", filename, "\n")
  cat("    Number of outcomes:", nrow(exact_correct$proportions), "\n")
}

# Generate plots for each scenario and N value
for (scenario in tripanel_scenarios) {
  cat("\nScenario:", scenario$name, "\n")
  cat("  Three-panel:", scenario$dispersed_label, "|",
      scenario$correct_label, "|", scenario$narrow_label, "\n")

  for (n in n_values) {
    create_tripanel_probability_comparison(
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
cat("Generated exact probability tripanel plots for:\n")
cat("  Scenarios:", length(tripanel_scenarios), "\n")
cat("  N values:", paste(n_values, collapse = ", "), "\n")
cat("  Total plots:", length(tripanel_scenarios) * length(n_values), "\n")
cat("\nAll plots saved to:", output_dir, "\n\n")
cat("Done!\n")
