# Ternary Plot Saving Implementation

## Summary

Added functionality to automatically save ternary plots (visualizing true vs wrong alpha distributions) along with simulation results in both the Quarto document and standalone R script.

## Changes Made

### 1. `simulation_1.qmd`

#### New Function: `save_ternary_plot()` (Line ~328)

Added a new function to save base R graphics plots (ternary plots use base R graphics, not ggplot2):

```r
save_ternary_plot <- function(plot_func, scenario_name, alpha, alpha_wrong,
                              n_samp = 100,
                              results_dir = "../results/plots",
                              width = 800, height = 800) {
  # Create alpha-specific subdirectory
  alpha_str <- paste(alpha, collapse = "_")
  alpha_dir <- file.path(results_dir, paste0("alpha_", alpha_str))

  if (!dir.exists(alpha_dir)) {
    dir.create(alpha_dir, recursive = TRUE)
  }

  # Create filename with alpha (no timestamp - overwrites existing)
  filename <- file.path(alpha_dir, paste0(scenario_name, "_ternary_alpha_", alpha_str, ".png"))

  # Save using base R graphics
  png(filename, width = width, height = height, res = 150)
  plot_func(alpha, alpha_wrong, n_samp)
  dev.off()

  cat("Ternary plot saved to:", filename, "\n")
  return(filename)
}
```

**Key Features**:
- Uses `png()` and `dev.off()` instead of `ggsave()` (for base R graphics)
- Follows the same organization pattern: saves to `results/plots/alpha_X_Y_Z/` subdirectories
- Overwrites existing files (no timestamps)
- Filename pattern: `{scenario_name}_ternary_alpha_{alpha_values}.png`

#### Updated: `run_and_save_scenario()` (Line ~496)

Modified to save ternary plots when `save_results_flag = TRUE`:

```r
# Save if requested
if (save_results_flag) {
  save_results(df_results, scenario_name, alpha)
  save_plot(plots$comparison, "comparison", scenario_name, alpha)
  save_plot(plots$differences_method, "differences_method", scenario_name, alpha)
  save_ternary_plot(plot_alpha_ternary, scenario_name, alpha, alpha_wrong, n_samp = n_samp)  # NEW
}
```

### 2. `run_all_simulations_parallel.R`

#### Added: `plot_alpha_ternary()` Function (Line ~181)

Added the complete ternary plot visualization function (was not present in standalone script):

```r
plot_alpha_ternary <- function(alpha, alpha_wrong, n_samp = 100) {
  # Generate samples
  thetas <- lapply(FUN = draw_one_dirichlet, X = rep(1, n_samp), alpha = alpha)
  thetas_wrong <- lapply(FUN = draw_one_dirichlet, X = rep(1, n_samp), alpha = alpha_wrong)

  # Calculate mean and mode
  theta_true <- alpha / sum(alpha)
  theta_mode <- (alpha - 1) / (sum(alpha) - length(alpha))

  # Create ternary plot with density background
  # ... (full implementation)

  # Add legend
  legend("topright",
         legend = c("True alpha samples", "Wrong alpha samples", "True mean"),
         pch = c(20, 20, 21),
         col = c(rgb(0, 0, 1, 0.5), rgb(1, 0, 0, 0.5), "blue"),
         pt.bg = c(NA, NA, "lightblue"),
         pt.cex = c(1, 1, 1.5),
         bty = "n")
}
```

#### Added: `save_ternary_plot()` Function (Line ~224)

Same implementation as in `simulation_1.qmd`.

#### Updated: `run_scenario()` Function (Line ~312)

Modified to save ternary plots for each scenario:

```r
# Save results
save_results(df_results, scenario_config$name, scenario_config$alpha)

# Save ternary plot (NEW)
save_ternary_plot(plot_alpha_ternary, scenario_config$name,
                 scenario_config$alpha, scenario_config$alpha_wrong,
                 n_samp = n_samp)
```

## File Organization

Ternary plots are now saved to:

```
results/
└── plots/
    ├── alpha_2_3_15/
    │   ├── scenario1_misspecified_ternary_alpha_2_3_15.png
    │   ├── scenario2_narrow_calibrated_ternary_alpha_2_3_15.png
    │   ├── scenario3_dispersed_calibrated_ternary_alpha_2_3_15.png
    │   ├── scenario4_less_misspecified_ternary_alpha_2_3_15.png
    │   ├── scenario1_misspecified_comparison_alpha_2_3_15.png
    │   ├── scenario1_misspecified_differences_method_alpha_2_3_15.png
    │   └── ... (other plots)
    └── alpha_X_Y_Z/
        └── ... (plots for other alpha values)
```

## Usage Examples

### Interactive in RStudio (simulation_1.qmd)

When running scenarios, ternary plots are now automatically saved:

```r
scenario1 <- run_and_save_scenario(
  scenario_name = "scenario1_misspecified",
  alpha = c(2, 3, 15),
  alpha_wrong = c(15, 3, 2),
  n_samp = 100,
  nsim = 1000,
  base_seed = 42,
  save_results_flag = TRUE  # Ternary plot will be saved
)
```

**Output**:
```
Results saved to: ../results/data/alpha_2_3_15/scenario1_misspecified_alpha_2_3_15.rds
Plot saved to: ../results/plots/alpha_2_3_15/scenario1_misspecified_comparison_alpha_2_3_15.png
Plot saved to: ../results/plots/alpha_2_3_15/scenario1_misspecified_differences_method_alpha_2_3_15.png
Ternary plot saved to: ../results/plots/alpha_2_3_15/scenario1_misspecified_ternary_alpha_2_3_15.png
```

### Standalone Script

When running the parallel script, ternary plots are saved for each scenario automatically:

```bash
cd src/
Rscript run_all_simulations_parallel.R
```

**Output** (per scenario):
```
Running: scenario1_misspecified
Alpha (true): 2, 3, 15
Alpha (wrong): 15, 3, 2
========================================
  n_counts = 1 ... done
  n_counts = 2 ... done
  ...
Results saved to: ../results/data/alpha_2_3_15/scenario1_misspecified_alpha_2_3_15.rds
Ternary plot saved to: ../results/plots/alpha_2_3_15/scenario1_misspecified_ternary_alpha_2_3_15.png
Completed: scenario1_misspecified
```

## Technical Details

### Why Separate Save Function?

Ternary plots use **base R graphics** (`TernaryPlot()`, `AddToTernary()`, etc.), which differ from ggplot2:

- **ggplot2 plots**: Use `ggsave()` which works with plot objects
- **Base R graphics**: Require `png()` → plot commands → `dev.off()` workflow

The `save_ternary_plot()` function handles the base R graphics workflow by:
1. Opening a PNG device with `png(filename, ...)`
2. Calling the plot function to draw to the device
3. Closing the device with `dev.off()`

### Plot Specifications

- **Resolution**: 150 DPI
- **Dimensions**: 800 × 800 pixels
- **Format**: PNG
- **Overwrite behavior**: Existing files with the same name are overwritten (consistent with other results)

## Benefits

1. **Complete archival**: All visualization outputs are saved, not just summary statistics
2. **Reproducibility**: Can review alpha distributions without re-running simulations
3. **Documentation**: Visual documentation of each scenario's parameter choices
4. **Consistency**: Follows same organizational structure as other results (alpha subdirectories)
5. **Publication-ready**: High-resolution plots ready for inclusion in papers/presentations

## Notes

- Ternary plots are generated with the same `n_samp` parameter used for the simulation (default: 100)
- The visualization shows:
  - **Blue points**: Samples from true alpha distribution
  - **Red points**: Samples from wrong alpha distribution
  - **Blue circle**: Mean of true alpha distribution
  - **Background**: Density heatmap of true alpha distribution
- Files are saved to the same alpha-specific subdirectories as other results for easy organization

## Verification

To verify ternary plots are being saved, check the results directory:

```r
# List all ternary plots for a specific alpha
list.files("../results/plots/alpha_2_3_15", pattern = "*ternary*.png")

# Expected output:
# [1] "scenario1_misspecified_ternary_alpha_2_3_15.png"
# [2] "scenario2_narrow_calibrated_ternary_alpha_2_3_15.png"
# [3] "scenario3_dispersed_calibrated_ternary_alpha_2_3_15.png"
# [4] "scenario4_less_misspecified_ternary_alpha_2_3_15.png"
```
