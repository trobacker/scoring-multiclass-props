# Plotting Scripts Overview

This directory contains comprehensive plotting scripts for analyzing simulation results from both asymmetric and symmetric alpha scenarios.

## Scripts

### 1. `plot_asymmetric_alpha_all_n.R`

Analyzes scenarios with **asymmetric alpha = (2, 3, 15)** - skewed distribution.

**Scenarios included**: 4
- Scenario 1: Misspecified location
- Scenario 2: Narrow/Calibrated (overconfident)
- Scenario 3: Dispersed/Calibrated (underconfident)
- Scenario 4: Less Misspecified location

**Outputs**:
- `asymmetric_alpha_energy_scores_vs_n.png` - Mean ES for correct vs incorrect forecasts
- `asymmetric_alpha_discrimination_vs_n.png` - Win percentage across n_counts
- `asymmetric_alpha_es_differences_vs_n.png` - Mean ES differences
- `asymmetric_alpha_summary_stats.csv` - Complete numerical summary

**Total simulations analyzed**: 80,000

---

### 2. `plot_symmetric_alpha_all_n.R`

Analyzes scenarios with **symmetric alpha = (10, 10, 10)** - uniform distribution.

**Scenarios included**: 7
- S1: Mislocation (location error)
- S2: Moderate Dispersed (5× underconfident)
- S3: Extreme Dispersed (20× underconfident)
- S4: Moderate Narrow (5× overconfident)
- S5: Extreme Narrow (20× overconfident)
- S6: Ultra Dispersed (100× underconfident) ✨ NEW
- S7: Ultra Narrow (100× overconfident) ✨ NEW

**Outputs**:
- `symmetric_alpha_energy_scores_vs_n.png` - Mean ES for correct vs incorrect forecasts
- `symmetric_alpha_discrimination_vs_n.png` - Win percentage across n_counts
- `symmetric_alpha_es_differences_vs_n.png` - Mean ES differences
- `symmetric_alpha_dose_response.png` - Discrimination vs concentration error magnitude
- `symmetric_alpha_summary_stats.csv` - Complete numerical summary

**Total simulations analyzed**: 140,000

**Special feature**: Includes dose-response plot showing how discrimination changes with 5×, 20×, and 100× concentration errors.

---

### 3. `combine_scenarios_small_n.R` (Existing)

Focuses on small sample sizes (n_counts = 1-4) for asymmetric alpha scenarios.

**Outputs**:
- `combined_scenarios_small_n_comparison.png` - Dotplot comparison
- `combined_scenarios_small_n_differences.png` - Differences with mean overlays

---

## Key Visualizations

### Energy Scores vs n_counts
Shows how mean energy scores change with sample size for both correct and incorrect forecasts.
- **Lower is better** (better forecasts have lower ES)
- Faceted by method (Multinomial vs Non-Multinomial)
- Different colors for different scenarios
- Solid lines = correct forecast, dashed lines = incorrect forecast

### Discrimination (Win %) vs n_counts
Shows the proportion of simulations where the correct forecast scored better (lower ES) than the incorrect forecast.
- **Higher is better** (100% = perfect discrimination)
- Horizontal dotted line at 50% = random chance
- Solid lines = Multinomial, dashed lines = Non-Multinomial
- Log scale x-axis shows all n_counts from 1 to 500

### Energy Score Differences vs n_counts
Shows mean(ES_correct - ES_incorrect):
- **Negative values** = correct forecast scores better (desired)
- **Positive values** = incorrect forecast scores better (problematic)
- Horizontal line at 0 = no discrimination
- Same styling as discrimination plot

### Dose-Response Plot (Symmetric only)
Shows how discrimination changes with increasing concentration error magnitude:
- X-axis: Concentration ratio (5×, 20×, 100×) on log scale
- Y-axis: Proportion correct wins
- Faceted by error type (Underconfident vs Overconfident)
- Different colors for different n_counts values
- Reveals whether discrimination scales monotonically

---

## Usage

### Running the Scripts

```bash
# From project root
Rscript results/plotting_code/plot_asymmetric_alpha_all_n.R
Rscript results/plotting_code/plot_symmetric_alpha_all_n.R

# Or from plotting_code directory
cd results/plotting_code
Rscript plot_asymmetric_alpha_all_n.R
Rscript plot_symmetric_alpha_all_n.R
```

### Viewing Results

All plots are saved as high-resolution PNG files (300 dpi) in `./results/plotting_code/`:

```bash
ls -la results/plotting_code/*.png
```

### Summary Statistics

CSV files contain detailed statistics for all scenario × method × n_counts combinations:
- Mean energy scores (correct and incorrect)
- Mean ES differences
- Proportion of correct wins (discrimination)

Load in R:
```r
asym_stats <- read.csv("results/plotting_code/asymmetric_alpha_summary_stats.csv")
sym_stats <- read.csv("results/plotting_code/symmetric_alpha_summary_stats.csv")
```

---

## Color Schemes

### Asymmetric Alpha (4 scenarios)
- Red: Scenario 1 (Misspecified)
- Blue: Scenario 2 (Narrow/Calibrated)
- Green: Scenario 3 (Dispersed/Calibrated)
- Purple: Scenario 4 (Less Misspecified)

### Symmetric Alpha (7 scenarios)
- Purple: S1 (Mislocation)
- Light Green → Green → Dark Green: S2, S3, S6 (Dispersed 5×, 20×, 100×)
- Light Orange → Orange → Red: S4, S5, S7 (Narrow 5×, 20×, 100×)

**Rationale**: Color gradients show increasing severity of dispersion errors.

---

## Key Patterns to Look For

### In Discrimination Plots

1. **Scenario 1 (Misspecified location)**: Both methods show high discrimination (~90-100%)
2. **Scenario 2 (Overconfident)**: Non-MN shows 100% discrimination at small n, MN struggles
3. **Scenario 3 (Underconfident)**: MN shows good discrimination, Non-MN fails (0% at n=1-2)
4. **Convergence**: All methods converge to similar discrimination at n ≥ 100

### In Energy Score Plots

1. **Separation**: Large vertical gap between correct/incorrect = good discrimination
2. **Method differences**: Multinomial scores are often smoother across n
3. **Small n behavior**: Most dramatic differences occur at n < 10

### In Dose-Response Plots (Symmetric only)

1. **Monotonicity**: Does discrimination increase monotonically with concentration error?
2. **Saturation**: Do extreme errors (100×) show plateau effects?
3. **Asymmetry**: Compare underconfident vs overconfident patterns

---

## Troubleshooting

### Script fails to find data files
- Check that simulations have completed
- Verify files exist in `./results/data/alpha_2_3_15/` or `./results/data/alpha_10_10_10/`
- Ensure working directory is project root

### Missing scenarios in symmetric plot
- Check that all 7 scenario files exist
- Scenarios 6-7 require running `run_simulations_alpha_symmetric_additional.R`
- Script will warn about missing files but plot available data

### Plots look crowded
- 7 scenarios create busy plots - this is expected
- Use dose-response plot to focus on dispersion scenarios only
- Consider creating custom subsets for presentations

---

## Next Steps

1. **Examine plots**: Look for patterns in discrimination across scenarios
2. **Compare alpha settings**: How do symmetric vs asymmetric results differ?
3. **Identify critical n**: At what n_counts do methods diverge/converge?
4. **Dose-response analysis**: Is discrimination linear in log(concentration error)?
5. **Write up findings**: Use plots and summary stats for manuscript

---

## Summary Statistics Format

CSV files contain these columns:
- `scenario_label`: Scenario name
- `method_label`: "Multinomial" or "Non-Multinomial"
- `n_counts`: Sample size (1, 2, 3, 4, 5, 10, 25, 50, 100, 500)
- `mean_es_correct`: Mean energy score for correct forecast
- `mean_es_incorrect`: Mean energy score for incorrect forecast
- `mean_es_diff`: Mean(ES_correct - ES_incorrect)
- `prop_correct_wins`: Proportion where ES_correct < ES_incorrect

Each row represents 1,000 simulation iterations summarized.

---

## Plot Specifications

All plots use:
- **Theme**: `theme_bw()` for clean, publication-ready appearance
- **Resolution**: 300 DPI (high quality for papers)
- **Sizes**: 10-12 inches wide, 6-8 inches tall
- **Font sizes**: Base 11pt, titles 13pt bold
- **Log scale x-axis**: Shows all n_counts values clearly
- **Color-blind friendly**: Uses ColorBrewer palettes where possible

---

## Advanced Usage

### Customizing Plots

To modify plots, edit the R scripts:

1. **Change colors**: Update `scale_color_manual()` values
2. **Add/remove scenarios**: Modify scenario lists at top of scripts
3. **Filter n_counts**: Add `filter()` before plotting
4. **Change plot dimensions**: Modify `width` and `height` in `ggsave()`

### Combining Results

To analyze both alpha settings together:

```r
asym <- read.csv("results/plotting_code/asymmetric_alpha_summary_stats.csv")
sym <- read.csv("results/plotting_code/symmetric_alpha_summary_stats.csv")

asym$alpha_type <- "Asymmetric"
sym$alpha_type <- "Symmetric"

combined <- rbind(asym, sym)

# Now plot with facet_wrap(~ alpha_type)
```

---

## References

These plotting scripts complement:
- `writing/brainstorming/simulations/interpretation_patterns.md` - Detailed findings
- `writing/brainstorming/simulations/results_summary.md` - Publication summary
- `src/ADDITIONAL_SCENARIOS_README.md` - Scenario descriptions
