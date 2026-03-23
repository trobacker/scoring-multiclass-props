# Additional Extreme Scenarios for Alpha (10, 10, 10)

## Overview

This document describes two **ultra-extreme** scenarios added to extend the dose-response investigation of dispersion errors.

## Existing Scenarios (Previously Completed)

The main symmetric script (`run_simulations_alpha_symmetric.R`) includes 5 scenarios:

| Scenario | Alpha Wrong | Concentration Ratio | Type |
|----------|-------------|-------------------|------|
| 1: Mislocation | (5, 10, 15) | 1× (same) | Location error |
| 2: Moderate Dispersed | (2, 2, 2) | 5× less | Underconfident |
| 3: Extreme Dispersed | (0.5, 0.5, 0.5) | 20× less | Very underconfident |
| 4: Moderate Narrow | (50, 50, 50) | 5× more | Overconfident |
| 5: Extreme Narrow | (200, 200, 200) | 20× more | Very overconfident |

All scenarios use **true alpha = (10, 10, 10)** with sum = 30.

---

## New Scenarios (Additional Script)

The new script (`run_simulations_alpha_symmetric_additional.R`) adds 2 **ultra-extreme** scenarios:

### Scenario 6: Ultra-Dispersed (Extremely Underconfident)

**Alpha wrong**: (0.1, 0.1, 0.1)
- **Sum**: 0.3 (vs. true sum = 30)
- **Concentration ratio**: **100× less concentrated**
- **Mean**: Still (1/3, 1/3, 1/3) - same as true
- **Interpretation**: Forecast is so underconfident it's nearly uniform across the entire simplex
- **Expected behavior**:
  - Non-MN method will likely favor this forecast even more strongly at small n
  - MN method should still detect it as wrong, but discrimination may be weaker than scenario 3

### Scenario 7: Ultra-Narrow (Extremely Overconfident)

**Alpha wrong**: (1000, 1000, 1000)
- **Sum**: 3000 (vs. true sum = 30)
- **Concentration ratio**: **100× more concentrated**
- **Mean**: Still (1/3, 1/3, 1/3) - same as true
- **Interpretation**: Forecast is so overconfident samples cluster extremely tightly around the mean
- **Expected behavior**:
  - Non-MN method will strongly penalize this (perfect discrimination at small n)
  - MN method will show even worse discrimination than scenario 5 (noise dominates signal)

---

## Research Questions

These additional scenarios help answer:

1. **Is there a dose-response relationship?**
   - Does discrimination increase/decrease monotonically with concentration error?
   - Scenarios 2 → 3 → 6 (dispersed): 5× → 20× → 100×
   - Scenarios 4 → 5 → 7 (narrow): 5× → 20× → 100×

2. **Are there saturation effects?**
   - Does discrimination plateau at extreme misspecification?
   - E.g., if scenario 3 shows 100% discrimination, will scenario 6 also show 100%?

3. **Does the asymmetry persist at extreme levels?**
   - Non-MN vs. MN performance gap may widen or narrow at ultra-extreme levels
   - Tests whether the fundamental trade-off scales

4. **What are the practical limits of each method?**
   - At what level of misspecification does each method completely fail/succeed?

---

## Comparison to Previous Scenarios

### Dispersion Gradient (Underconfident)

| Scenario | Alpha Wrong | Sum | Ratio | Non-MN Expected | MN Expected |
|----------|-------------|-----|-------|-----------------|-------------|
| 2: Moderate | (2, 2, 2) | 6 | 5× less | May favor wrong at n<5 | Correct discrimination |
| 3: Extreme | (0.5, 0.5, 0.5) | 1.5 | 20× less | Strongly favors wrong at n<10 | Good discrimination |
| **6: Ultra** | **(0.1, 0.1, 0.1)** | **0.3** | **100× less** | **Extreme bias to wrong?** | **Weakening discrimination?** |

### Concentration Gradient (Overconfident)

| Scenario | Alpha Wrong | Sum | Ratio | Non-MN Expected | MN Expected |
|----------|-------------|-----|-------|-----------------|-------------|
| 4: Moderate | (50, 50, 50) | 150 | 5× more | Good discrimination | Some noise interference |
| 5: Extreme | (200, 200, 200) | 600 | 20× more | Excellent discrimination | Substantial noise |
| **7: Ultra** | **(1000, 1000, 1000)** | **3000** | **100× more** | **Perfect discrimination?** | **Severe noise?** |

---

## Files and Output

### Script
- **Filename**: `src/run_simulations_alpha_symmetric_additional.R`
- **Based on**: `src/run_simulations_alpha_symmetric.R`
- **Scenarios**: Only scenarios 6 and 7 (not re-running scenarios 1-5)

### Output Location
All results save to the same directory as the original symmetric scenarios:
```
./results/data/alpha_10_10_10/
├── scenario1_mislocation_symmetric_alpha_10_10_10.rds
├── scenario2_moderate_dispersed_alpha_10_10_10.rds
├── scenario3_extreme_dispersed_alpha_10_10_10.rds
├── scenario4_moderate_narrow_alpha_10_10_10.rds
├── scenario5_extreme_narrow_alpha_10_10_10.rds
├── scenario6_ultra_dispersed_alpha_10_10_10.rds          ← NEW
├── scenario7_ultra_narrow_alpha_10_10_10.rds             ← NEW
└── additional_scenarios_combined.rds                      ← NEW
```

### Plots
Ternary plots and dotplots save to:
```
./results/plots/alpha_10_10_10/
```

---

## Expected Runtime

**2 scenarios × 10 sample sizes × 1000 iterations × 2 methods = 40,000 simulations**

| Environment | Strategy | Est. Time |
|-------------|----------|-----------|
| **Terminal** | Parallel (full cores) | ~3-5 minutes |
| **Positron/RStudio** | Sequential + 4 cores | ~8-12 minutes |

Significantly faster than the full 5-scenario script since only running 2 scenarios.

---

## Analysis Plan

After completion, compare across all 7 scenarios:

1. **Plot discrimination vs. concentration ratio**
   - X-axis: log(concentration ratio) [-2, -1.3, -0.7, 0, 0.7, 1.3, 2]
   - Y-axis: Proportion correct wins
   - Separate lines for each n_counts value
   - Separate panels for MN vs non-MN

2. **Identify saturation points**
   - Where does discrimination reach 100% (or 0%)?
   - Does doubling concentration error double discrimination difference?

3. **Test for monotonicity**
   - Does discrimination always increase with |concentration error|?
   - Or are there non-monotonic regions?

4. **Method comparison**
   - At what concentration error magnitude do the methods show maximum divergence?
   - Is there a crossover point where MN becomes worse than non-MN?

---

## Usage

### Running the Additional Scenarios

```bash
# From project root (fastest)
Rscript src/run_simulations_alpha_symmetric_additional.R

# From src/ directory
cd src
Rscript run_simulations_alpha_symmetric_additional.R
```

### Re-running All Scenarios

To regenerate all 7 scenarios from scratch:
```bash
# First run scenarios 1-5
Rscript src/run_simulations_alpha_symmetric.R

# Then run scenarios 6-7
Rscript src/run_simulations_alpha_symmetric_additional.R
```

---

## Next Steps

1. ✅ Run additional scenarios (scenarios 6-7)
2. Load all 7 scenario results
3. Create combined analysis plots
4. Compare dose-response patterns
5. Update `interpretation_patterns.md` with findings
6. Consider whether even more extreme scenarios are needed (diminishing returns likely)

---

## Notes

- The ultra-dispersed scenario (α = 0.1, 0.1, 0.1) is approaching the lower limit of meaningful Dirichlet concentration. Values much lower may have numerical issues.
- The ultra-narrow scenario (α = 1000, 1000, 1000) represents extreme overconfidence where all forecast samples will be nearly identical to the mean.
- These scenarios test the absolute limits of each scoring method's sensitivity to dispersion errors.
