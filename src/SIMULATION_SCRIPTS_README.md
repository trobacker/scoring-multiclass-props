# Simulation Scripts Overview

This directory contains two main simulation scripts that explore scoring discrimination under different true forecast distributions.

## Scripts

### 1. `run_simulations_alpha_asymmetric.R`

**True distribution**: $\boldsymbol{\alpha}_{\text{true}} = (2, 3, 15)$
- **Type**: Asymmetric/skewed distribution
- **Mean proportions**: (0.1, 0.15, 0.75) — heavily skewed toward third category
- **Concentration**: $\sum \alpha_k = 20$ (moderate)
- **Use case**: Represents real-world scenarios like variant surveillance where one variant dominates

#### Scenarios Tested

| Scenario | Wrong Alpha | Type of Error | Description |
|----------|-------------|---------------|-------------|
| 1: Misspecified | (15, 3, 2) | Location (severe) | Completely reversed proportions |
| 2: Narrow/Calibrated | (200, 300, 1500) | Overconfident | 100× more concentrated, same mean |
| 3: Dispersed/Calibrated | (0.4, 0.6, 3) | Underconfident | 5× less concentrated, same mean |
| 4: Less Misspecified | (4, 5, 11) | Location (mild) | Moderately shifted proportions |

**Results location**: `./results/data/alpha_2_3_15/`

---

### 2. `run_simulations_alpha_symmetric.R` ⭐ NEW

**True distribution**: $\boldsymbol{\alpha}_{\text{true}} = (10, 10, 10)$
- **Type**: Symmetric/uniform distribution
- **Mean proportions**: (1/3, 1/3, 1/3) — equal prevalence across categories
- **Concentration**: $\sum \alpha_k = 30$ (moderate)
- **Use case**: Baseline case with no prior bias toward any category; useful for isolating dispersion effects

#### Scenarios Tested

| Scenario | Wrong Alpha | Type of Error | Description |
|----------|-------------|---------------|-------------|
| 1: Mislocation (symmetric) | (5, 10, 15) | Location | Shifted mean, same concentration (sum=30) |
| 2: Moderate Dispersed | (2, 2, 2) | Underconfident (5×) | 5× less concentrated (sum=6), same mean |
| 3: Extreme Dispersed | (0.5, 0.5, 0.5) | Underconfident (20×) | 20× less concentrated (sum=1.5), same mean |
| 4: Moderate Narrow | (50, 50, 50) | Overconfident (5×) | 5× more concentrated (sum=150), same mean |
| 5: Extreme Narrow | (200, 200, 200) | Overconfident (20×) | 20× more concentrated (sum=600), same mean |

**Results location**: `./results/data/alpha_10_10_10/`

---

## Key Differences Between Scripts

### Asymmetric Script (α = 2, 3, 15)
- **Focus**: Realistic skewed distributions
- **Scenarios**: Mix of location and dispersion errors
- **Dispersion range**: Single level each (5× dispersed, 100× narrow)
- **Best for**: Understanding real-world variant nowcasting

### Symmetric Script (α = 10, 10, 10)
- **Focus**: Balanced/uniform baseline
- **Scenarios**: Systematic exploration of dispersion
- **Dispersion range**: Two levels each (5× and 20× for both under/overconfidence)
- **Best for**: Isolating dispersion effects without location confounds

---

## Usage

### Running Individual Scripts

**From project root** (recommended):
```bash
# Asymmetric scenarios (original)
Rscript src/run_simulations_alpha_asymmetric.R

# Symmetric scenarios (new)
Rscript src/run_simulations_alpha_symmetric.R
```

**From src/ directory** (also supported):
```bash
cd src
Rscript run_simulations_alpha_asymmetric.R
Rscript run_simulations_alpha_symmetric.R
cd ..
```

**Important**: Both scripts automatically detect their execution context and source helper functions accordingly. They work from either the project root or the src/ directory.

### Running from RStudio or Positron

Both scripts automatically detect GUI environments (RStudio, Positron) and use safe settings to avoid forking issues on macOS:
- Scenarios run **sequentially** (one at a time)
- Uses **4 cores maximum** per scenario (conservative)
- More reliable, slightly slower

To run from RStudio/Positron Console:
```r
# Make sure working directory is project root
source("src/run_simulations_alpha_asymmetric.R")
# or
source("src/run_simulations_alpha_symmetric.R")
```

For maximum speed, run from terminal (see above).

**Important**: If you previously encountered errors with parallel processing in Positron/RStudio, the scripts have been updated to fix this. See `src/PARALLEL_PROCESSING_GUIDE.md` for details.

### Parameters

Both scripts use the same experimental parameters:
- **n_samp**: 100 (forecast ensemble size)
- **nsim**: 1000 (simulation iterations per scenario × n_counts combination)
- **n_counts_vec**: [1, 2, 3, 4, 5, 10, 25, 50, 100, 500] (sample sizes tested)
- **N_multinomial**: 100 (multinomial resamples for MN method)

---

## Output Structure

Each script generates organized output in `./results/`:

```
results/
├── data/
│   ├── alpha_2_3_15/          # Asymmetric results
│   │   ├── scenario1_misspecified_alpha_2_3_15.rds
│   │   ├── scenario2_narrow_calibrated_alpha_2_3_15.rds
│   │   ├── scenario3_dispersed_calibrated_alpha_2_3_15.rds
│   │   └── scenario4_less_misspecified_alpha_2_3_15.rds
│   └── alpha_10_10_10/        # Symmetric results (NEW)
│       ├── scenario1_mislocation_symmetric_alpha_10_10_10.rds
│       ├── scenario2_moderate_dispersed_alpha_10_10_10.rds
│       ├── scenario3_extreme_dispersed_alpha_10_10_10.rds
│       ├── scenario4_moderate_narrow_alpha_10_10_10.rds
│       └── scenario5_extreme_narrow_alpha_10_10_10.rds
└── plots/
    ├── alpha_2_3_15/          # Asymmetric plots
    └── alpha_10_10_10/        # Symmetric plots (NEW)
```

---

## Research Questions Addressed

### Asymmetric Script
1. How do methods perform when the true distribution is skewed?
2. Can scores discriminate severe vs. mild location errors?
3. How does extreme overconfidence (100×) affect discrimination?

### Symmetric Script
1. How do methods perform when no category is privileged?
2. What is the dose-response relationship between concentration error and discrimination?
3. Are moderate vs. extreme dispersion errors detected differently?
4. Is there symmetry in how over- and underconfidence are penalized?

---

## Computational Time

**Asymmetric script** (4 scenarios):
- Terminal (parallel scenarios): ~5-10 minutes
- RStudio (sequential scenarios): ~15-25 minutes

**Symmetric script** (5 scenarios):
- Terminal (parallel scenarios): ~7-12 minutes
- RStudio (sequential scenarios): ~20-30 minutes

*Times assume 10+ cores available for parallelization*

---

## Next Steps

After running both scripts, compare results to understand:
1. Whether discrimination patterns differ between symmetric vs. asymmetric true distributions
2. Whether dispersion effects scale monotonically (moderate → extreme)
3. Whether location errors interact with distribution symmetry

Use `results/plotting_code/combine_scenarios_small_n.R` and related scripts to create comparative visualizations.
