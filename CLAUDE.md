# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This repository contains research code for exploring scoring methods for multi-class proportions, inspired by the Variant Nowcast Hub. The work examines energy scores for evaluating probabilistic forecasts of variant proportions using Dirichlet distributions, multinomial sampling, and Dirichlet-multinomial models.

## Core Research Question

The project investigates whether scoring on counts versus proportions makes a difference, and explores the role of multinomial sampling in the scoring procedure for multi-class proportion forecasts. It specifically examines how energy scores behave under different sampling schemes and sample sizes.

## Key Files and Architecture

### Source Files (`src/`)

- **`helper-functions.R`**: Core utility functions for the entire project
  - Dirichlet density functions (`dd_func`)
  - Multinomial density functions (`mn_func`)
  - Dirichlet sampling (`draw_one_dirichlet`)
  - Sample space enumeration (`enumerate_multinomial_sample_space`)
  - Ternary plotting functions for visualizing distributions and mixtures:
    - `plot_implied_multinomial`: Visualizes single multinomial distributions
    - `plot_implied_multinomial_mixture`: Visualizes mixtures of multinomial distributions
    - `plot_dirichlet_multinomial`: Visualizes Dirichlet-multinomial distributions
  - Color utilities (`col_r`) for transparent plotting

- **`simulation_1.qmd`**: Main simulation experiments (Quarto document)
  - Defines two key simulation functions:
    - `es_sampling`: Scores WITHOUT multinomial sampling
    - `es_mn_sampling`: Scores WITH multinomial sampling (100 multinomial samples per theta)
  - Compares "correct" forecasts (using true alpha) vs "incorrect" forecasts (using wrong alpha)
  - Varies `n_counts` parameter (5, 25, 50, 100, 500) to examine effect of sample size
  - Contains plotting support via `plot_sampling_results`
  - Fixed parameters: 100 submitted samples (thetas), 1000 simulation runs

- **`scoring_counts_and_props.qmd`**: Theoretical exploration and proof
  - Demonstrates that energy scores on counts = N × (energy scores on proportions)
  - Mathematical proof that scoring methods are equivalent up to a scaling factor
  - Explores why multinomial sampling is used despite this equivalence
  - Shows that multinomial sampling captures variation dependent on N without requiring modelers to forecast N

- **`variant_scoring.Rmd`**: Original exploratory document (R Markdown)
  - Contains initial conceptual introduction to variant nowcast scoring challenges
  - Visualizes implied multinomial distributions for different N values (5, 50, 500)
  - Compares mixture of implied multinomials vs true Dirichlet-multinomial
  - Contains early simulation code that was later refined in `simulation_1.qmd`

### Output Files

- `simulation_1.html`: Rendered output from `simulation_1.qmd`
- `simulation_1_files/`: Supporting files (plots, libraries) for rendered HTML

## Key Concepts and Notation

- **θ (theta)**: K-vector of true proportions for K classes/variants
- **C**: Vector of observed counts for each class
- **N**: Total number of sequences/observations (N = Σ Cₖ)
- **α (alpha)**: Dirichlet distribution parameters
- **Energy Score (ES)**: Proper scoring rule for multivariate probabilistic forecasts
- Models submit 100 samples from predictive distribution of θ
- Hub generates predictions for counts by sampling from Multinomial(N, θ̂⁽ˢ⁾)

## Development Workflow

### Running Quarto Documents

```bash
# Render a single Quarto document to HTML
quarto render src/simulation_1.qmd

# Render all Quarto documents
quarto render src/

# Preview with live reload during development
quarto preview src/simulation_1.qmd
```

### Running R Markdown Documents

Open in RStudio and use the "Knit" button, or from R console:
```r
rmarkdown::render("src/variant_scoring.Rmd")
```

### Working in R

Source the helper functions at the start of any R session:
```r
source("src/helper-functions.R")
```

## Key R Packages

The project depends on:
- `scoringRules`: Energy score and variogram score calculations
- `brms`: Dirichlet distribution functions (`rdirichlet`, `ddirichlet`)
- `Ternary`: Ternary plots for visualizing 3-class proportions
- `mc2d`: Multinomial distribution functions (`dmultinomial`, `rmultinomial`)
- `extraDistr`: Dirichlet-multinomial distribution (`ddirmnom`, `rdirmnom`)
- `ggplot2`, `ggsimplex`: Alternative visualization methods
- `tidyr`: Data reshaping for comparative plots

## Common Patterns

### Running a Simulation Experiment

Simulations follow this general pattern:
1. Define true alpha parameters (e.g., `alpha <- c(2, 4, 4)`)
2. Define wrong/comparison alpha (e.g., `alpha_wrong <- c(20, 5, 5)`)
3. Pre-generate forecast samples from Dirichlet distributions
4. Loop over simulation runs (typically 1000):
   - Generate observation from true distribution
   - Score both correct and incorrect forecasts
   - Store energy scores
5. Compare distributions of energy scores via histograms and pairwise differences

### Energy Score Calculation

```r
# dat should be a matrix with K rows and S columns (K classes, S samples)
thetas_mat <- simplify2array(thetas_list)

# Calculate energy score for observation y
es <- es_sample(y = observed_proportions, dat = thetas_mat)

# For counts: ES_counts = N * ES_proportions
es_counts <- n * es_sample(y = proportions, dat = thetas_mat)
```

## Repository State

The codebase has uncommitted changes:
- Modified: `src/simulation_1.html`, `src/simulation_1.qmd`, `src/variant_scoring.Rmd`
- Untracked: `scoring-multiclass-props.Rproj`, `src/simulation_1_files/`

The main development branch is `main`.
