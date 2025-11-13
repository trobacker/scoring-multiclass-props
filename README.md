# Scoring Multi-Class Proportions

Research code for exploring scoring methods for multi-class proportions, inspired by the [Variant Nowcast Hub](https://github.com/reichlab/variant-nowcast-hub). This work examines energy scores for evaluating probabilistic forecasts of variant proportions using Dirichlet distributions, multinomial sampling, and Dirichlet-multinomial models.

## Research Questions

This project investigates:
- **Does scoring on counts vs. proportions make a difference?** We show mathematically that energy scores on counts equal N × energy scores on proportions.
- **What is the role of multinomial sampling in the scoring procedure?** We explore how multinomial sampling captures variation dependent on sample size N without requiring modelers to forecast N explicitly.
- **How do energy scores behave under different sampling schemes and sample sizes?** We examine performance across varying N (1 to 500 sequences) and different Dirichlet parameter configurations.

## Key Features

- **Comprehensive Scenarios**: Examines misspecified, calibrated (narrow/dispersed), and moderate misspecification cases
- **Multiple Dimensionalities**: 3-class and 5-class proportion experiments
- **Rich Visualizations**: Ternary plots, dotplots comparing scoring methods, and distribution comparisons
- **Reproducible Results**: All simulation results, plots, and data files are version-controlled
- **Parallel Processing**: Simulation experiments support parallel execution for efficiency

## Repository Structure

```
.
├── src/
│   ├── simulation_1.qmd              # Main simulation experiments (Quarto)
│   ├── helper-functions.R            # Core utility functions
│   ├── variant_scoring.Rmd           # Original exploratory analysis
│   ├── scoring_counts_and_props.qmd  # Theoretical proof: counts vs props
│   └── run_all_simulations_parallel.R # Standalone parallel runner
├── results/
│   ├── data/                         # Simulation results (.rds files)
│   │   ├── alpha_2_3_15/            # 3-class scenarios
│   │   ├── alpha_2_5_10_20_40/      # 5-class scenarios
│   │   └── alpha_10_10_10_10_10/    # Uniform 5-class scenario
│   ├── plots/                        # Visualizations (.png files)
│   ├── all_scenarios_combined.rds    # Combined 3-class results
│   └── all_scenarios_highdim_combined.rds  # Combined 5-class results
├── docs/
│   ├── CHANGELOG.md                         # Project changelog
│   ├── PARALLELIZATION_IMPLEMENTATION_NOTES.md
│   ├── PARALLEL_SIMULATIONS_README.md
│   └── TERNARY_PLOT_SAVING_NOTES.md
└── README.md                         # This file
```

## Requirements

### R Packages

```r
install.packages(c(
  "scoringRules",    # Energy score calculations
  "brms",            # Dirichlet distribution functions
  "Ternary",         # Ternary plots
  "mc2d",            # Multinomial distribution
  "extraDistr",      # Dirichlet-multinomial
  "ggplot2",         # Visualization
  "tidyr",           # Data manipulation
  "parallel"         # Parallel processing
))
```

### System Requirements

- R >= 4.0
- For parallel processing: macOS/Linux (uses `mclapply`; Windows users need `parLapply`)
- Quarto (optional, for rendering .qmd documents)

## Getting Started

### 1. Clone the Repository

```bash
git clone https://github.com/trobacker/scoring-multiclass-props.git
cd scoring-multiclass-props
```

### 2. Run Interactive Simulations

Open `src/simulation_1.qmd` in RStudio and run chunks interactively.

### 3. Render Quarto Documents

```bash
cd src/
quarto render simulation_1.qmd
```

### 4. Run All Simulations with Parallel Processing

To run all scenarios with parallelization:

```bash
cd src/
Rscript run_all_simulations_parallel.R
```

## Key Simulation Functions

### `es_sampling()`
Scores **without** multinomial sampling - samples from Dirichlet, generates multinomial observation, scores directly on proportions.

### `es_mn_sampling()`
Scores **with** multinomial sampling - for each Dirichlet sample, draws 100 multinomial samples before scoring.

### `run_and_save_scenario()`
Complete scenario runner that:
- Pre-generates forecast samples
- Runs simulations across n_counts values (1, 2, 3, 4, 5, 10, 25, 50, 100, 500)
- Generates comparison plots and ternary visualizations
- Saves results to `results/data/` and `results/plots/`

## Simulation Scenarios

### 3-Class Scenarios (alpha_2_3_15)

1. **Scenario 1 - Misspecified**: True (2,3,15) vs Wrong (15,3,2)
2. **Scenario 2 - Narrow Calibrated**: True (2,3,15) vs Wrong (200,300,1500)
3. **Scenario 3 - Dispersed Calibrated**: True (2,3,15) vs Wrong (0.4,0.6,3)
4. **Scenario 4 - Less Misspecified**: True (2,3,15) vs Wrong (4,5,11)

### 5-Class Scenarios

- **alpha_2_5_10_20_40**: Ascending proportions experiments
- **alpha_10_10_10_10_10**: Uniform vs peaked distributions

## Key Findings

1. **Energy scores on counts = N × energy scores on proportions** (proven mathematically in `scoring_counts_and_props.qmd`)

2. **Multinomial sampling captures N-dependent variation** without requiring modelers to forecast N, which is valuable for variant nowcasting applications

3. **Energy scores distinguish well-calibrated vs. misspecified forecasts** across all sample sizes examined

## Documentation

- **[CHANGELOG.md](docs/CHANGELOG.md)**: Project history and recent changes
- **[PARALLELIZATION_IMPLEMENTATION_NOTES.md](docs/PARALLELIZATION_IMPLEMENTATION_NOTES.md)**: Technical details on parallel processing
- **[PARALLEL_SIMULATIONS_README.md](docs/PARALLEL_SIMULATIONS_README.md)**: User guide for running parallel simulations
- **[TERNARY_PLOT_SAVING_NOTES.md](docs/TERNARY_PLOT_SAVING_NOTES.md)**: Ternary plot implementation details

## Citation

If you use this code in your research, please cite:

```bibtex
@misc{scoring-multiclass-props,
  author = {Robacker, Thomas},
  title = {Scoring Multi-Class Proportions},
  year = {2024},
  url = {https://github.com/trobacker/scoring-multiclass-props}
}
```

## Acknowledgments

This research was inspired by the [Variant Nowcast Hub](https://github.com/reichlab/variant-nowcast-hub) and explores scoring methodology for probabilistic variant proportion forecasts.

## License

MIT License - See LICENSE file for details.

## Contact

For questions or issues, please open an issue on GitHub.
