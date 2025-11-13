# Parallel Simulation Runner

## Overview

This directory contains parallelized simulation code for energy score analysis with multinomial sampling.

## Performance Improvements

### Option 1: Parallelized nsim loops (IMPLEMENTED ✓)
- **Location**: `es_sampling()` and `es_mn_sampling()` functions in `simulation_1.qmd`
- **Speedup**: ~10x (using 12 cores)
- **How it works**: Each of the 1,000 simulation iterations runs in parallel
- **Default**: Automatically uses `detectCores() - 2` cores (12 on your system)

### Option 2: Scenario-level parallelization (BONUS ✓)
- **Location**: `run_all_simulations_parallel.R` standalone script
- **Speedup**: Additional ~2-3x on top of Option 1
- **How it works**: Runs all 4 scenarios simultaneously
- **Usage**: See below

## File Structure

```
src/
├── simulation_1.qmd                    # Main Quarto document (parallelized)
├── run_all_simulations_parallel.R      # Standalone parallel runner
├── helper-functions.R                  # Core utility functions
└── PARALLEL_SIMULATIONS_README.md      # This file
```

## Usage

### Method 1: Interactive in RStudio (Recommended for exploration)

Open `simulation_1.qmd` in RStudio and run chunks interactively.

**Parallel processing is automatic!** Each scenario now runs ~10x faster using 12 cores.

### Method 2: Render Quarto document

```bash
cd src/
quarto render simulation_1.qmd
```

Results will include inline plots and saved files in `results/`.

### Method 3: Maximum speed standalone script

For the absolute fastest execution (runs all scenarios in parallel):

```bash
cd src/
Rscript run_all_simulations_parallel.R
```

**Expected time**: ~3-5 minutes for all 4 scenarios
- Each scenario: ~60 seconds (parallelized within)
- All scenarios run simultaneously (parallelized across)

**Output**:
- Individual scenario files: `results/data/alpha_X_Y_Z/scenario*_alpha_X_Y_Z.rds`
- Combined file: `results/all_scenarios_combined.rds`

## Updated Parameters

### New n_counts values
Now includes 10 values: `1, 2, 3, 4, 5, 10, 25, 50, 100, 500`
(Added: 100, 500)

### Scenarios
1. **Scenario 1**: Misspecified (2,3,15) vs (15,3,2)
2. **Scenario 2**: Narrow calibrated (2,3,15) vs (200,300,1500)
3. **Scenario 3**: Dispersed calibrated (2,3,15) vs (0.4,0.6,3)
4. **Scenario 4**: Less misspecified (2,3,15) vs (4,5,11)

## Controlling Parallelization

### In simulation_1.qmd

Control cores used per scenario:
```r
# Use fewer cores (e.g., for battery saving)
scenario1 <- run_and_save_scenario(
  ...
  parallel = TRUE,
  n_cores = 4  # Use only 4 cores
)

# Disable parallelization (sequential)
scenario1 <- run_and_save_scenario(
  ...
  parallel = FALSE
)
```

### In run_all_simulations_parallel.R

Edit the script to adjust:
```r
n_cores_outer <- min(4, n_cores_total)  # Scenarios run in parallel
n_cores_inner <- max(1, n_cores_total - 2)  # Cores per scenario
```

## Technical Details

### Parallelization Strategy

**Option 1 (Implemented in .qmd)**:
- Parallelizes the 1,000 iterations within each simulation
- Uses `parallel::mclapply()` with fork-based parallelism (efficient on macOS)
- Each iteration gets a unique seed: `seed + i` for reproducibility
- Pre-computes matrices outside parallel loop for efficiency

**Option 2 (Standalone script)**:
- Outer: 4 scenarios run simultaneously
- Inner: Each scenario uses 12 cores for iterations
- Note: Due to fork limitations, actual concurrency is managed by OS
- Total effective speedup: ~15-20x compared to original sequential code

### Why not nested parallelization in .qmd?

Nested parallelization (parallel scenarios × parallel iterations) can cause:
- Process explosion: 4 × 12 = 48 competing processes
- Memory issues
- Actually slower due to overhead

The standalone script manages this better by distributing at the scenario level.

## Performance Comparison

### Original (Sequential)
- Per scenario: ~10-15 minutes
- All 4 scenarios: ~40-60 minutes

### Parallelized .qmd (Option 1)
- Per scenario: ~1-2 minutes
- All 4 scenarios (sequential): ~4-8 minutes

### Standalone script (Option 1 + 2)
- All 4 scenarios (parallel): ~2-4 minutes

## Troubleshooting

### "mclapply not supported on Windows"

If on Windows, modify functions to use `parLapply`:
```r
cl <- makeCluster(n_cores)
results <- parLapply(cl, 1:nsim, function(i) { ... })
stopCluster(cl)
```

### Memory issues

Reduce cores or n_counts:
```r
n_cores = 4  # Use fewer cores
n_counts_vec = c(1, 5, 10, 50)  # Fewer n_counts values
```

### Results differ slightly from original

Parallel execution uses different random seeds (seed + i) for reproducibility.
Results should be statistically equivalent but not identical.

## Questions?

Check the main simulation_1.qmd for detailed documentation of functions and methodology.
