# Changelog

All notable changes to the scoring-multiclass-props project are documented here.

## 2025-11-13

### Added

- **Parallel Processing Implementation**: Implemented two-level parallelization strategy achieving 10-20x speedup
  - Added parallel processing to `es_sampling()` and `es_mn_sampling()` functions using `parallel::mclapply()`
  - Functions now use `detectCores() - 2` cores by default with configurable parallel parameters
  - Unique seed management (`seed + i`) ensures reproducibility across parallel executions
  - Sequential fallback mode available via `parallel = FALSE` parameter

- **Extended Sample Size Range**: Added larger n_counts values (100, 500) to examine asymptotic behavior
  - Updated default `n_counts_vec` from `c(1, 2, 3, 4, 5, 10, 25, 50)` to `c(1, 2, 3, 4, 5, 10, 25, 50, 100, 500)`

- **Ternary Plot Auto-Saving**: Implemented automatic saving of ternary plot visualizations
  - Created `save_ternary_plot()` function for base R graphics plots
  - Ternary plots now saved alongside comparison and boxplot results
  - Files organized in alpha-specific subdirectories: `results/plots/alpha_X_Y_Z/`
  - Filename pattern: `{scenario_name}_ternary_alpha_{alpha_values}.png`

- **Comprehensive Documentation**: Created three implementation notes documents
  - `src/PARALLELIZATION_IMPLEMENTATION_NOTES.md`: Detailed technical documentation of parallelization changes
  - `src/PARALLEL_SIMULATIONS_README.md`: User-facing guide for running parallel simulations
  - `src/TERNARY_PLOT_SAVING_NOTES.md`: Documentation of ternary plot saving implementation

- **Simulation Artifacts**: Generated simulation output files
  - `src/simulation_1.rmarkdown`: Intermediate R markdown rendering artifact
  - `src/simulation_1_files/`: Supporting files for rendered HTML output (figures, libraries)

### Changed

- **Function Signatures**: Updated simulation functions with parallel processing parameters
  - `es_sampling()`: Added `parallel` and `n_cores` parameters
  - `es_mn_sampling()`: Added `parallel` and `n_cores` parameters
  - `run_simulation_experiment()`: Added parallel parameter pass-through
  - `run_and_save_scenario()`: Added parallel configuration support

- **Results Organization**: Updated plot saving to include ternary visualizations
  - Modified `run_and_save_scenario()` to call `save_ternary_plot()`
  - Ternary plots saved with 150 DPI resolution at 800×800 pixels

### Performance

- **Sequential Execution** (original): ~40-60 minutes for all 4 scenarios
- **Parallelized Functions** (new): ~4-8 minutes for all 4 scenarios (~10x speedup)
- **Standalone Parallel Script**: ~2-4 minutes for all 4 scenarios (~15-20x speedup)

### Technical Details

- Uses fork-based parallelism via `mclapply()` (efficient on macOS/Linux)
- Pre-computes matrices outside parallel loops for efficiency
- Automatic load balancing handled by OS
- Memory-efficient: no explicit data export required with forking
- Windows compatibility: Would require `parLapply()` instead of `mclapply()`

---

## Previous Changes

See git commit history for changes prior to 2025-11-13.
