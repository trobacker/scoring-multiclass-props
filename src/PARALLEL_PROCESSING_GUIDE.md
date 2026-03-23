# Parallel Processing Guide for Simulation Scripts

## The Problem: Forking Issues in GUI Environments

When running R simulations with parallel processing on **macOS** in **GUI environments** (RStudio, Positron), you may encounter errors or warnings related to forking. This is because:

1. **Fork safety**: macOS has restrictions on forking processes from GUI applications
2. **Memory issues**: Forked processes in GUIs can cause memory corruption or crashes
3. **Silent failures**: Sometimes parallel operations fail silently, returning NULL results

### Symptoms
- No data saved despite script completion
- Warnings about forking or process failures
- Error messages mentioning `mclapply` or `fork()`
- Results contain NULL or incomplete data

---

## The Solution

Both simulation scripts (`run_simulations_alpha_asymmetric.R` and `run_simulations_alpha_symmetric.R`) now **automatically detect** your environment and adjust accordingly.

### Automatic Detection

The scripts detect:
- ✅ **RStudio**: `RSTUDIO=1` environment variable
- ✅ **Positron**: `POSITRON=1` environment variable
- ✅ **Terminal**: Neither variable set

### Behavior by Environment

| Environment | Scenario Processing | Cores Used | Safety Level |
|-------------|-------------------|------------|--------------|
| **Positron/RStudio** | Sequential | 4 cores max | High (conservative) |
| **Terminal (Rscript)** | Parallel | All available - 2 | Normal |

#### In Positron/RStudio
```
DETECTED POSITRON ENVIRONMENT
=========================================
Running scenarios SEQUENTIALLY to avoid forking issues on macOS
(Using parallelization within each scenario for speedup)

Using 4 cores per scenario (conservative setting for GUI)
```

- Scenarios run **one at a time** (no outer parallelization)
- Each scenario uses **up to 4 cores** internally (safe for forking)
- More reliable, slightly slower (~20-30 mins for 5 scenarios)

#### In Terminal
```
Running 5 scenarios in parallel...
```

- Scenarios run **in parallel** (up to 4 scenarios simultaneously)
- Each scenario uses **many cores** (total cores - 2)
- Fastest possible (~7-12 mins for 5 scenarios)

---

## Best Practices

### For Reliability (Recommended for Positron/RStudio)

✅ **Just run the script normally** - it will automatically use safe settings:

```r
# In Positron/RStudio R Console
source("src/run_simulations_alpha_symmetric.R")
```

The script detects Positron and:
- Runs scenarios sequentially
- Uses only 4 cores per scenario
- Avoids fork safety issues

### For Maximum Speed (Terminal)

If you need faster execution and can use the terminal:

```bash
# From project root
Rscript src/run_simulations_alpha_symmetric.R

# Or from src/ directory
cd src
Rscript run_simulations_alpha_symmetric.R
```

This uses full parallelization and all available cores.

---

## Manual Core Control (Advanced)

If you still encounter issues or want more control, you can manually edit the script:

### Option 1: Force Sequential Mode

Find this section in the script:
```r
in_gui <- in_rstudio || in_positron
```

Change to:
```r
in_gui <- TRUE  # Force safe mode
```

### Option 2: Reduce Cores Further

Find this line:
```r
n_cores_safe <- max(1, min(n_cores_inner, 4))
```

Change to use fewer cores:
```r
n_cores_safe <- 2  # Use only 2 cores
```

### Option 3: Disable Inner Parallelization

In the `run_scenario()` function call, set:
```r
run_scenario(scenario_config,
             n_samp = 100,
             nsim = 1000,
             n_cores_inner = 1)  # No parallelization
```

---

## Understanding the Parallelization Levels

The scripts use **two levels** of parallelization:

### Level 1: Outer (Scenario-level)
- Runs multiple scenarios simultaneously
- Uses `parallel::mclapply` with `mc.cores = n_cores_outer`
- **Disabled in GUI environments** (forking issues)

### Level 2: Inner (Iteration-level)
- Within each scenario, parallelizes the 1000 simulation iterations
- Uses `parallel::mclapply` with `mc.cores = n_cores_inner`
- **Limited to 4 cores in GUI environments** (safer)
- Used for both `es_sampling()` and `es_mn_sampling()` functions

### Core Allocation Examples

**Your system (14 cores):**

| Mode | Outer Cores | Inner Cores | Total Usage | Time Estimate |
|------|-------------|-------------|-------------|---------------|
| Positron (safe) | 1 (sequential) | 4 | 4 cores | ~25 mins |
| Terminal (fast) | 4 scenarios | 12 per scenario | Up to 14 | ~10 mins |

---

## Troubleshooting

### Problem: Script hangs or never completes
**Solution**: The script may be waiting for parallel processes. Restart and ensure you're in safe mode (Positron should auto-detect).

### Problem: "cannot fork" error
**Solution**: This confirms forking issues. The updated script should prevent this, but if it persists:
1. Update both scripts (they now have Positron detection)
2. Try restarting R/Positron
3. Consider running from terminal

### Problem: Results are empty or NULL
**Solution**: Parallel processes likely failed silently. Use the updated scripts which:
- Detect your environment automatically
- Use conservative settings in GUIs
- Run scenarios sequentially in Positron

### Problem: Still getting warnings after updates
**Solution**:
1. Check you're running the latest version of the scripts
2. Verify Positron is detected: Look for "DETECTED POSITRON ENVIRONMENT" message
3. If detection fails, manually force safe mode (see Option 1 above)

---

## Performance Comparison

**5 scenarios × 10 sample sizes × 1000 iterations × 2 methods = 100,000 simulations**

| Environment | Strategy | Est. Time | Reliability |
|-------------|----------|-----------|-------------|
| **Positron (auto)** | Sequential + 4 cores | 20-30 min | ⭐⭐⭐⭐⭐ High |
| **Terminal (auto)** | Parallel + all cores | 7-12 min | ⭐⭐⭐⭐⭐ High |
| **Old script (Positron)** | Parallel (buggy) | Failed ❌ | ⭐ Low |

---

## Technical Details: Why Forking Fails

### The Fork Problem
`parallel::mclapply` uses **forking** on Unix/macOS:
1. Creates child process by copying parent's memory
2. Child runs computation
3. Parent collects results

### Why GUIs Break This
- macOS restricts forking from GUI applications
- GUIs have complex state that doesn't fork cleanly
- UI frameworks (Cocoa) are not fork-safe
- Results in crashes, hangs, or silent failures

### Why Terminal Works
- No GUI framework to interfere
- Clean process environment
- macOS restrictions are less strict
- Forking works as intended

### Our Solution
- Detect GUI environment automatically
- Use sequential processing (avoid outer forking)
- Limit inner forking to safe level (4 cores)
- Maintain parallelization where safe

---

## When to Use Each Mode

### Use GUI Mode (Positron/RStudio) When:
- ✅ Developing/debugging code
- ✅ Interactive exploration
- ✅ Single scenario testing
- ✅ You need to see intermediate output
- ✅ Running on laptop with limited battery

### Use Terminal Mode When:
- ✅ Final production runs
- ✅ Running multiple experiments
- ✅ Maximum performance needed
- ✅ Running on server/cluster
- ✅ Batch processing scenarios

---

## Summary

**The fix is automatic!** The updated scripts:
1. ✅ Detect Positron and RStudio automatically
2. ✅ Use safe settings in GUI environments
3. ✅ Use full performance in terminal
4. ✅ Limit cores to prevent crashes
5. ✅ Provide clear status messages

You should now be able to run simulations successfully from Positron without any manual configuration. If you need maximum speed, run from terminal with `Rscript`.
