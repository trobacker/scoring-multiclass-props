# Pre-computing term2 for Energy Score Optimization

## Overview

This document explores a computational optimization for energy score calculations in simulation experiments that can reduce computational complexity by orders of magnitude while maintaining statistical validity.

## Key Insight: term2 is Independent of Observations

The energy score decomposes into two terms:
- **term1**: Mean distance from forecast samples to observation (depends on both forecast and observation)
- **term2**: 0.5 × mean pairwise distance among forecast samples (depends ONLY on forecast samples)

```r
# Term 1: depends on both y (observation) and dat (forecast samples)
diff_y <- dat - matrix(y, nrow = d, ncol = n)
term1 <- mean(sqrt(colSums(diff_y^2)))

# Term 2: depends ONLY on dat (forecast samples)
dmat <- as.matrix(dist(t(dat)))
term2 <- 0.5 * mean(dmat)

# Energy Score
ES <- term1 - term2
```

Since term2 depends only on the forecast distribution (not the observation), it represents an intrinsic property of the forecast's spread/sharpness.

## Proposed Optimization

### Current Approach
- **es_mn_sampling()**: n_samp=100 thetas, N_multinomial=100 per theta → 10,000 samples
- **Recalculate term2** for each of 1,000 simulation iterations
- **Complexity**: O(n²) × iterations = 10,000² × 1,000 = 10¹¹ operations

### Optimized Approach
- **Generate once**: n_samp=10,000 thetas, N_multinomial=1 per theta → 10,000 samples
- **Pre-compute term2 ONCE** before simulations
- **Each iteration**: Only compute term1 (which depends on observation)
- **Complexity**: O(n²) × 1 = 10,000² × 1 = 10⁸ operations

**Result: 1000× reduction in term2 computation**

## Computational Benefits

**term2 calculation complexity:** O(n²) for n samples

| Component | Current | Optimized | Reduction |
|-----------|---------|-----------|-----------|
| term2 calculations per scenario | 1,000 iterations | 1 time | 1000× |
| Total pairwise distances | 10¹¹ | 10⁸ | 1000× |
| Can scale to iterations | ~1,000 | ~10,000+ | 10× |

## Conceptual Validity

This optimization is valid from a proper scoring rule perspective:

### 1. term2 as Distributional Property
With a large sample (n=10,000), term2 converges to the **expected pairwise distance** in the forecast distribution—a fixed property of:
- Dirichlet(α) for non-multinomial sampling
- Dirichlet-Multinomial(α, n_counts) for multinomial sampling

### 2. term1 as Observation-Specific
Must be recalculated for each new observation since it measures forecast-observation alignment.

### 3. Maintaining Forecast Variation
We can still vary forecasts for term1 calculation while using the pre-computed expected term2:
- term2 represents the distributional property (expected spread)
- term1 samples can still vary across simulations
- Separation of concerns: distributional properties vs. observation-specific scoring

## Implementation Options

### Option 1: Pre-compute term2 as Scalar (Recommended)

Most efficient approach - pre-compute term2 once per (α, n_counts, method) combination:

```r
precompute_term2 <- function(alpha, n_counts, n_samp = 10000, method = "mn") {
  K <- length(alpha)

  if (method == "mn") {
    # For multinomial sampling method
    thetas <- lapply(FUN = draw_one_dirichlet, X = rep(1, n_samp), alpha = alpha)
    # Draw 1 multinomial sample per theta
    forecast_samples <- do.call(cbind,
      lapply(thetas, function(theta) rmultinom(n = 1, size = n_counts, prob = theta))
    )
    p_matrix <- forecast_samples / n_counts
  } else {
    # For non-multinomial method
    thetas <- lapply(FUN = draw_one_dirichlet, X = rep(1, n_samp), alpha = alpha)
    p_matrix <- matrix(unlist(thetas), nrow = K, byrow = FALSE)
  }

  # Compute term2 once
  dmat <- as.matrix(dist(t(p_matrix)))
  term2 <- 0.5 * mean(dmat)

  return(term2)
}
```

**Usage in simulations:**

```r
# Pre-computation phase (once per scenario)
term2_correct_nonmn <- precompute_term2(alpha, n_counts, n_samp=10000, method="nonmn")
term2_incorrect_nonmn <- precompute_term2(alpha_wrong, n_counts, n_samp=10000, method="nonmn")
term2_correct_mn <- precompute_term2(alpha, n_counts, n_samp=10000, method="mn")
term2_incorrect_mn <- precompute_term2(alpha_wrong, n_counts, n_samp=10000, method="mn")

# Simulation loop (can now scale to 10,000+ iterations)
for (i in 1:nsim) {
  # Generate observation
  prob_vec <- as.numeric(rdirichlet(1, alpha))
  C <- rmultinom(n = 1, size = n_counts, prob = prob_vec)
  p <- C / n_counts

  # Generate forecast sample for term1 calculation only
  thetas <- lapply(FUN = draw_one_dirichlet, X = rep(1, n_samp_term1), alpha = alpha)
  p_matrix <- matrix(unlist(thetas), nrow = K, byrow = FALSE)

  # Calculate only term1 (term2 is pre-computed)
  diff_y <- p_matrix - matrix(p, nrow = K, ncol = n_samp_term1)
  term1 <- mean(sqrt(colSums(diff_y^2)))

  # Use pre-computed term2
  ES <- term1 - term2_correct_nonmn
}
```

### Option 2: Pre-compute and Store Full Sample (Memory Intensive)

Store the entire forecast sample matrix for potential reuse:

```r
precompute_forecast_samples <- function(alpha, n_counts, n_samp = 10000, method = "mn") {
  # ... generate p_matrix as above ...

  dmat <- as.matrix(dist(t(p_matrix)))
  term2 <- 0.5 * mean(dmat)

  list(
    samples = p_matrix,
    term2 = term2
  )
}
```

**Pros**: Allows reusing samples for term1 calculation
**Cons**: Memory intensive for large n_samp

## Forecast Sampling Strategy Considerations

### Current Strategy (n_samp=100, N_multinomial=100)
- **Interpretation**: "100 forecasters, each submitting 100 multinomial draws"
- Better captures multinomial sampling variability per theta
- Each theta value gets more representation in the sample

### Proposed Strategy (n_samp=10,000, N_multinomial=1)
- **Interpretation**: "10,000 forecasters, each submitting 1 multinomial draw"
- Better explores the theta space (more diverse thetas)
- Relies on law of large numbers to capture multinomial variability

**Both are valid**, but they emphasize different aspects:
- **Current**: Better for small number of distinct forecast distributions
- **Proposed**: Better for capturing distributional properties with large samples

## Statistical Implications

### Benefits of n_samp=10,000 for term2:
- **Very stable estimate** of term2 (low variance)
- More representative of the **true** forecast distribution spread
- Reduces noise in term2, making term1 differences more interpretable
- Separates distributional properties from observation-specific variability

### Trade-offs:
- Pre-computed term2 uses expected spread across all simulations
- Loses some sample-to-sample variability in term2
- **Resolution**: This is actually desirable! term2 represents a property of the distribution, not individual samples

## Applies to Both Functions

**Important**: This optimization speeds up **BOTH** sampling methods:

1. **es_sampling()** (non-multinomial):
   - Pre-compute term2 from pure Dirichlet samples
   - Each iteration only computes term1

2. **es_mn_sampling()** (multinomial):
   - Pre-compute term2 from multinomial-sampled thetas
   - Each iteration only computes term1

Both functions calculate term2 the same way (pairwise distances), just on different types of samples. The computational savings apply equally to both.

## Implementation Recommendations

### For Production Simulations:

1. **Pre-computation phase**: Calculate term2 for each (α, n_counts, method) combination
   - α_correct, non-mn: term2_correct_nonmn[n_counts]
   - α_wrong, non-mn: term2_incorrect_nonmn[n_counts]
   - α_correct, mn: term2_correct_mn[n_counts]
   - α_wrong, mn: term2_incorrect_mn[n_counts]

2. **Simulation loop**: Use pre-computed term2, only calculate term1 per iteration

3. **Sample size for term1**: Can use smaller samples (e.g., 1,000) for further speedup since we only need it for term1

4. **Scale up simulations**: With these savings, can feasibly run 10,000+ simulation iterations

### Key Design Decisions:

1. **n_samp for term2 pre-computation**: 10,000 recommended for stability
2. **n_samp for term1 per iteration**: 1,000-5,000 recommended (faster, still accurate)
3. **N_multinomial**: 1 per theta for mn method (explores theta space better)
4. **Verification**: Create test comparing dynamic vs pre-computed term2 to validate

## Important Considerations

### 1. n_counts Dependency
term2 differs for each n_counts value (especially for multinomial method), so must pre-compute separately for each n_counts in the experimental design.

### 2. Compatibility with Current Work
This optimization is **compatible** with recent work on varying forecasts:
- We still vary forecasts across simulations (for term1)
- term2 is a distributional property (should be constant)
- Separates distributional properties from observation-specific scoring

### 3. Verification Strategy
Before deploying to production simulations:
1. Create test script comparing pre-computed vs dynamic term2
2. Verify that large-sample term2 is stable
3. Confirm ES results are equivalent
4. Validate speedup claims

## References

- Gneiting, T., & Raftery, A. E. (2007). Strictly proper scoring rules, prediction, and estimation. *Journal of the American Statistical Association*, 102(477), 359-378.

## Date
2024-03-24
