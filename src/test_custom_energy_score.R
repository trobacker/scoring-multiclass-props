# Script to test custom energy score vs. es_sample in scoringRules
# Result: small floating point arithmetic differences
# Step 1: Import dependencies
library(brms)
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
library(scoringRules)

# Step 2: Manual energy score function
manual_energy_score <- function(y, dat) {
  # y: vector of length d
  # dat: matrix with rows = dimensions, columns = samples (d × n)
  
  d <- nrow(dat)
  n <- ncol(dat)
  
  # Term 1: mean distance to observation
  diff_y <- dat - matrix(y, nrow = d, ncol = n)
  term1 <- mean(sqrt(colSums(diff_y^2)))
  
  # Term 2: mean pairwise distance among samples
  dmat <- as.matrix(dist(t(dat)))
  term2 <- mean(dmat)
  
  ES <- term1 - 0.5 * term2
  
  list(
    term1 = as.numeric(term1),
    term2 = as.numeric(term2),
    ES    = as.numeric(ES)
  )
}

# Step 3: Single-run Dirichlet simulation 
set.seed(123)

alpha_true      <- c(5, 2, 3)
alpha_correct   <- c(5, 2, 3)
alpha_incorrect <- c(2, 8, 1)

# One observed outcome
y <- as.vector(brms::rdirichlet(n = 1, alpha = alpha_true))

N <- 100
samples_correct   <- t(brms::rdirichlet(n = N, alpha_correct))
samples_incorrect <- t(brms::rdirichlet(n = N, alpha_incorrect))

es_correct   <- manual_energy_score(y, samples_correct)
es_incorrect <- manual_energy_score(y, samples_incorrect)

es_correct
es_incorrect

es_sample_correct <- scoringRules::es_sample(y = y, dat = samples_correct)
es_sample_incorrect <- scoringRules::es_sample(y = y, dat = samples_incorrect)

# Step 4: Are the custom and es_sample energy scores the same?
options(digits = 20) # increase digits to display
if( es_correct$ES == es_sample_correct && es_incorrect$ES == es_sample_incorrect){
  print("Scores match!")
} else{
  cat("ES custom correct: ", es_correct$ES, "\n")
  cat("es_sample corect: ", es_sample_correct, "\n")
  cat("Scores do not match! Likely a floating-point arithimetic error. 
  On my device, it's a difference in the 16th decimal position.
  Notice that scoringRules relies on C++ underneath, whereas this custom function is all R.")
}
options(digits = 7) # Set back to default value before script end