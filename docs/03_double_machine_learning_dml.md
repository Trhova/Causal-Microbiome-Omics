# Double Machine Learning (DML) for Causal Effects

This document explains Double Machine Learning, a method to estimate causal effects in high-dimensional settings using machine learning for nuisance parameter estimation.

Reference for the core idea:
- Chernozhukov et al. (2018), "Double/debiased machine learning for treatment and structural parameters": https://arxiv.org/abs/1608.00060

## DML in Our Toy Example

Estimate the effect of BugA on Tumor, adjusting for Diet, BugB, BugC, BugD.

## Steps

1. **Partial out confounders from X**: Predict BugA from [Diet, BugB, BugC, BugD] using ML. Compute residuals: BugA_residual = BugA - predicted_BugA.

2. **Partial out confounders from Y**: Predict Tumor from [Diet, BugB, BugC, BugD] using ML. Compute residuals: Tumor_residual = Tumor - predicted_Tumor.

3. **Final regression**: Regress Tumor_residual ~ BugA_residual. The coefficient is the causal effect estimate.

## Cross-Fitting

To avoid overfitting bias:
- Split data into folds (e.g., 2 folds for toy data).
- Train ML models on one fold, compute residuals on the other.
- Average estimates across folds.

## Why It Works

ML handles complex relationships with confounders, but the final simple regression ensures the effect estimate is unbiased under assumptions.

## Code

See `python/02_dml_example.py` for a runnable demo using the `doubleml` package:
- Fits a partially linear regression DML model via `DoubleMLPLR`
- Uses simple learners (`sklearn` linear regression) to keep the example easy to read
- Also generates a larger synthetic dataset with the same variable names, since n=4 is too small to be stable
