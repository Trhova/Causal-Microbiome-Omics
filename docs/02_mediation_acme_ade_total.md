# Mediation Analysis: ACME, ADE, and Total Effects

This document dives deeper into causal mediation analysis, using the R `mediation` package to decompose effects into direct and indirect components.

## Mediation in Our Toy Example

We hypothesize: BugA → Metabolite → Tumor

- **Total effect**: Overall impact of BugA on Tumor.
- **Indirect effect (ACME)**: Portion through Metabolite.
- **Direct effect (ADE)**: Portion not through Metabolite.

## The Two Models

1. **Mediator Model**: Metabolite ~ BugA + Diet
   - Estimates how BugA affects Metabolite, adjusting for confounder Diet.

2. **Outcome Model**: Tumor ~ BugA + Metabolite + Diet
   - Estimates how Metabolite and BugA affect Tumor, adjusting for Diet.

## Running Mediation in R

See `r/02_mediation_mediation_pkg.R` for code.

Key output:
- ACME: Indirect effect via Metabolite.
- ADE: Direct effect of BugA on Tumor.
- Total: ACME + ADE.

## Interpretation

Assuming assumptions hold (no unmeasured confounding affecting Metabolite-Tumor link), ACME tells us the effect transmitted through the metabolite pathway.

## Assumptions

Mediation needs stronger assumptions than “just estimating a total effect”. In plain language, you need to believe:

- After adjusting for measured confounders (here: `Diet`), there are no hidden common causes of `BugA` and `Metabolite` (X→M confounding).
- After adjusting for measured confounders **and** `BugA`, there are no hidden common causes of `Metabolite` and `Tumor` (M→Y confounding).
- You are not accidentally conditioning on a collider (e.g., selection into “has metabolomics measured”) when doing the analysis.

If these assumptions fail, ACME/ADE can be a systematically wrong estimate (biased).
