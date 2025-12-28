# Mediation analysis using the `mediation` package (ACME / ADE / total effect)
#
# Run from the repo root:
#   Rscript r/02_mediation_mediation_pkg.R
#
# Install if needed (outside this repo): install.packages("mediation")

if (!requireNamespace("mediation", quietly = TRUE)) {
  stop("Package 'mediation' is not installed. Install it with install.packages('mediation').")
}

library(mediation)

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else ""
  if (nzchar(script_path)) dirname(normalizePath(script_path)) else getwd()
}

load_toy_csv <- function() {
  script_dir <- get_script_dir()
  csv_path <- file.path(script_dir, "..", "data", "toy_microbiome_metabolite_outcome.csv")
  read.csv(csv_path)
}

run_mediation <- function(df, treat_value, control_value, sims = 1000) {
  # Two regression models:
  # 1) mediator model: M ~ X + Z
  # 2) outcome model:  Y ~ X + M + Z
  m_fit <- lm(Metabolite ~ BugA + Diet, data = df)
  y_fit <- lm(Tumor ~ BugA + Metabolite + Diet, data = df)

  med <- mediate(
    model.m = m_fit,
    model.y = y_fit,
    treat = "BugA",
    mediator = "Metabolite",
    treat.value = treat_value,
    control.value = control_value,
    sims = sims
  )

  return(med)
}

df_small <- load_toy_csv()
print(df_small)

# With continuous treatments, `mediate()` compares two chosen values of BugA.
# Here we use "low" vs "high" from the observed toy data.
control_value <- as.numeric(quantile(df_small$BugA, probs = 0.25))
treat_value <- as.numeric(quantile(df_small$BugA, probs = 0.75))

set.seed(1)
med_small <- run_mediation(df_small, treat_value = treat_value, control_value = control_value, sims = 1000)

cat("\nMediation on the 4-patient toy dataset (purely illustrative; n is tiny):\n")
print(summary(med_small))

# Optional: a larger synthetic dataset with the same story (for a more stable demo).
rdirichlet_simple <- function(n, alpha) {
  draws <- matrix(rgamma(n * length(alpha), shape = alpha, rate = 1), nrow = n)
  draws / rowSums(draws)
}

simulate_dataset <- function(n = 400, seed = 1) {
  set.seed(seed)
  Diet <- rbinom(n, size = 1, prob = 0.5)

  alpha_base <- c(2, 4, 3, 2)
  alpha_hi <- alpha_base * c(3.0, 0.7, 1.0, 1.0)
  alpha <- t(sapply(Diet, function(z) if (z == 1) alpha_hi else alpha_base))
  bugs <- t(sapply(1:n, function(i) rdirichlet_simple(1, alpha[i, ])))

  BugA <- bugs[, 1]
  BugB <- bugs[, 2]
  BugC <- bugs[, 3]
  BugD <- bugs[, 4]

  logAoverD <- log((BugA + 1e-6) / (BugD + 1e-6))

  Metabolite <- 1.5 + 0.8 * logAoverD + 0.3 * Diet + rnorm(n, 0, 0.15)
  Tumor <- 80 - 12 * Metabolite - 3 * logAoverD - 6 * Diet + rnorm(n, 0, 3)

  data.frame(
    Patient = paste0("S", seq_len(n)),
    Diet = Diet,
    BugA = BugA,
    BugB = BugB,
    BugC = BugC,
    BugD = BugD,
    Metabolite = Metabolite,
    Tumor = Tumor
  )
}

df <- simulate_dataset(n = 400, seed = 1)
control_value <- as.numeric(quantile(df$BugA, probs = 0.25))
treat_value <- as.numeric(quantile(df$BugA, probs = 0.75))

set.seed(1)
med <- run_mediation(df, treat_value = treat_value, control_value = control_value, sims = 500)
cat("\nMediation on a synthetic dataset (n=400) with the same story:\n")
print(summary(med))
