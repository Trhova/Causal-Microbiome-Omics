# Load and inspect the toy dataset

# Robust path: works when run as `Rscript r/01_toy_data.R` from repo root.
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else ""
script_dir <- if (nzchar(script_path)) dirname(normalizePath(script_path)) else getwd()
csv_path <- file.path(script_dir, "..", "data", "toy_microbiome_metabolite_outcome.csv")

df <- read.csv(csv_path)

# View the data
print(df)

# Summary statistics
summary(df)
