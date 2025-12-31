"""
Double Machine Learning (DML) demo using the `doubleml` Python package.

This is intentionally written as a simple, readable script.

What this script does
1) Loads the 4-patient toy CSV in `data/`.
2) Fits a partially linear regression (PLR) DML model with DoubleML:
   - Y = Tumor
   - D = BugA (or a log-ratio for compositional sanity)
   - X = measured confounders (Diet, and optionally other bugs as log-ratios)
3) Also runs the same setup on a larger synthetic dataset so you can see a more stable estimate.
"""

import os

import numpy as np
import pandas as pd

try:
    from doubleml import DoubleMLData, DoubleMLPLR
except Exception as exc:  # pragma: no cover
    raise SystemExit(
        "This example requires the `doubleml` package.\n"
        "Install from the repo root with: `pip install -e .`\n"
        f"Import error: {exc}"
    )

from sklearn.linear_model import LinearRegression


PSEUDOCOUNT = 1e-6


def log_ratio(df: pd.DataFrame, numerator: str, denominator: str) -> np.ndarray:
    num = df[numerator].to_numpy(dtype=float) + PSEUDOCOUNT
    den = df[denominator].to_numpy(dtype=float) + PSEUDOCOUNT
    return np.log(num / den)


def load_toy_csv() -> pd.DataFrame:
    script_dir = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(script_dir, "..", "data", "toy_microbiome_metabolite_outcome.csv")
    return pd.read_csv(path)


def simulate_dataset(n: int = 400, seed: int = 1) -> pd.DataFrame:
    rng = np.random.default_rng(seed)

    diet = rng.integers(0, 2, size=n)

    alpha_base = np.array([2.0, 4.0, 3.0, 2.0], dtype=float)
    alpha = np.where(diet[:, None] == 1, alpha_base * np.array([3.0, 0.7, 1.0, 1.0]), alpha_base)
    bugs = rng.gamma(shape=alpha, scale=1.0)
    bugs = bugs / bugs.sum(axis=1, keepdims=True)

    bug_a, bug_b, bug_c, bug_d = bugs.T
    log_a_over_d = np.log((bug_a + PSEUDOCOUNT) / (bug_d + PSEUDOCOUNT))

    metabolite = 1.5 + 0.8 * log_a_over_d + 0.3 * diet + rng.normal(0, 0.15, size=n)
    tumor = 80.0 - 12.0 * metabolite - 3.0 * log_a_over_d - 6.0 * diet + rng.normal(0, 3.0, size=n)

    return pd.DataFrame(
        {
            "Patient": [f"S{i+1}" for i in range(n)],
            "Diet": diet.astype(int),
            "BugA": bug_a,
            "BugB": bug_b,
            "BugC": bug_c,
            "BugD": bug_d,
            "Metabolite": metabolite,
            "Tumor": tumor,
        }
    )


def fit_doubleml_plr(df: pd.DataFrame, d_col: str, x_cols: list[str], n_folds: int) -> None:
    data = DoubleMLData(df, y_col="Tumor", d_cols=d_col, x_cols=x_cols)

    ml_l = LinearRegression()
    ml_m = LinearRegression()

    dml_plr = DoubleMLPLR(data, ml_l=ml_l, ml_m=ml_m, n_folds=n_folds)
    dml_plr.fit()

    print(dml_plr.summary)


def main() -> None:
    print("\nToy CSV (4 patients):")
    df_small = load_toy_csv()
    print(df_small.to_string(index=False))

    # In compositional data, a raw abundance change is easiest to misread.
    # A simple ALR-style log-ratio makes the treatment "BugA relative to BugD".
    df_small = df_small.copy()
    df_small["log_BugA_over_BugD"] = log_ratio(df_small, "BugA", "BugD")

    print("\nDoubleML PLR on the 4-patient toy dataset (purely illustrative; n is tiny):")
    fit_doubleml_plr(df_small, d_col="log_BugA_over_BugD", x_cols=["Diet"], n_folds=2)

    print("\nSynthetic dataset (n=400) with the same story:")
    df = simulate_dataset(n=400, seed=1)
    df["log_BugA_over_BugD"] = log_ratio(df, "BugA", "BugD")
    df["log_BugB_over_BugD"] = log_ratio(df, "BugB", "BugD")
    df["log_BugC_over_BugD"] = log_ratio(df, "BugC", "BugD")

    print("\nDoubleML PLR with only Diet as the confounder:")
    fit_doubleml_plr(df, d_col="log_BugA_over_BugD", x_cols=["Diet"], n_folds=5)

    print("\nDoubleML PLR with Diet + other bugs (as log-ratios) as confounders:")
    fit_doubleml_plr(df, d_col="log_BugA_over_BugD", x_cols=["Diet", "log_BugB_over_BugD", "log_BugC_over_BugD"], n_folds=5)


if __name__ == "__main__":
    main()
