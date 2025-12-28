"""
Double Machine Learning (DML) intuition demo (no extra dependencies).

What this script does
1) Loads the 4-patient toy CSV in `data/`.
2) Demonstrates *cross-fitted residualization*:
   - predict X from Z on train fold, compute X residuals on held-out fold
   - predict Y from Z on train fold, compute Y residuals on held-out fold
   - regress Y_residual ~ X_residual

Notes
- With only 4 patients, estimates are noisy. For a more stable demo, this script also
  generates a larger synthetic dataset that follows the same story and variable names.
- This is educational code. For real analyses, consider libraries like `doubleml`/`econml`,
  and think hard about the adjustment set using a DAG.
"""

from __future__ import annotations

import os
from dataclasses import dataclass

import numpy as np
import pandas as pd


def _script_dir() -> str:
    return os.path.dirname(os.path.abspath(__file__))


def load_toy_csv() -> pd.DataFrame:
    path = os.path.join(_script_dir(), "..", "data", "toy_microbiome_metabolite_outcome.csv")
    df = pd.read_csv(path)
    return df


def alr_log_ratio(df: pd.DataFrame, numerator: str, denominator: str, pseudocount: float = 1e-6) -> np.ndarray:
    """
    ALR-style log-ratio: log(numerator/denominator).

    Why: microbiome relative abundances are compositional (they sum to 1), so a "BugA increase"
    is only meaningful *relative* to something else.
    """
    num = df[numerator].to_numpy(dtype=float) + pseudocount
    den = df[denominator].to_numpy(dtype=float) + pseudocount
    return np.log(num / den)


def add_intercept(X: np.ndarray) -> np.ndarray:
    return np.column_stack([np.ones(X.shape[0]), X])


@dataclass(frozen=True)
class OLSFit:
    coef: np.ndarray

    def predict(self, X: np.ndarray) -> np.ndarray:
        return X @ self.coef


def fit_ols(X: np.ndarray, y: np.ndarray) -> OLSFit:
    """
    Ordinary Least Squares via least squares.

    This is a stand-in for "ML": in real DML you could replace this with random forests,
    gradient boosting, etc. Cross-fitting is the key idea.
    """
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    return OLSFit(coef=coef)


def kfold_indices(n: int, k: int, seed: int = 0) -> list[np.ndarray]:
    rng = np.random.default_rng(seed)
    perm = rng.permutation(n)
    return [perm[i::k] for i in range(k)]


def cross_fitted_residuals(y: np.ndarray, Z: np.ndarray, folds: list[np.ndarray]) -> np.ndarray:
    """
    Cross-fitted residuals: each point's prediction is made by a model that did *not* train on that point.
    """
    n = y.shape[0]
    y_hat = np.empty(n, dtype=float)
    all_idx = np.arange(n)
    for test_idx in folds:
        train_mask = np.ones(n, dtype=bool)
        train_mask[test_idx] = False
        train_idx = all_idx[train_mask]

        model = fit_ols(add_intercept(Z[train_idx]), y[train_idx])
        y_hat[test_idx] = model.predict(add_intercept(Z[test_idx]))
    return y - y_hat


def effect_from_residuals(y_resid: np.ndarray, x_resid: np.ndarray) -> float:
    """
    Final stage of (partially linear) DML: regress y_resid on x_resid.
    """
    X = add_intercept(x_resid.reshape(-1, 1))
    fit = fit_ols(X, y_resid)
    # coef[0] is intercept, coef[1] is effect estimate
    return float(fit.coef[1])


def run_dml(df: pd.DataFrame, x: np.ndarray, y: np.ndarray, Z: np.ndarray, k_folds: int, seed: int) -> float:
    folds = kfold_indices(len(df), k=k_folds, seed=seed)
    x_resid = cross_fitted_residuals(x, Z, folds)
    y_resid = cross_fitted_residuals(y, Z, folds)
    return effect_from_residuals(y_resid, x_resid)


def simulate_dataset(n: int = 200, seed: int = 0) -> pd.DataFrame:
    """
    Synthetic data with the same variable names and the same intended story:
      Diet (Z) -> BugA (X)
      Diet (Z) -> Tumor (Y)
      BugA (X) -> Metabolite (M) -> Tumor (Y)

    We also keep microbiome compositional (A+B+C+D=1) and use a log-ratio internally.
    """
    rng = np.random.default_rng(seed)

    diet = rng.integers(0, 2, size=n)

    # Dirichlet composition; Diet shifts mass toward BugA.
    alpha_base = np.array([2.0, 4.0, 3.0, 2.0], dtype=float)
    alpha = np.where(diet[:, None] == 1, alpha_base * np.array([3.0, 0.7, 1.0, 1.0]), alpha_base)
    bugs = rng.gamma(shape=alpha, scale=1.0)
    bugs = bugs / bugs.sum(axis=1, keepdims=True)

    bug_a, bug_b, bug_c, bug_d = bugs.T

    log_a_over_d = np.log((bug_a + 1e-6) / (bug_d + 1e-6))

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


def residualize_by_group_mean(values: np.ndarray, groups: np.ndarray) -> np.ndarray:
    """
    Toy residualization used in the README: predict by group mean, then take leftovers.
    """
    df = pd.DataFrame({"v": values, "g": groups})
    mean_by_group = df.groupby("g", observed=True)["v"].transform("mean").to_numpy(dtype=float)
    return values - mean_by_group


def main() -> None:
    df_small = load_toy_csv()
    print("Toy CSV (4 patients):")
    print(df_small.to_string(index=False))
    print()

    # With n=4, we focus on the "leftovers" intuition using the Diet-group mean residuals
    # (this matches the table in the README). Cross-fitting is hard to demonstrate with n=4.
    diet = df_small["Diet"].to_numpy(dtype=int)
    x_raw = df_small["BugA"].to_numpy(dtype=float)
    y = df_small["Tumor"].to_numpy(dtype=float)

    x_resid = residualize_by_group_mean(x_raw, diet)
    y_resid = residualize_by_group_mean(y, diet)
    eff = effect_from_residuals(y_resid, x_resid)

    demo = df_small[["Patient", "Diet", "BugA", "Tumor"]].copy()
    demo["BugA_leftover"] = x_resid
    demo["Tumor_leftover"] = y_resid

    print("Leftovers after adjusting for Diet by group means (toy intuition):")
    print(demo.to_string(index=False))
    print(f"\nEffect from leftovers (Tumor_leftover ~ BugA_leftover): {eff:.3f}")
    print()

    df = simulate_dataset(n=400, seed=1)
    x_raw = df["BugA"].to_numpy(dtype=float)
    x_logratio = alr_log_ratio(df, "BugA", "BugD")
    y = df["Tumor"].to_numpy(dtype=float)
    Z = np.column_stack(
        [
            df["Diet"].to_numpy(dtype=float),
            alr_log_ratio(df, "BugB", "BugD"),
            alr_log_ratio(df, "BugC", "BugD"),
        ]
    )

    eff_raw = run_dml(df, x_raw, y, Z, k_folds=5, seed=123)
    eff_lr = run_dml(df, x_logratio, y, Z, k_folds=5, seed=123)
    print("Synthetic data (n=400) with the same story:")
    print(f"DML-style estimate (raw BugA): {eff_raw:.3f}")
    print(f"DML-style estimate (log(BugA/BugD)): {eff_lr:.3f}")


if __name__ == "__main__":
    main()
