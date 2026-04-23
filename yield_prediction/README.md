# Machine Learning for Reaction Yield Prediction

This project applies linear regression models to predict reaction yields from quantum-chemically computed energies of reaction intermediates.

## Overview

- **Input (X.csv)**: Absolute energies of 52 reaction intermediates for each reaction condition (backbone–pyridone combination)
- **Output (y.csv)**: Reaction yield (`beta_av`, %) for each backbone–pyridone combination
- **Dataset size**: 104 samples (21 backbone variants × up to 13 pyridone variants)

## Data Preprocessing (`dataset.ipynb`)

1. **Relative energies**: For each sample, the absolute energy of intermediate `0` (reference state) is subtracted from all other intermediates to obtain relative energies (`*_zero` features).
2. **Logit transformation**: The yield `beta_av` is transformed as `log(beta_av / (100 - beta_av))` to map the bounded percentage onto the real line (`beta_av_logit`), which is used as the regression target.
3. The combined dataset is saved to `dataset.parquet`.

## Models

Four linear regression models are evaluated, all using the same nested cross-validation and recursive feature elimination (RFE) framework:

| Directory | Model | Hyperparameter |
|-----------|-------|---------------|
| `OLS/` | Ordinary Least Squares | — |
| `Ridge/` | Ridge Regression | `alpha` |
| `Lasso/` | Lasso Regression | `alpha` |
| `PLS/` | Partial Least Squares Regression | `n_components` |

## Methodology

### Feature Engineering

Relative energies (intermediate energy − reference state energy) are used as features. Features are standardized with `StandardScaler` applied to both X and y before fitting.

### Recursive Feature Elimination (RFE)

Starting from all 52 relative energy features, features are eliminated one by one (the feature with the smallest absolute coefficient is removed at each step) until 4 features remain. Models are evaluated at each step to track performance vs. number of features.

### Nested Cross-Validation

- **Outer loop**: Leave-one-group-out cross-validation over pyridone groups (each pyridone held out as the test set in turn).
- **Inner loop**: Leave-one-group-out cross-validation over the remaining pyridone groups for hyperparameter optimization (Ridge/Lasso `alpha`, PLS `n_components`).
- Predictions and coefficients are recorded at each RFE step.

### Evaluation

Model performance is measured by **RMSE** on held-out test sets. Results are saved as:
- `*_train_val_preds.parquet`: predictions on training/validation data
- `*_test_preds.parquet`: predictions on test data
- `*_coefs.parquet`: model coefficients

Scatter plots and RMSE vs. number of features plots are generated for each model.

## Ridge Regression — Additional Experiments (`Ridge/fit/`)

| Notebook | Description |
|----------|-------------|
| `001_fit` | Standard nested CV, outer fold over pyridone groups |
| `004_fit_MPAA` | Outer fold over backbone groups (MPAA-style split) |
| `007_fit_random` | Negative control: target labels shuffled randomly (seed 42) |
| `010_fit_int0` | Single-feature model using only the reference intermediate energy |

## Ridge Regression — Refitting on Full Dataset (`Ridge/refit/`)

After cross-validation, the Ridge model is refit on the entire dataset for interpretability analysis:

| Notebook | Description |
|----------|-------------|
| `001_refit` | Refit on all data with CV-based hyperparameter selection |
| `002_RMSE` | RMSE of the refitted model |
| `003_coefs` | Visualization of model coefficients (bar and heatmap) |
| `004_feature_correlation` | Pairwise correlation among selected features |
| `005_feature_coefficient_correlation` | Correlation between feature values and fitted coefficients |
| `006_refit_noise` | Robustness test: Gaussian noise added to labels (σ = 0.1 – 1.0) |
| `007_RMSE_noise` | RMSE summary across noise levels |
| `008_coefs_noise` | Coefficient stability across noise levels |

## Model Comparison (`all_models/fit/`)

| Notebook | Description |
|----------|-------------|
| `001_RMSE` | Comparison of test RMSE across OLS, Ridge, Lasso, and PLS |
| `002_coefs` | Comparison of model coefficients across all models |

## Project Structure

```
machine-learning/
├── X.csv                    # Intermediate energies (features)
├── y.csv                    # Reaction yields (target)
├── dataset.ipynb            # Data preprocessing → dataset.parquet
├── dataset.parquet          # Preprocessed dataset
├── OLS/fit/                 # Ordinary least squares
├── Ridge/
│   ├── fit/                 # Ridge CV experiments
│   └── refit/               # Ridge refit & interpretability
├── Lasso/fit/               # Lasso regression
├── PLS/fit/                 # Partial least squares
├── all_models/fit/          # Cross-model comparison
├── requirements.txt
├── Dockerfile
└── docker-compose.dev.yml
```

## Environment

Python 3.10 with the following key packages:

- `scikit-learn` 1.5.2
- `polars` 1.19.0
- `pandas` 1.5.3
- `numpy` 1.24.2
- `jupyterlab` 4.0.10
- `matplotlib` 3.10.7
- `seaborn` 0.13.2

### Running with Docker

```bash
docker compose -f docker-compose.dev.yml up
```

JupyterLab will be available at `http://localhost:8888`.
