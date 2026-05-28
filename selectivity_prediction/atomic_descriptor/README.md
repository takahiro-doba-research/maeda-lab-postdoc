# Selectivity prediction with atomic descriptor

A machine learning pipeline for predicting palladium-catalyzed cross-coupling reaction selectivity (β/α ratio) using atom-based molecular descriptors.

## Overview

This project builds predictive models on 104 reaction data points (8 MPAA × 13 pyridones) by:

1. Featurizing molecules with xTB (semi-empirical quantum chemistry) and morfeus (steric/electronic descriptors)
2. Training Ridge regression models with Recursive Feature Elimination (RFE)
3. Evaluating models via Leave-One-Group-Out cross-validation in two settings: pyridone-based and MPAA-based
4. Visualizing model coefficients to identify the most important features

## Directory Structure

```
atomic_descriptor/
├── backbone/             # Gaussian input files (*.com) for MPAA molecules
├── pyridone/             # Gaussian input files (*.com) for pyridone molecules
├── 000_dataset.ipynb     # Feature computation and dataset construction
├── 001_fit.ipynb         # Model training (Leave-One-Pyridone-Out CV)
├── 002_RMSE.ipynb        # RMSE evaluation and scatter plots (Pyridone CV)
├── 003_coefs.ipynb       # Coefficient visualization (Pyridone CV)
├── 004_fit_MPAA.ipynb    # Model training (Leave-One-MPAA-Out CV)
├── 005_RMSE_MPAA.ipynb   # RMSE evaluation and scatter plots (MPAA CV)
├── 006_coefs_MPAA.ipynb  # Coefficient visualization (MPAA CV)
├── dataset.parquet       # Combined dataset
├── features_backbone.parquet  # MPAA features
├── features_pyridone.parquet  # Pyridone features
├── y.csv                 # Target variable (MPAA, pyridone, beta_av)
├── requirements.txt
└── Dockerfile
```

## Notebook Execution Order

```
000_dataset.ipynb
    ↓ generates dataset.parquet
001_fit.ipynb   →  002_RMSE.ipynb   →  003_coefs.ipynb
004_fit_MPAA.ipynb  →  005_RMSE_MPAA.ipynb  →  006_coefs_MPAA.ipynb
```

## Features

The following descriptors are computed for each molecule using xTB and morfeus.

### Molecular-level descriptors
| Descriptor | Description |
|------------|-------------|
| HOMO / LUMO | Frontier orbital energies |
| IP / EA | Ionization potential and electron affinity |
| Cone Angle | Palladium only. Steric cone angle |

### Atom-level descriptors (computed at specific labeled atoms)
| Descriptor | Description |
|------------|-------------|
| charge | Partial charge |
| fukui_minus / zero / plus | Fukui functions (nucleophilicity / radical / electrophilicity) |
| polarizability | Atomic polarizability |
| dipole | Atomic dipole moment |
| p_int | Dispersion interaction parameter (morfeus Dispersion) |
| vbur | Fraction of buried volume (%) |
| B1 / B5 / L | Sterimol parameters (minimum width / maximum width / length) |

- **MPAA**: labeled atoms 1, 3, 7, 9 → 49 features in total
- **Pyridone**: labeled atoms 1, 2, 3 → 37 features in total
- **Combined**: 86 features

## Target Variable

| Variable | Description |
|----------|-------------|
| Target | `beta_alpha` = `log(β_av / α_av)` |

## Machine Learning

- **Model**: Ridge regression (features and target standardized with `StandardScaler`)
- **Feature selection**: Recursive Feature Elimination (RFE) — reduces from 86 features down to 4
- **Hyperparameter tuning**: Regularization coefficient α optimized via inner cross-validation
- **Model selection**: The feature count that minimizes the outer CV test RMSE is selected

### Cross-validation settings

| Notebook | Outer CV | Inner CV |
|----------|----------|----------|
| `001_fit` | Leave-One-Pyridone-Out | Leave-One-Pyridone-Out |
| `004_fit_MPAA` | Leave-One-MPAA-Out | Leave-One-MPAA-Out |

## Setup

### Using Docker (recommended)

```bash
docker compose -f docker-compose.dev.yml up
# Access JupyterLab at http://localhost:8888
```

### Local installation

```bash
pip install -r requirements.txt
pip install -e grrmlib/
```

Install xTB (v6.7.1) separately and ensure the `xtb` binary is on your `PATH`.

## Dependencies

| Library | Purpose |
|---------|---------|
| morfeus-ml | Descriptor computation (BuriedVolume, ConeAngle, Dispersion, Sterimol, XTB) |
| polars | DataFrame operations |
| scikit-learn | Ridge regression, StandardScaler, RMSE |
| numpy | Numerical computation |
| matplotlib / seaborn | Visualization |
| grrmlib (local) | Reading Gaussian/GRRM input files |

See `requirements.txt` for exact versions.

## Output Files

| File | Description |
|------|-------------|
| `001_fit_train_val_preds.parquet` | Train/validation predictions (Pyridone CV) |
| `001_fit_test_preds.parquet` | Test predictions (Pyridone CV) |
| `001_fit_coefs.parquet` | Ridge coefficients at each RFE step (Pyridone CV) |
| `004_fit_MPAA_*.parquet` | Same as above (MPAA CV) |
| `002_RMSE_plot.png` | RMSE vs. number of features |
| `002_RMSE_scatter_all.png` | Predicted vs. experimental yield (all data) |
| `002_RMSE_scatter_each.png` | Predicted vs. experimental yield (per pyridone) |
| `003_coefs_imshow.png` | Heatmap of model coefficients |
| `003_coefs_bar.png` | Bar chart of mean ± std coefficients |
| `005_*` / `006_*` | Equivalent plots for MPAA CV |
