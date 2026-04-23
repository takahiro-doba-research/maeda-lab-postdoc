# Selectivity Prediction

A linear regression pipeline for predicting chemical reaction selectivity (β/α ratio) using quantum-chemical energy descriptors and Ridge regression with Recursive Feature Elimination (RFE).

## Overview

This project predicts the logarithmic β/α selectivity (`log(β_av / α_av)`) of reactions between pyridone variants and backbone (MPAA ligand) variants. Features are energy differences between α- and β-intermediate states computed by quantum chemistry calculations. Two cross-validation scenarios are compared to assess generalization to unseen substrates.

## Dataset

| Variable | Count | Description |
|---|---|---|
| Pyridone | 13 (IDs 0–12) | Pyridone substrate variants |
| Backbone | 8 (IDs 0–5, 19, 20) | Chiral backbone (MPAA ligand) variants |
| Samples | 104 | All pyridone–backbone combinations |
| Features | 40 | Energy differences (α − β intermediates) at each reaction-path position |
| Target | `beta_alpha` | `log(β_av / α_av)` |

**Input files:**
- `energy_descriptor_alpha.csv` — energies of α-pathway intermediates
- `energy_descriptor_beta.csv` — energies of β-pathway intermediates
- `combined.xlsx` — experimentally measured α and β selectivity values

## Notebooks

| Notebook | Description |
|---|---|
| `000_dataset.ipynb` | Constructs the dataset: computes α − β energy differences and joins with experimental selectivity |
| `001_fit.ipynb` | Ridge + RFE with **leave-one-pyridone-out** nested cross-validation |
| `002_RMSE.ipynb` | RMSE vs. number of features curve and scatter plots (pyridone-out CV) |
| `003_coefs.ipynb` | Ridge coefficient analysis (pyridone-out CV) |
| `004_fit_MPAA.ipynb` | Ridge + RFE with **leave-one-MPAA-out** nested cross-validation |
| `005_RMSE_MPAA.ipynb` | RMSE vs. number of features curve and scatter plots (MPAA-out CV) |
| `006_coefs_MPAA.ipynb` | Ridge coefficient analysis (MPAA-out CV) |

## Methodology

### Feature engineering

Energy descriptors are computed as the difference between the α-pathway and β-pathway intermediate energies:

```
ΔE = E_alpha − E_beta
```

This yields 40 scalar features per pyridone–MPAA pair.

### Model

**Ridge regression** (`sklearn.linear_model.Ridge`) with feature standardization (`StandardScaler` applied to both X and y).

**Recursive Feature Elimination (RFE):** Features are eliminated one at a time by removing the feature with the smallest absolute Ridge coefficient, from 40 down to 4.

**Nested cross-validation:** Hyperparameter `alpha` is tuned by an inner leave-one-group-out loop at each RFE step; the outer loop evaluates generalization on the held-out group.

## Outputs

| File | Description |
|---|---|
| `dataset_selectivity.parquet` | Merged feature–target dataset |
| `001_fit_train_val_preds.parquet` | Train/val predictions for all pyridone folds and RFE steps |
| `001_fit_test_preds.parquet` | Test predictions for pyridone-out CV |
| `001_fit_coefs.parquet` | Ridge coefficients for pyridone-out CV |
| `004_fit_MPAA_*.parquet` | Corresponding outputs for MPAA-out CV |
| `002_RMSE_plot.png` | RMSE vs. n_features (pyridone-out) |
| `002_RMSE_scatter_all.png` | Predicted vs. experimental scatter (pyridone-out, all data) |
| `002_RMSE_scatter_each.png` | Per-pyridone scatter plots |
| `003_coefs_imshow.png` | Heatmap of coefficients per pyridone test fold |
| `003_coefs_bar.png` | Mean ± std of coefficients across pyridone folds |
| `005_RMSE_MPAA_*.png` | Corresponding figures for MPAA-out CV |
| `006_coefs_MPAA_*.png` | Coefficient figures for MPAA-out CV |

## Setup

### With Docker (recommended)

```bash
docker compose -f docker-compose.dev.yml up
```

JupyterLab will be available at `http://localhost:8888`.

### Without Docker

```bash
pip install -r requirements.txt
jupyter lab
```

**Python version:** 3.11

**Key dependencies:**

| Package | Version |
|---|---|
| polars | 1.19.0 |
| scikit-learn | 1.5.2 |
| numpy | 1.24.2 |
| matplotlib | 3.10.7 |
| seaborn | 0.13.2 |
| pandas | 1.5.3 |
| fastexcel | 0.13.0 |

## Execution order

Run notebooks in numerical order:

```
000_dataset.ipynb
001_fit.ipynb
002_RMSE.ipynb
003_coefs.ipynb
004_fit_MPAA.ipynb
005_RMSE_MPAA.ipynb
006_coefs_MPAA.ipynb
```
