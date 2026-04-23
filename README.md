# maeda-lab-postdoc

## Description

This repository contains codes and data for the paper:  
Doba, T.\*; Harabuchi, Y.; Nagata, Y.; Maeda, S.\* Construction of an Interpretable Regression Model for Yield Prediction and Mechanistic Insight Enabled by Automated Reaction Path Exploration. **2026**

The repository covers the full data pipeline for Pd-catalyzed C–H functionalization, from quantum chemical calculations and experimental screening to linear regression–based yield prediction. The code has been modified from its initial version, while maintaining the core functionality. It has only been tested in the computer environment listed below and may require minor adjustments to run on other systems.

## Author

Takahiro Doba (Kyoto University)

## Acknowledgments

This study was supported by JST-ERATO (grant number JPMJER1903) to S.M., JSPS KAKENHI (grant number 22KJ0043) to T.D., and JSPS KAKENHI (grant number 23B203) to Y.N. and Y.H.

---

## Repository Structure

```
.
├── calculation/                # SC-AFIR reaction path search and energy descriptor generation
├── experiment/                 # Experimental yield data visualization and tabulation
├── yield_prediction/           # Linear regression models for yield prediction
├── selectivity_prediction/     # Regioselectivity (β/α) prediction models
└── atom_based_featurization/   # Atom-based descriptor featurization and yield prediction
```

---

## `calculation/` — Quantum Chemical Calculations

SC-AFIR calculations using GRRM/ONIOM(B3LYP/LanL2DZ:UFF) to map the reaction path network of a Pd-catalyzed C–H functionalization reaction. The reaction involves a pyridone and an MPAA (mono-protected amino acid) as ligands.

### Reaction pathway (6 key intermediates)

| Step | Intermediate |
|------|-------------|
| 0 | Benzene complex |
| 1 | C–H activation complex |
| 2 | Olefin coordination complex |
| 3 | Proton transfer complex |
| 4 | Carbopalladation complex |
| 5 | Product complex |

### Workflow

1. **SC-AFIR search** (6 sequential runs, `td_c1_193_SCAFIR*.com`): explores the potential energy surface and produces lists of equilibrium structures (`*_EQ_list.log`) and approximate transition states (`*_PT_list.log`).
2. **Reaction path analysis** (`td_c1_193_SCAFIR*.ipynb`): filters EQs by interatomic distances using `grrmlib` to identify each intermediate; SCAFIR6 notebook builds the full reaction path network and extracts all key EQ/PT structures.
3. **Repath** (`td_c1_193_Repath1.com`): refines reaction paths from the SCAFIR6 results.
4. **Fragment separation**: splits each intermediate into fragments and reoptimizes (B3LYP/Def2SVP) with `jobcontroller.SeparationJob`. These functions are also implemented in `grrmlib.SeparationWorkflow`.
5. **Substituent connection**: connects 8 MPAA variants (`backbone/`) and 13 pyridone variants (`pyridone/`) combinatorially and reoptimizes with `jobcontroller.ConnectionJob`. These functions are also implemented in `grrmlib.ConnectionWorkflow`.
6. **Energy descriptors**: SCF energies for all 104 backbone–pyridone combinations across 52 intermediates are exported to `energy_descriptor.csv` (= `X.csv` used in linear regression).

See [`calculation/README.md`](calculation/README.md) for details.

---

## `experiment/` — Experimental Yield Data

Experimental yield measurements for 8 MPAA variants (backbone 0–5, 19, 20) × 13 pyridone variants (pyridone 0–12), producing both alpha and beta regioisomers (each measured in duplicate).

### Contents

| Notebook | Output |
|----------|--------|
| `001_yield.ipynb` | Heatmaps of beta yield, alpha yield, and alpha/beta selectivity; stacked yield histograms by pyridone and MPAA |
| `002_yield.ipynb` | Formatted yield summary tables (split by pyridone groups 0–4, 5–9, 10–12) |

Raw data: `combined.xlsx` (`backbone`, `pyridone`, `beta1`, `beta2`, `beta_av`, `alpha1`, `alpha2`, `alpha_av`).

See [`experiment/README.md`](experiment/README.md) for details.

---

## `yield_prediction/` — Yield Prediction Models

Linear regression models that predict the beta product yield (`beta_av`) from quantum-chemically computed reaction intermediate energies.

- **Features (X.csv)**: absolute SCF energies of 52 intermediates → converted to relative energies (each intermediate − reference state)
- **Target (y.csv)**: `beta_av` (%) → logit-transformed for regression
- **Dataset**: 104 samples (8 backbone × 13 pyridone combinations)

### Models compared

| Model | Key hyperparameter |
|-------|--------------------|
| OLS (Ordinary Least Squares) | — |
| Ridge Regression | `alpha` |
| Lasso Regression | `alpha` |
| PLS (Partial Least Squares) | `n_components` |

All models share the same framework: recursive feature elimination (RFE, from 52 down to 4 features) with nested leave-one-group-out cross-validation (LOGOCV). Ridge regression receives additional analysis including alternative data splits, a shuffled-label negative control, feature correlation analysis, and noise robustness testing.

See [`yield_prediction/README.md`](yield_prediction/README.md) for details.

---

## `selectivity_prediction/` — Regioselectivity Prediction Models

Ridge regression pipeline for predicting the logarithmic β/α regioselectivity (`log(β_av / α_av)`) of Pd-catalyzed C–H functionalization reactions.

- **Features**: 40 energy-difference descriptors per sample, computed as `ΔE = E_alpha − E_beta` (difference between α- and β-pathway intermediate energies)
- **Target**: `log(β_av / α_av)` (log ratio of experimentally measured β and α yields)
- **Dataset**: 104 samples (8 MPAA × 13 pyridone combinations)
- **Input files**: `energy_descriptor_alpha.csv`, `energy_descriptor_beta.csv`, `combined.xlsx`

### Methodology

Ridge regression with `StandardScaler` normalization and **Recursive Feature Elimination (RFE)** (from 40 down to 4 features). Two nested LOGOCV scenarios are compared:

| Notebooks | Cross-validation | Description |
|-----------|-----------------|-------------|
| `001_fit` / `002_RMSE` / `003_coefs` | Leave-one-pyridone-out | Generalization to unseen pyridone substrates |
| `004_fit_MPAA` / `005_RMSE_MPAA` / `006_coefs_MPAA` | Leave-one-MPAA-out | Generalization to unseen MPAA ligands |

Run notebooks in numerical order starting from `000_dataset.ipynb`.

See [`selectivity_prediction/README.md`](selectivity_prediction/README.md) for details.

---

## `atom_based_featurization/` — Atom-Based Descriptor Featurization

Linear regression pipeline for yield prediction (`beta_av`) using atom-level and molecular-level descriptors computed by semi-empirical quantum chemistry (xTB) and the morfeus library, as an alternative to the reaction-intermediate energy descriptors used in `yield_prediction/`.

- **Features**: 86 descriptors per sample — 49 from MPAA (labeled atoms 1, 3, 7, 9) and 37 from pyridone (labeled atoms 1, 2, 3)
- **Target**: `beta_av_logit` (logit-transformed beta product yield)
- **Dataset**: 104 samples (8 MPAA × 13 pyridone combinations)

### Feature types

| Category | Descriptors |
|----------|-------------|
| Molecular-level | HOMO/LUMO energies, ionization potential (IP), electron affinity (EA), cone angle (Pd only) |
| Atom-level | partial charge, Fukui functions (nucleophilic/radical/electrophilic), polarizability, dipole, dispersion parameter (`p_int`), buried volume (`vbur`), Sterimol parameters (B1, B5, L) |

### Methodology

Ridge regression with `StandardScaler` normalization and **RFE** (from 86 down to 4 features). Hyperparameter `alpha` is tuned via an inner LOGOCV loop. Two CV settings are evaluated:

| Notebooks | Outer CV | Description |
|-----------|----------|-------------|
| `001_fit` / `002_RMSE` / `003_coefs` | Leave-one-pyridone-out | Generalization to unseen pyridone substrates |
| `004_fit_MPAA` / `005_RMSE_MPAA` / `006_coefs_MPAA` | Leave-one-MPAA-out | Generalization to unseen MPAA ligands |

Run notebooks in order: `000_dataset.ipynb` → `001_fit` → `002_RMSE` → `003_coefs` → `004_fit_MPAA` → `005_RMSE_MPAA` → `006_coefs_MPAA`.

> **Note**: Requires xTB v6.7.1 installed separately with the `xtb` binary on `PATH`.

See [`atom_based_featurization/README.md`](atom_based_featurization/README.md) for details.

---

## Environment

Each subdirectory has its own `Dockerfile` and `docker-compose.dev.yml`. All environments run JupyterLab on `http://localhost:8888` via:

```bash
docker compose -f docker-compose.dev.yml up
```

All environments use Python 3.10–3.11 with `polars`, `numpy`, `matplotlib`, and `jupyterlab` as common dependencies. Additional requirements per module:

| Module | Extra dependencies |
|--------|--------------------|
| `yield_prediction/` | `scikit-learn`, `pandas`, `seaborn` |
| `selectivity_prediction/` | `scikit-learn`, `pandas`, `seaborn`, `fastexcel` |
| `atom_based_featurization/` | `scikit-learn`, `morfeus-ml`, `grrmlib` (local); **xTB v6.7.1** (external binary) |
| `calculation/` | `grrmlib`, `networkx`, `cclib`; external **GRRM** and **Gaussian** installations |
