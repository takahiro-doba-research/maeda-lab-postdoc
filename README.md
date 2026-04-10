# maeda-lab-postdoc

## Description

This repository contains codes and data for the paper:  
Doba, T.\*; Harabuchi, Y.; Nagata, Y.; Maeda, S.\* Construction of an Interpretable Regression Model for Yield Prediction and Mechanistic Elucidation Enabled by Automated Reaction Path Exploration. **2026**

The repository covers the full data pipeline for Pd-catalyzed C–H functionalization, from quantum chemical calculations and experimental screening to machine learning–based yield prediction. The code has been modified from its initial version, while maintaining the core functionality. It has only been tested in the computer environment listed below and may require minor adjustments to run on other systems.

## Author

Takahiro Doba (Kyoto University)

## Acknowledgments

This study was supported by JST-ERATO (grant number JPMJER1903) to S.M., JSPS KAKENHI (grant number 22KJ0043) to T.D., and JSPS KAKENHI (grant number 23B203) to Y.N. and Y.H.

---

## Repository Structure

```
.
├── calculation/       # SC-AFIR reaction path search and energy descriptor generation
├── experiment/        # Experimental yield data visualization and tabulation
└── machine-learning/  # Linear regression models for yield prediction
```

---

## `calculation/` — Quantum Chemical Calculations

SC-AFIR (Stochastic Conformer-based AFIR) calculations using GRRM/ONIOM(B3LYP/LanL2DZ:UFF) to map the reaction path network of a Pd-catalyzed C–H functionalization reaction. The reaction involves a pyridone directing group and an MPAA (mono-protected amino acid) ligand.

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
3. **Repath** (`td_c1_193_Repath1.com`): refines transition state paths from the SCAFIR6 results.
4. **Fragment separation**: splits each intermediate into fragments and reoptimizes (B3LYP/Def2SVP) with `jobcontroller.SeparationJob`.
5. **Substituent connection**: connects 8 MPAA variants (`backbone/`) and 13 pyridone variants (`pyridone/`) combinatorially and reoptimizes with `jobcontroller.ConnectionJob`.
6. **Energy descriptors**: SCF energies for all 104 backbone–pyridone combinations across 52 intermediates are exported to `energy_descriptor.csv` (= `X.csv` used in machine learning).

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

## `machine-learning/` — Yield Prediction Models

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

All models share the same framework: recursive feature elimination (RFE, from 52 down to 4 features) with nested leave-one-pyridone-out cross-validation. Ridge regression receives additional analysis including alternative data splits, a shuffled-label negative control, feature correlation analysis, and noise robustness testing.

See [`machine-learning/README.md`](machine-learning/README.md) for details.

---

## Data Flow

```
calculation/
  SC-AFIR search
      ↓
  Fragment separation & substituent connection
      ↓
  energy_descriptor.csv ──────────────────────────────→ machine-learning/X.csv
                                                                  ↓
experiment/                                              Linear regression
  combined.xlsx ──────────────────────────────────────→ machine-learning/y.csv
  (beta_av yield)                                                  ↓
                                                         Yield prediction model
```

## Environment

Each subdirectory has its own `Dockerfile` and `docker-compose.dev.yml`. All three environments run JupyterLab on `http://localhost:8888` via:

```bash
docker compose -f docker-compose.dev.yml up
```

All environments use Python 3.10–3.11 with `polars`, `numpy`, `matplotlib`, and `jupyterlab` as common dependencies. The `machine-learning/` environment additionally requires `scikit-learn`; the `calculation/` environment requires `grrmlib`, `networkx`, and `cclib`, as well as external **GRRM** and **Gaussian** installations.
