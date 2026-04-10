# Experiment — Reaction Yield Data Analysis

This project visualizes and tabulates experimental yield data from a reaction screening campaign across combinations of MPAA catalysts (backbone) and pyridone substrates.

## Overview

- **Raw data**: `combined.xlsx` — yield measurements for each backbone (MPAA) × pyridone combination
- **Products measured**: alpha and beta regioisomers
- **Dataset size**: 8 MPAA variants (backbone 0–5,19,20) × 13 pyridone variants (pyridone 0–12)

## Data Structure (`combined.xlsx`)

| Column | Description |
|--------|-------------|
| `backbone` | MPAA catalyst index (0–5,19,20) |
| `pyridone` | Pyridone substrate index (0–12) |
| `beta1` | Beta product yield, run 1 (%) |
| `beta2` | Beta product yield, run 2 (%) |
| `beta_av` | Beta product average yield (%) |
| `alpha1` | Alpha product yield, run 1 (%) |
| `alpha2` | Alpha product yield, run 2 (%) |
| `alpha_av` | Alpha product average yield (%) |

## Notebooks

### `001_yield.ipynb` — Heatmaps and Histograms

Visualizes the yield landscape across the MPAA × pyridone matrix.

**Outputs:**

| File | Description |
|------|-------------|
| `001_yield_imshow_beta.png` | Heatmap of beta product average yield (MPAA × pyridone) |
| `001_yield_imshow_alpha.png` | Heatmap of alpha product average yield (MPAA × pyridone) |
| `001_yield_imshow_selectivity.png` | Heatmap of alpha/beta selectivity ratio (MPAA × pyridone) |
| `001_yield_hist_pyridone.png` | Stacked histogram of beta yield distribution, colored by pyridone |
| `001_yield_hist_MPAA.png` | Stacked histogram of beta yield distribution, colored by MPAA |

All heatmaps display numerical values in each cell, with a color scale from jet colormap.

### `002_yield.ipynb` — Summary Tables

Generates formatted yield summary tables for publication or reporting.

- MPAA indices are remapped (19 → 6, 20 → 7) for display.
- All yield values are formatted to one decimal place.
- Tables are split by pyridone group for readability.

**Outputs:**

| File | Description |
|------|-------------|
| `002_yield_table.png` | Full summary table (all pyridones) |
| `002_yield_table0-4.png` | Summary table for pyridones 0–4 |
| `002_yield_table5-9.png` | Summary table for pyridones 5–9 |
| `002_yield_table10-12.png` | Summary table for pyridones 10–12 |

## Project Structure

```
experiment/
├── combined.xlsx               # Raw yield data
├── 001_yield.ipynb             # Heatmaps and histograms
├── 001_yield_imshow_beta.png
├── 001_yield_imshow_alpha.png
├── 001_yield_imshow_selectivity.png
├── 001_yield_hist_pyridone.png
├── 001_yield_hist_MPAA.png
├── 002_yield.ipynb             # Summary tables
├── 002_yield_table.png
├── 002_yield_table0-4.png
├── 002_yield_table5-9.png
├── 002_yield_table10-12.png
├── requirements.txt
├── Dockerfile
└── docker-compose.dev.yml
```

## Environment

Python 3.10 with the following key packages:

- `polars` 1.19.0
- `matplotlib` 3.10.7
- `numpy` 1.24.2
- `fastexcel` 0.13.0 (for reading `.xlsx` files with Polars)
- `jupyterlab` 4.0.10

### Running with Docker

```bash
docker compose -f docker-compose.dev.yml up
```

JupyterLab will be available at `http://localhost:8888`.
