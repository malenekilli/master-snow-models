# Master Snow Models: Evaluation in Complex Alpine Terrain

Master thesis – Snow model comparison (FSM2, seNorge, SNOWPACK) in the Jotunheimen region, southern Norway

## Abstract

Seasonal snow is a key component of the hydrological and hazard regimes in Norwegian mountain regions, yet its representation in numerical models remains challenging due to complex energy-balance processes, sparse forcing data and strong spatial variability. This thesis investigates how three contrasting snow model paradigms represent snow depth in complex alpine terrain and how model evaluation depends on reference choice and spatial scale.

The study focuses on the Jotunheimen region in southern Norway during the 2024–2025 snow season and uses **SNOWPACK** as a physics-based process benchmark alongside two operationally relevant gridded products, **FSM2** and **seNorge**.

### Key Findings

- **All three models** reproduced the overall seasonal evolution of snow depth, including main accumulation and melt phases
- **FSM2** generally exhibited lower RMSE than seNorge when compared to SNOWPACK, with closer agreement in temporal evolution
- **seNorge** achieved the lowest RMSE when validated against independent station observations
- **Model ranking is not absolute** but depends critically on the chosen reference dataset, spatial scale and terrain context
- **Largest discrepancies** occurred during the melt season and at high elevations in western parts of the domain

## Project Structure

```
├── README.md                          # This file
├── data/                              # Raw and processed datasets
│   ├── snowdepth.csv                  # Station observation reference data
│   ├── FSM2/                          # Flexible Snow Model v2 outputs
│   ├── seNorge/                       # Norwegian seNorge snow products
│   ├── SNOWPACK/                      # SNOWPACK model simulations
│   ├── GIS/                           # Geospatial data (regions, pixels)
│   └── processed/                     # Processed datasets and results
├── notebooks/                         # Jupyter notebooks for analysis
│   ├── build_nc_files/                # Dataset compilation scripts
│   ├── 01-08_*                        # Data exploration and preprocessing
│   ├── 09-11_*                        # Validation and comparison analysis
│   ├── 12-18_*                        # Thesis visualization and results
│   └── README.md                      # Notebook guide
├── scripts/                           # R and Python analysis scripts
├── figures_thesis/                    # Generated figures for thesis
└── build_nc_files/                    # NetCDF dataset construction
```

## Datasets

### Input Data

- **FSM2** (Flexible Snow Model v2): Grid-based snow depth (SD) and snow water equivalent (SWE), 2022–2025
- **seNorge**: Operational snow depth and SWE products, 2022–2025
- **SNOWPACK**: Point simulations at 4,484 virtual station locations converted to gridded format
- **Station Observations**: Daily snow depth from 11 meteorological stations (2022–2025)
- **GIS Data**: Forecasting regions (WGS84), pixel coordinates (UTM)

### Processed Data

All key datasets are stored in `data/processed/`:

- `fsm2_sd_all.nc` / `senorge_all.nc` / `snowpack_all.nc` – Complete model datasets
- `*_validation_*.csv` – Comparison metrics and validation results
- `seasonal_*.csv` – Seasonal aggregations and analysis
- `validation_stations_spatial.geojson` – Station locations with spatial context
- `model_comparison_metrics.csv` – Summary statistics across all models

## Analysis Overview

### 1. Data Compilation (`build_nc_files/`)
- Merge annual FSM2, seNorge, and SNOWPACK outputs into continuous NetCDF files
- Standardize grids and coordinate systems

### 2. Exploratory Analysis (`01_-_08_*.ipynb`)
- Examine temporal and spatial patterns in each model
- Clip domain to Jotunheimen region
- Create baseline statistics and seasonal summaries

### 3. Model-to-Model Comparison (`04_-_08_*.ipynb`)
- Compare FSM2 and seNorge directly against SNOWPACK benchmark
- Calculate bias, RMSE, and temporal correlation
- Analyze seasonal and spatial patterns of disagreement

### 4. Observation-Based Validation (`09_-_11_*.ipynb`)
- Validate all three models against independent station observations
- Assess influence of terrain complexity and representativeness
- Multi-station validation across elevation and terrain types

### 5. Thesis Visualization (`12_-_18_*.ipynb`)
- Create publication-quality figures for thesis
- Spatial maps and distribution plots
- Model comparison summaries

## Methods

### Evaluation Metrics
- **Bias**: Mean difference between model and reference
- **RMSE**: Root mean square error (primary metric)
- **Temporal Correlation**: Pearson correlation of time series
- **Median Error**: Robust measure less sensitive to outliers

### Spatial Scale
- Grid cell resolution: varies by model (FSM2: ~1 km, seNorge: ~1 km, SNOWPACK: 4,484 virtual stations gridded)
- Domain: Jotunheimen region, southern Norway
- Elevation range: ~800–2,400 m

### Reference Datasets
1. **SNOWPACK**: Physics-based energy-balance model (validation baseline)
2. **Station Observations**: Independent measurements from 11 meteorological stations
3. **Multi-station Metrics**: Cross-validation to examine representativeness

## Getting Started

### Requirements
- Python 3.8+
- Jupyter Notebook or JupyterLab
- Key packages: `numpy`, `xarray`, `pandas`, `netCDF4`, `matplotlib`, `cartopy`, `geopandas`, `scipy`

### Installation

```bash
# Clone repository
git clone <repository-url>
cd Master_snow_models

# Install dependencies (create virtual environment first)
python -m venv venv
source venv/bin/activate  # or 'venv\Scripts\activate' on Windows
pip install -r requirements.txt

# Launch Jupyter
jupyter notebook
```

### Quick Start

1. Start with `notebooks/01_explore_fsm2.ipynb` and `01_explore_seNorge.ipynb` to understand the datasets
2. Review `notebooks/04_compare_models.ipynb` for model comparison methodology
3. Check `notebooks/10_multi_station_validation.ipynb` for validation results
4. See `notebooks/14_seasonal_model_comparison.ipynb` for seasonal analysis
5. Explore `notebooks/16_thesis_result_plots.ipynb` for final results

## Key Results Files

| File | Description |
|------|-------------|
| `observation_validation_summary.csv` | Metrics for all three models vs. station observations |
| `model_comparison_metrics.csv` | FSM2 vs. seNorge comparison statistics |
| `seasonal_model_metrics.csv` | Seasonal breakdown of model performance |
| `multi_station_validation_results.csv` | Station-by-station validation metrics |
| `station_distance_statistics.csv` | Analysis of model-station representativeness |

## Figures

Generated figures are saved in `figures_thesis/`:
- Model comparison maps and time series
- Spatial error distributions
- Seasonal evolution plots
- Validation summary figures

## Technical Notes

- **NetCDF Format**: All gridded models stored in NetCDF4 format with standard CF conventions
- **Coordinate System**: UTM for pixel data, WGS84 for GIS overlays
- **Temporal Resolution**: Daily data for all models and observations
- **Missing Data**: Handled with xarray masking; see individual notebooks for specifics

## References

- **FSM2**: Essery et al. (2013) – Flexible Snow Model
- **seNorge**: Norwegian Centre for Water Resources and Energy (NVE) operational product
- **SNOWPACK**: Lehning et al. (1999, 2002) – Physics-based snow cover model
- **Study Region**: Jotunheimen National Park, southern Norway (61°N, 8°E)

## Author

*Master Thesis, University of NTNU*  
*Malene Lien Killi*  
*June 2026*



**Questions or suggestions?** Please open an issue or contact the author.
