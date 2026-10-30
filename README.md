# FTRF-SIR influenza forecasting

This repository contains the code and archived input data used for the FTRF-SIR influenza forecasting study. The analysis combines forecasts from a susceptible--infectious--removed model with random-forest residual learning and information transferred from source regions.

## Contents

- `FTRF_SIR_Rebuilt/code/` contains the main Python analysis.
- `FTRF_SIR_Rebuilt/data/raw/` contains the four CSV files read by the Python and R scripts.
- `data/` contains a second, central copy of the same four input files.
- `R_scientific_revision/` contains the additional R analysis and its Python summary scripts.
- `ILI/`, `CLs/`, `CPT/`, and `PHLs/` contain the R entry scripts for the four surveillance outcomes.

The files in `data/` and `FTRF_SIR_Rebuilt/data/raw/` are identical. The second copy is intentional: `data/` keeps the original inputs together, while `FTRF_SIR_Rebuilt/data/raw/` preserves the paths used by the analysis. After cloning the repository, no data files need to be moved.

## Data

The input files are archived exports of aggregate U.S. state-level influenza surveillance data from the [CDC FluView portal](https://gis.cdc.gov/grasp/fluview/fluportaldashboard.html):

- `ILINet.csv`
- `WHO_NREVSS_Clinical_Labs.csv`
- `WHO_NREVSS_Combined_prior_to_2015_16.csv`
- `WHO_NREVSS_Public_Health_Labs.csv`

These files contain surveillance summaries rather than patient-level records. Because FluView data may be revised retrospectively, the copies in this repository should be used when reproducing the reported analysis.

## Python setup

Python 3.10 or later is recommended. From the repository root:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

On Linux or macOS, activate the environment with `source .venv/bin/activate`.

## Main analysis

Run the Python scripts from `FTRF_SIR_Rebuilt/code/` in the following order:

```powershell
python prepare_data.py
python run_forecasts.py
python summarize.py
python extension.py
python summarize_extension.py
python plot_extension.py
python validate.py
python validate_extension.py
```

The scripts write predictions, summary tables, and figures to `FTRF_SIR_Rebuilt/results/`. Plotting scripts read the saved tables, so figures can be redrawn without fitting the forecasting models again.

## Additional R analysis

R 4.4 or later and the `randomForest` package are required. From `R_scientific_revision/`, run:

```powershell
Rscript setup.R
Rscript run.R
Rscript forecast.R
Rscript validate.R
python summarize.py
python plot_paper.py
```

This R analysis is a separate sensitivity analysis. Its sample limits, numerical SIR settings, optimization settings, and comparison models differ from those used in the main Python analysis.

## Reproducibility notes

- Model settings and random seeds are defined in the source files.
- Forecasting methods in a given comparison use the same historical cutoffs and eligible target tasks.
- Generated results and local software environments are not stored in Git; they are recreated by the scripts above.

## License

The code is released under the MIT License. The CDC data retain the terms and attribution requirements of their original source.
