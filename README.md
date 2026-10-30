# FTRF-SIR influenza forecasting code

This repository contains the analysis code for the manuscript **"Integrating Transfer Learning with Epidemic Dynamics: A Hybrid FTRF-SIR Framework for Robust Infectious Disease Forecasting."** The implementation combines susceptible--infectious--removed forecasts, random forest residual correction, and predictions transferred from source regions.

## Repository structure

- `FTRF_SIR_Rebuilt/code/`: primary Python forecasting, scoring, validation, summarization, and plotting scripts.
- `FTRF_SIR_Rebuilt/data/raw/`: original CDC CSV files at the paths used directly by the Python and R analyses.
- `data/`: a central copy of the complete original-data set for convenient inspection and archiving.
- `R_scientific_revision/`: additional R analyses and their Python summarization and plotting scripts.
- `ILI/`, `CLs/`, `CPT/`, and `PHLs/`: R entry scripts for the corresponding surveillance outcomes.

The repository includes the source code and the archived input data used in the study. The same four original CSV files are stored in both data locations listed above so that the analysis scripts run without path changes while the complete input set remains available in one central directory.

## Data preparation

The analyses use aggregate U.S. state-level influenza surveillance records from the [CDC FluView portal](https://gis.cdc.gov/grasp/fluview/fluportaldashboard.html). No patient-level data are used. The archived CSV files required by the scripts are already included in `FTRF_SIR_Rebuilt/data/raw/`:

- `ILINet.csv`
- `WHO_NREVSS_Clinical_Labs.csv`
- `WHO_NREVSS_Combined_prior_to_2015_16.csv`
- `WHO_NREVSS_Public_Health_Labs.csv`

Identical copies are provided in the top-level `data/` directory. CDC surveillance records may be revised retrospectively; the included snapshots preserve the inputs used for the reported analysis.

## Python environment

Python 3.10 or later is recommended. From the repository root, create an isolated environment and install the required packages:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

On Linux or macOS, activate the environment with `source .venv/bin/activate`.

## Primary Python analysis

Run the following commands from `FTRF_SIR_Rebuilt/code/`:

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

The scripts use fixed random seeds and historical response cutoffs. Forecasting methods within each comparison are evaluated on the same eligible target tasks. Generated predictions, summaries, and figures are written under `FTRF_SIR_Rebuilt/results/`.

`plot_results.py` and `plot_extension.py` read saved result tables, so figures can be regenerated without refitting the forecasting models.

## Additional R analysis

R 4.4 or later is recommended. The R analysis requires the `randomForest` package. From `R_scientific_revision/`, run:

```powershell
Rscript setup.R
Rscript run.R
Rscript forecast.R
Rscript validate.R
python summarize.py
python plot_paper.py
```

This is a separate sensitivity implementation. Its sample limits, numerical SIR settings, optimization settings, and comparison models differ from those of the primary Python benchmark.

## Reproducibility notes

- Models use only outcomes available by the stated historical fitting cutoff.
- Random seeds and model settings are defined in the source files.
- Source-code comments are provided in English.
- The validation scripts check numerical refinement, information cutoffs, common evaluation tasks, and saved-output consistency.

## License

The analysis code is released under the MIT License. See `LICENSE`. Data remain subject to the terms and attribution requirements of their original provider.
