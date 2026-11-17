# FTRF-SIR influenza forecasting

This repository contains the code, archived input data, and result figures for the FTRF-SIR influenza forecasting study. The model combines a susceptible–infectious–removed (SIR) forecast with random-forest residual learning and predictions transferred from source regions.

## Contents

- `FTRF_SIR_Rebuilt/code/`: Python data preparation, forecasting, simulation, scoring, statistical comparisons, and plotting.
- `R_scientific_revision/`: the five conditional-regression tasks, additional weekly forecasts, and full-predictor simulations, with Python summary and plotting scripts.
- `ILI/`, `CLs/`, `CPT/`, and `PHLs/`: the original R entry scripts, which call the shared analysis.
- `data/` and `FTRF_SIR_Rebuilt/data/raw/`: the four archived input CSV files.
- [`figures/main/`](figures/main/): six result figures from the main text, Figures 2–7.
- [`figures/supplement/`](figures/supplement/): five supplementary result figures, Figures S1–S5.
- [`FIGURE_SOURCES.md`](FIGURE_SOURCES.md): the data, code, and plotting choices behind each result figure.
- [`FIGURE_PROVENANCE.json`](FIGURE_PROVENANCE.json): file checksums linking the eleven figures to their plotting scripts and ten saved input tables.

The files in `data/` and `FTRF_SIR_Rebuilt/data/raw/` are identical. The second copy is intentional: `data/` keeps the original inputs together, while `FTRF_SIR_Rebuilt/data/raw/` preserves the paths used by the analysis. After cloning the repository, no data files need to be moved.

## Data

The input files are archived exports of aggregate U.S. state-level influenza surveillance data from the [CDC FluView portal](https://gis.cdc.gov/grasp/fluview/fluportaldashboard.html):

- `ILINet.csv`
- `WHO_NREVSS_Clinical_Labs.csv`
- `WHO_NREVSS_Combined_prior_to_2015_16.csv`
- `WHO_NREVSS_Public_Health_Labs.csv`

These files contain surveillance summaries rather than patient-level records. Because FluView data may be revised retrospectively, the copies in this repository should be used when reproducing the reported analysis.

## Python setup

Use Python 3.13 and the pinned dependencies, which record the saved analysis environment. From the repository root:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
```

On Linux or macOS, activate the environment with `source .venv/bin/activate`.

## Main analysis

Run all commands from the repository root. Before the first run, create the output directories used by the scripts:

```powershell
python -c "from pathlib import Path; [Path(p).mkdir(parents=True, exist_ok=True) for p in ['FTRF_SIR_Rebuilt/results/tables','FTRF_SIR_Rebuilt/results/figures','FTRF_SIR_Rebuilt/audit','latex/English/figures','latex/Chinese/figures']]"
```

The `latex/.../figures` directories are image output paths retained by the plotting code; the analysis does not require a TeX installation or manuscript source files.

The principal sample-budget comparison uses:

```powershell
python FTRF_SIR_Rebuilt/code/prepare_data.py
python FTRF_SIR_Rebuilt/code/extension.py --workers 4
python FTRF_SIR_Rebuilt/code/summarize_extension.py
python FTRF_SIR_Rebuilt/code/validate_extension.py
```

The earlier nine-model comparison and the separate parameter-transfer simulation use:

```powershell
python FTRF_SIR_Rebuilt/code/run_forecasts.py
python FTRF_SIR_Rebuilt/code/posterior_experiment.py
python FTRF_SIR_Rebuilt/code/summarize.py
python FTRF_SIR_Rebuilt/code/posterior_rank_diagnostics.py
python FTRF_SIR_Rebuilt/code/validate.py
```

Run `prepare_data.py` before either sequence. The parameter simulation must precede `summarize.py`, which reads both forecast and simulation results. These analyses use different evaluation designs; their results are reported separately.

Predictions and summaries are written to `FTRF_SIR_Rebuilt/results/extension/` for the principal comparison and `FTRF_SIR_Rebuilt/results/tables/` for the earlier comparison and parameter simulation.

## Additional R analysis

R 4.4 or later and the `randomForest` package are required. From the repository root, run:

```powershell
Rscript R_scientific_revision/setup.R
Rscript R_scientific_revision/run.R
Rscript R_scientific_revision/forecast.R
Rscript R_scientific_revision/save_simulation_truth.R
Rscript R_scientific_revision/validate.R
python R_scientific_revision/summarize.py
```

`setup.R` installs `randomForest` in a project-local library. The scripts save results in `R_scientific_revision/results/`. This analysis has different outcomes, sample limits, numerical SIR settings, and comparison models from the principal Python analysis.

## Reproducibility notes

- Model settings and random seeds are defined in the source files.
- Forecasting methods in a given comparison use the same historical cutoffs and eligible target tasks.
- The repository includes the raw snapshots, eleven final result figures, and the ten saved CSV tables read by their plotting scripts. These tables preserve the plotted results without requiring model fitting to obtain them again. Other intermediate files and local software environments remain outside Git.
- Each figure's source script, input tables, and plotting choices are identified in [FIGURE_SOURCES.md](FIGURE_SOURCES.md). [FIGURE_PROVENANCE.json](FIGURE_PROVENANCE.json) records SHA-256 checksums for the raw inputs, saved plotting tables, scripts, and archived figures.
- Text checksums normalize CRLF line endings to LF, allowing the same content to be checked across operating systems. PDF checksums use the original file bytes.
- `python verify_figure_sources.py` checks those archived files using the Python standard library. It reads file contents and CSV headers; it does not execute analysis or plotting scripts.

The saved plotting tables were archived from the completed study. This repository update checked file integrity and table schemas; the numerical experiments and figures were not rerun.

## License

The code is released under the MIT License. The CDC data retain the terms and attribution requirements of their original source.
