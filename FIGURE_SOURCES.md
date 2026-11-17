# Sources of the manuscript result figures

The six main-text result figures are in `figures/main/`, and the five supplementary figures are in `figures/supplement/`. The CSV paths below identify outputs produced by the analysis scripts. The repository includes the raw inputs and final figure PDFs; intermediate result tables are generated during analysis.

RMSE denotes root mean squared error, and WIS denotes weighted interval score. The plotting code reads the calculated scores or saved predictions identified below.

## Figure 2: Model comparisons

- Figure: [`ext_ablation.pdf`](figures/main/ext_ablation.pdf)
- Plotting code: [`plot_extension.py`](FTRF_SIR_Rebuilt/code/plot_extension.py)
- Input: `FTRF_SIR_Rebuilt/results/extension/paired.csv`

The plot selects WIS comparisons against the target random forest with 52 target training samples and Gaussian residual intervals. Each panel shows model-minus-reference WIS differences for one dataset. Error bars are descriptive 95% percentile intervals from 5,000 paired state–season bootstrap resamples, calculated in `summarize_extension.py`.

## Figure 3: Learning curves

- Figure: [`ext_learning.pdf`](figures/main/ext_learning.pdf)
- Plotting code: [`plot_extension.py`](FTRF_SIR_Rebuilt/code/plot_extension.py)
- Input: `FTRF_SIR_Rebuilt/results/extension/summary.csv`

For each dataset, the plot connects scores at 24, 52, and 104 target training samples for five models. The upper row shows RMSE and the lower row shows Gaussian-interval WIS. The lines connect aggregate scores and do not represent uncertainty bands.

## Figure 4: Interval coverage and width

- Figure: [`ext_intervals.pdf`](figures/main/ext_intervals.pdf)
- Plotting code: [`plot_extension.py`](FTRF_SIR_Rebuilt/code/plot_extension.py)
- Input: `FTRF_SIR_Rebuilt/results/extension/summary.csv`

At 52 target training samples, mean 95% interval width is plotted against empirical coverage for the target random forest, local SIR residual model, and FTRF-SIR. Filled circles represent Gaussian intervals and open squares represent order-statistic intervals. Lines join the two rules for the same model; the dashed line marks nominal 95% coverage.

## Figure 5: Five conditional-regression tasks

- Figure: [`original_five_tasks.pdf`](figures/main/original_five_tasks.pdf)
- Plotting code: [`plot_paper.py`](R_scientific_revision/plot_paper.py)
- Inputs: `R_scientific_revision/results/conditional_block_errors.csv` and `conditional_predictions.csv`

The boxplots summarize test-block mean squared errors by task and model across target regions and chronological cutoffs. Model order is read from the saved predictions. Boxes show the interquartile range, the centre line is the median, and whiskers use the 1.5-interquartile-range rule; outlier markers are omitted. Units are squared percentage points for ILI and clinical influenza B positivity, and squared counts for the remaining tasks.

## Figure 6: Additional weekly forecasts

- Figure: [`forecast_comparison.pdf`](figures/main/forecast_comparison.pdf)
- Plotting code: [`plot_paper.py`](R_scientific_revision/plot_paper.py)
- Input: `R_scientific_revision/results/forecast_summary.csv`

`R_scientific_revision/summarize.py` calculates RMSE from common evaluation rows in `forecast_predictions.csv`. The three panels show influenza-like illness, clinical influenza B positivity, and historical combined-laboratory influenza B positivity. Within each panel, horizontal bars are ordered by increasing RMSE, measured in percentage points.

## Figure 7: Full-predictor simulations

- Figure: [`simulation_comparison.pdf`](figures/main/simulation_comparison.pdf)
- Plotting code: [`plot_paper.py`](R_scientific_revision/plot_paper.py)
- Input: `R_scientific_revision/results/simulation_summary.csv`

`forecast.R` generates the simulated data and predictions; `summarize.py` calculates the aggregate scores. The two panels show matched and shifted source conditions, with models ordered by increasing RMSE in percentage points. The bars show aggregate RMSE; paired bootstrap intervals are reported separately in the supplementary tables.

## Figure S1: Forecast horizons

- Figure: [`ext_horizons.pdf`](figures/supplement/ext_horizons.pdf)
- Plotting code: [`plot_extension.py`](FTRF_SIR_Rebuilt/code/plot_extension.py)
- Input: `FTRF_SIR_Rebuilt/results/extension/horizons.csv`

The plot selects 52 target training samples and Gaussian records, then connects RMSE values at one-, two-, and four-week forecast horizons for five models. Each panel corresponds to a dataset, with errors measured in percentage points.

## Figure S2: Source selection

- Figure: [`ext_sources.pdf`](figures/supplement/ext_sources.pdf)
- Plotting code: [`plot_extension.py`](FTRF_SIR_Rebuilt/code/plot_extension.py)
- Input: `FTRF_SIR_Rebuilt/results/extension/summary.csv`

The bars compare Gaussian-interval WIS at 52 target training samples when FTRF-SIR uses the two nearest source states, all five sources, or the two farthest sources. The nearest-versus-farthest comparison holds the number of sources fixed; comparisons with all five sources also change that number.

## Figure S3: Posterior approximation and inference cost

- Figure: [`posterior_fidelity.pdf`](figures/supplement/posterior_fidelity.pdf)
- Plotting code: [`plot_results.py`](FTRF_SIR_Rebuilt/code/plot_results.py)
- Input: `FTRF_SIR_Rebuilt/results/tables/posterior_metrics.csv`

`posterior_experiment.py` calculates approximation discrepancies against a numerical-grid reference using the same data and prior. The three boxplot panels compare Laplace and Metropolis–Hastings approximations by standardized marginal Wasserstein distance, biased squared maximum mean discrepancy, and inference time, pooling both source conditions and priors. Timing excludes prediction and scoring.

## Figure S4: Marginal posterior distributions

- Figure: [`posterior_marginals.pdf`](figures/supplement/posterior_marginals.pdf)
- Plotting code: [`plot_results.py`](FTRF_SIR_Rebuilt/code/plot_results.py)
- Input: `FTRF_SIR_Rebuilt/results/tables/posterior_samples.csv`

Samples generated by `posterior_experiment.py` are restricted to matched sources, the first replicate (`rep=0`), and the transfer prior. The two panels show 45-bin density-histogram outlines for the transmission and removal parameters, comparing Grid, Laplace, and Metropolis–Hastings samples. Parameter units are per week.

## Figure S5: Observed and predicted trajectories

- Figure: [`original_forecast_trajectories.pdf`](figures/supplement/original_forecast_trajectories.pdf)
- Plotting code: [`plot_paper.py`](R_scientific_revision/plot_paper.py)
- Input: `R_scientific_revision/results/forecast_predictions.csv`

The plot selects Alabama and a one-week forecast horizon, orders records by test time, and compares observations with predictions. Rows correspond to the three outcomes, and columns to SIR, target random forest, and FTRF-SIR. The horizontal axis indexes ordered forecast origins rather than consecutive calendar weeks; the vertical axis shows percentages.
