"""Draw final manuscript figures from saved results without fitting models."""
from pathlib import Path
import shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'R_scientific_revision/results'
FIG = ROOT / 'latex/English/figures'
plt.rcParams.update({'font.size': 11, 'font.family': 'DejaVu Sans',
                     'pdf.fonttype': 42, 'axes.spines.top': False,
                     'axes.spines.right': False, 'axes.titlesize': 11})
NAMES = {'RF_tar': 'Target RF', 'FTRF': 'Transfer RF', 'SR': 'Ridge stacking',
         'GMA': 'Greedy averaging', 'RF_agg': 'Pooled RF',
         'TL_RF': 'RF residual transfer', 'TL_KRR': 'Kernel ridge transfer',
         'SIR': 'SIR', 'FTRF-SIR': 'FTRF-SIR', 'Persistence': 'Persistence',
         'Persistence_residual': 'Persistence residual',
         'SIR_residual': 'SIR residual', 'Half_SIR_RF': 'Half SIR / half RF',
         'All_SIR_mean': 'All-region SIR average',
         'Auxiliary_SIR_mean': 'Source SIR average'}
TASKS = {'CLS_B': 'Clinical influenza B', 'CPT_B': 'Combined influenza B',
         'CPT_B_rate': 'Combined influenza B', 'ILI': 'Influenza-like illness',
         'PHL_A': 'Public health influenza A', 'PHL_B': 'Public health influenza B'}

def save(fig, name):
    fig.savefig(FIG / (name + '.pdf'), bbox_inches='tight')
    plt.close(fig)
    shutil.copy2(FIG / (name + '.pdf'), ROOT / 'latex/Chinese/figures' / (name + '.pdf'))
    research_figures = ROOT / 'R_scientific_revision/figures'
    research_figures.mkdir(exist_ok=True)
    shutil.copy2(FIG / (name + '.pdf'), research_figures / (name + '.pdf'))

def main():
    blocks = pd.read_csv(DATA / 'conditional_block_errors.csv')
    methods = list(pd.read_csv(DATA / 'conditional_predictions.csv', usecols=['model']).model.unique())
    fig, axs = plt.subplots(2, 3, figsize=(12, 8.5))
    for ax, task in zip(axs.flat, sorted(blocks.task.unique())):
        v = blocks[blocks.task == task]
        labels = [NAMES[m].replace(' ', '\n', 1) for m in methods]
        ax.boxplot([v[v.model == m].SE for m in methods], tick_labels=labels, showfliers=False)
        ax.tick_params(axis='x', rotation=60, labelsize=10)
        ax.set_title(TASKS[task])
        ax.set_ylabel('MSE (squared percentage points)' if task in ['ILI', 'CLS_B'] else 'MSE (squared counts)')
    axs.flat[-1].axis('off')
    fig.tight_layout()
    save(fig, 'original_five_tasks')
    for group, column in [('forecast', 'task'), ('simulation', 'condition')]:
        data = pd.read_csv(DATA / (group + '_summary.csv'))
        values = sorted(data[column].unique())
        fig, axs = plt.subplots(1, len(values), figsize=(11, 4.2), squeeze=False)
        for ax, value in zip(axs.flat, values):
            v = data[data[column] == value].sort_values('RMSE')
            ax.barh([NAMES[m] for m in v.model], v.RMSE, color='#3876a8')
            ax.set_title(TASKS[value] if group == 'forecast' else value.capitalize() + ' sources')
            ax.set_xlabel('RMSE (percentage points)')
            ax.invert_yaxis()
        fig.tight_layout()
        save(fig, group + '_comparison')
    data = pd.read_csv(DATA / 'forecast_predictions.csv')
    data = data[(data.target == 'Alabama') & (data.h == 1)]
    fig, axs = plt.subplots(3, 3, figsize=(10, 8.3))
    for i, task in enumerate(sorted(data.task.unique())):
        for j, method in enumerate(['SIR', 'RF_tar', 'FTRF-SIR']):
            v = data[(data.task == task) & (data.model == method)].sort_values('test_time')
            ax = axs[i, j]
            ax.plot(np.arange(len(v)), v.observed, 'o-', color='#3876a8', label='Observed', ms=3)
            ax.plot(np.arange(len(v)), v.predicted, 's--', color='#b4872c', label='Predicted', ms=3)
            ax.set_title(NAMES[method] if i == 0 else '')
            if j == 0:
                ax.set_ylabel(TASKS[task] + '\n(%)')
            if i == 2:
                ax.set_xlabel('Ordered forecast origin')
    fig.legend(*axs[0, 0].get_legend_handles_labels(), loc='upper center',
               ncol=2, frameon=False, bbox_to_anchor=(.5, 1.025))
    fig.tight_layout()
    save(fig, 'original_forecast_trajectories')

if __name__ == '__main__':
    main()
