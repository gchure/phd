# %%
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import sys
import glob
import phd.stats
import phd.viz
colors, palette = phd.viz.phd_style()

# -------------
data = pd.read_csv('../../data/ch5_mscl/mscl_survival_data.csv')
grouped = data.groupby('survival')

# Load the MLG910 data and compute the area
mlg910 = pd.read_csv('../../data/ch5_mscl/mlg910_calibration_data.csv')

#%%
# Set up the axis
fig, ax = plt.subplots(2, 2, figsize=((6, 5)))
phd.viz.despine(ax.ravel())
ax = ax.ravel()
for i, a in enumerate(ax):
    a.xaxis.set_tick_params(labelsize=8)
    a.yaxis.set_tick_params(labelsize=8)
    a.set_xlabel('projected cell area [µm$^2$]', fontsize=8)
    if i >= 2: 
        a.set_yscale('log')
        a.set_ylim([5, 1E4])
ax[0].axis('off')

# Labels and titles
ax[1].set_ylabel('cumulative distribution', fontsize=8)
phd.viz.titlebox(ax[2], 'total channel count', size=8, bgcolor='white', 
                color=colors['black'], pad=0.05, boxsize='12%')
phd.viz.titlebox(ax[3], 'effective channel count', size=8, bgcolor='white', 
                color=colors['black'], pad=0.05, boxsize='12%')
ax[2].set_ylabel('total channels per cell', fontsize=8)
ax[3].set_ylabel('effective channels per cell', fontsize=8)
fig.text(0, 1, '(A)', fontsize=8)
fig.text(0.5, 1, '(B)', fontsize=8)
fig.text(0, 0.5, '(C)', fontsize=8)
fig.text(0.5, 0.5, '(D)', fontsize=8)

# Plot the cumulative distributions for the standard candle and experimental
# cells
mlg910_x, mlg910_y = np.sort(mlg910['area'].values), np.arange(0, len(mlg910)) / len(mlg910)
data_x, data_y = np.sort(data['area'].values), np.arange(0, len(data)) / len(data)
_ = ax[1].step(mlg910_x, mlg910_y, color=colors['orange'], lw=1.5, label='MLG910')
_ = ax[1].step(data_x, data_y, color=colors['purple'], lw=1.5, label='SD mutants')
_ = ax[1].legend(loc='lower right', fontsize=8)

color_dict = {True: colors['blue'], False: colors['green']}
label_dict = {True: 'survival', False: 'death'}
for g, d in grouped:
    _ = ax[2].plot(d['area'], d['scaled_intensity'] * d['area'] / d['calibration_factor'], '.',
                   color=color_dict[g], ms=5, alpha=0.5, label=label_dict[g],
                   markeredgecolor='white', markeredgewidth=0.5)
    _ = ax[3].plot(d['area'], d['effective_channels'], '.', color=color_dict[g],
                   ms=5, alpha=0.5, label=label_dict[g], markeredgecolor='white',
                   markeredgewidth=0.5)
_ = ax[2].legend(fontsize=6)
plt.tight_layout()
plt.savefig('../figs/figS5_plots.svg', bbox_inches='tight')

# %%
