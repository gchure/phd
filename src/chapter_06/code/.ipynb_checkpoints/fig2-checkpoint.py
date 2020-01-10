# %%
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import glob
import phd.viz
import phd.stats

colors, palettes = phd.viz.phd_style()
MAX_EXP = 100
MEAN_AREA = 4.86
MEAN_CAL = 3272

data = pd.read_csv("../../../data/csv/mscl_survival_data.csv")
# mlg910_files = glob.glob('../processing/*mlg910*/output/*.csv')
# mlg910 = pd.concat([pd.read_csv(f, comment='#') for f in mlg910_files])
# mlg910.loc[mlg910['date']==20180410, 'area'] *= 0.5
# mlg910['scaled_intensity'] = (mlg910['intensity'] - mlg910['mean_bg']) * MAX_EXP / mlg910['exposure_ms']
# mlg910['experiment'] = 'cal'
# mlg910['effective_channels'] = mlg910['scaled_intensity'] * MEAN_AREA / MEAN_CAL
# data['experiment'] = 'shock'

# data = pd.concat([data, mlg910], ignore_index=True)

#%%
# Keep only the shock data and mlg910
intensity_order = (
    data.groupby(["rbs"])["scaled_intensity"].mean().sort_values()[::-1].index
)

# Set up the figureself.
fig, ax = plt.subplots(1, 2, figsize=(6, 3))
ax[0].axis("off")
# # Plot the rescaled intensity
# _ = sns.boxplot(data['rbs'], data['scaled_intensity'] / data[data['rbs'] == 'mlg910']['scaled_intensity'].mean(), order=intensity_order,
#                 fliersize=0, linewidth=0.75, palette='Greens', ax=ax[0])
# _ = sns.stripplot(data['rbs'], data['scaled_intensity'] / data[data['rbs'] == 'mlg910']['scaled_intensity'].mean(), order=intensity_order,
#                   jitter=True, marker='.', size=2.5, alpha=0.5, color='k',
#                   ax=ax[0])
# ax[0].vlines(1, 0, 3, color=colors['pale_yellow'], lw=30, zorder=0, alpha=0.75)
ax[1].vlines(1, 0, 1000, color=colors["pale_yellow"], lw=30, zorder=0, alpha=0.75)

# Plot the channel number
channel_order = (
    data.groupby(["rbs"])["effective_channels"].mean().sort_values()[::-1].index
)

_ = sns.boxplot(
    "rbs",
    "effective_channels",
    data=data,
    order=intensity_order,
    fliersize=0,
    linewidth=0.75,
    palette="Greens",
    ax=ax[1],
)
_ = sns.stripplot(
    "rbs",
    "effective_channels",
    data=data,
    order=intensity_order,
    jitter=True,
    marker=".",
    size=2.5,
    alpha=0.5,
    color="k",
    ax=ax[1],
    zorder=1001,
)

# Add the marker for the standard candle strain
ax[0].set_ylim(0, 3)
ax[1].set_ylim(0, 1000)
labels = [i.upper() for i in intensity_order]

# Format the axes and save
# ax[0].set_ylabel('relative intensity', fontsize=8)
ax[1].set_ylabel("effective channel number", fontsize=8)
ax[0].set_xlabel("RBS modification", fontsize=8)
ax[1].set_xlabel("RBS modification", fontsize=8)

for a in ax:
    a.set_xticklabels(labels)
    a.yaxis.set_tick_params(labelsize=8)
    a.xaxis.set_tick_params(labelsize=8, rotation=45)

fig.text(0.01, 0.95, "(A)", fontsize=8)
fig.text(0.5, 0.95, "(B)", fontsize=8)
plt.tight_layout()
plt.savefig("../figs/ch6_fig2_plots.svg".format(FIG_NO), bbox_inches="tight", dpi=300)
# plt.savefig('../../figs/fig{}_plots.png'.format(FIG_NO), bbox_inches='tight', dpi=300)


# %%
