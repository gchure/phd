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

data = pd.read_csv("../../data/ch6_mscl/mscl_survival_data.csv")
cal_data = pd.read_csv("../../data/ch6_mscl/mlg910_calibration_data.csv")
data['experiment'] = 'shock'
data = pd.concat([data, cal_data], ignore_index=True, sort=False)

#%%
# Keep only the shock data and mlg910
intensity_order = (
    data.groupby(["rbs"])["scaled_intensity"].mean().sort_values()[::-1].index
)

# Set up the figureself.
fig, ax = plt.subplots(1, 2, figsize=(6.5, 3))
ax[0].axis("off")
ax[1].vlines(1, 0, 1000, color='white', lw=30, zorder=0, alpha=0.75)

# Plot the channel number
channel_order = (
    data.groupby(["rbs"])["effective_channels"].mean().sort_values()[::-1].index
)

sns.boxplot()
_ = sns.boxplot(
    "rbs",
    "effective_channels",
    data=data,
    order=intensity_order,
    fliersize=0,
    linewidth=0.75,
    color=colors['black'],
    palette= sns.light_palette((210, 90, 60), n_colors=7, input="husl"), #"Greens",
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
    alpha=0.75,
    color=colors['black'],
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
    a.xaxis.set_tick_params(labelsize=8, rotation=0)

fig.text(0.01, 0.99, "(A)", fontsize=8)
fig.text(0.45, 0.99, "(B)", fontsize=8)
plt.tight_layout()
plt.savefig("../figs/ch6_fig2_plots.svg", bbox_inches="tight", dpi=300)
# plt.savefig('../../figs/fig{}_plots.png'.format(FIG_NO), bbox_inches='tight', dpi=300)


# %%
