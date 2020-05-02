# -*- coding: utf-8 -*-
# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import phd.viz
import phd.stats

colors, palette = phd.viz.phd_style()

# Load in the data and isolate to onlly the shock experiments.
shock_data = pd.read_csv("../../data/ch5_mscl/mscl_survival_data.csv")

# Define the colors
color_dict = {True: colors["purple"], False: colors["orange"]}
label_dict = {True: "survival", False: "death"}
zorder_dict = {True: 1, False: 2}
alpha_dict = {True: 0.9, False: 0.65}

# Group by survival and set the colors.
grouped = shock_data.groupby(["survival"])

#%% Generate the same plots but separated by shock rate.
grouped = shock_data.groupby(["shock_class", "survival"])

# Set up the figure canvas.
fig, ax = plt.subplots(1, 2, figsize=(6, 3))
phd.viz.despine(ax)
fig.text(0.01, 0.99, "(A)", fontsize=8)
fig.text(0.5, 0.99, "(B)", fontsize=8)
axes = {"slow": ax[0], "fast": ax[1]}
for a in ax:
    a.tick_params(labelsize=8)
for g, d in grouped:
    # Compute the ecdf.
    y = np.arange(0, len(d), 1) / len(d)
    x_median = np.sort(d["effective_channels"])
    x_min = np.sort(d["minimum_channels"])
    x_max = np.sort(d["maximum_channels"])

    # Plot the ecdf.
    axes[g[0]].plot(
        x_median, y, ".", ms=4, alpha=0.75, color=color_dict[g[1]], label="__nolegend__")
    axes[g[0]].fill_betweenx(
        y, x_min, x_max, color=color_dict[g[1]], alpha=0.3, label="__nolegend__"
    )
    axes[g[0]].plot(x_max, y, "-", lw=1, color=color_dict[g[1]], label=label_dict[g[1]])
    axes[g[0]].plot(x_min, y, "-", lw=1, color=color_dict[g[1]], label="__nolegend__")
    min_surv = np.min(x_median)
    print(g, min_surv)

    if g[1] == 1:
        axes[g[0]].fill_betweenx(
            np.linspace(0, 1, 300), 0, min_surv, color=colors['black'], alpha=0.25
        )

ax[0].legend(fontsize=8)
for a in ax:
    a.set_xlim(
        [shock_data["effective_channels"].min(), shock_data["effective_channels"].max()]
    )
    a.set_ylim([0, 1])
    a.set_xlabel("effective channels per cell", fontsize=8)
    a.set_ylabel("cumulative distribution", fontsize=8)
phd.viz.titlebox(ax[0], "slow shock ($<$ 1.0 Hz)", color=colors['black'], bgcolor=colors['grey'],
                fontsize=8, boxsize="10%")
phd.viz.titlebox(ax[1], "fast shock ($\geq$ 1.0 Hz)", color=colors['black'], bgcolor=colors['grey'],
                fontsize=8, boxsize="10%")

plt.tight_layout()
plt.savefig("../figs/ch5_fig4.pdf", bbox_inches="tight", dpi=300)
plt.savefig("../figs/ch5_fig4.png", bbox_inches="tight", dpi=300)


# %%
