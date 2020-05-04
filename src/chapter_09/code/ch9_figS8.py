# %%
# -*- coding: utf-8 -*-
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import phd.stats
import phd.viz
import scipy.stats
colors, palette = phd.viz.phd_style()


# Load the data
data = pd.read_csv('../../data/ch9_mscl_si/complete_mcmc_traces.csv')
fig, ax = plt.subplots(1, 2, figsize=(6, 3))
phd.viz.despine(ax.ravel())
fig.text(0, 0.95, '(A)', fontsize=8)
fig.text(0.5, 0.95, '(B)', fontsize=8)
for i, a in enumerate(ax.ravel()):
    a.set_ylabel(r'$\propto$ probability')
    # a.set_yticks([])
    if i < 1:
        a.set_xlabel(r'$\beta_0$')
    else:
        a.set_xlabel(r'$\beta_1$')

color_key = ['purple', 'blue']
labels = ['slow shock', 'fast shock']
for i in range(2):
    # Evaluate the KDE of the samples.
    beta0_range = np.linspace(-20, 0, 1000)
    beta1_range = np.linspace(0.5, 3, 1000)
    beta0_kernel = scipy.stats.gaussian_kde(data['beta_0__{}'.format(i)])
    beta0_kde = beta0_kernel(beta0_range)
    beta1_kernel = scipy.stats.gaussian_kde(data['beta_1__{}'.format(i)])
    beta1_kde = beta1_kernel(beta1_range)

    _ = ax[0].plot(beta0_range, beta0_kde / beta0_kde.sum(), color=colors[color_key[i]],
                   lw=1, zorder=i+1, label=labels[i])
    _ = ax[0].fill_between(beta0_range, beta0_kde/beta0_kde.sum(),
                           color=colors['light_' + color_key[i]], zorder=i+1, alpha=0.5)
    _ = ax[1].plot(beta1_range, beta1_kde / beta1_kde.sum(), color=colors[color_key[i]],
                   lw=1, zorder=i+1, label=labels[i])
    _ = ax[1].fill_between(beta1_range, beta1_kde / beta1_kde.sum(),
                           color=colors['light_' + color_key[i]], zorder=i+1, alpha=0.5)

for a in ax:
    a.legend(fontsize=6)
# _ = ax[1].legend(fontsize=6, handlelength=1, bbox_to_anchor=(1.6, 1.0))
plt.tight_layout()
plt.savefig('../figs/ch9_figS8.pdf', bbox_inches='tight')
plt.savefig('../figs/ch9_figS8.png', bbox_inches='tight')

# %%
