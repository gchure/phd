#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats
import phd.stats
import phd.viz
colors, palette = phd.viz.phd_style()

# Load the samples
sbc_samples = pd.read_csv('../../data/ch7_mutants_si/epRA_sbc.csv')


fig, ax = plt.subplots(2, 2, figsize=(6, 4))
phd.viz.despine(ax.ravel())

# Add labels
for i in range(2):
    ax[i, 0].set_xlabel('DNA binding energy [$k_BT$]')
    ax[i, 1].set_xlabel('$\sigma$')
    ax[0, i].set_ylabel('$\propto$ probability')
    ax[1, i].set_ylabel('cumulative distribution')

# Add title
phd.viz.titlebox(ax[0,0], r'DNA binding energy $\Delta\varepsilon_{RA}$', size=8, 
                color=colors['black'], bgcolor='white', pad=0.05, boxsize='12%')
phd.viz.titlebox(ax[0,1], r'standard deviation $\sigma$', size=8, 
                color=colors['black'], bgcolor='white', pad=0.05, boxsize='12%')

axes = {'ep_RA':0, 'sigma':1}
fig.text(0, 0.95, '(A)', fontsize=8)
fig.text(0.5, 0.95, '(B)', fontsize=8)


#  Plot the ground truth distributions
for g, d in sbc_samples.groupby(['param']):
    hist_ax = ax[0, axes[g]]
    ecdf_ax = ax[1, axes[g]]

    # Histogram
    hist, bins = np.histogram(d['ground_truth'], bins=15, density=True)
    hist_ax.step(bins[:-1], hist, color=colors['purple'], lw=0.75, label='ground truth')
    hist_ax.fill_between(bins[:-1], hist, step='pre', color=colors['light_purple'], 
                         alpha=0.4)
    
    # ECDF
    x, y = np.sort(d['ground_truth']), np.arange(0, len(d), 1) / len(d)
    ecdf_ax.step(x, y, color=colors['purple'], lw=0.75, label='ground truth')

# Plot the SBC distributions
for g, d in sbc_samples.groupby(['param']):
    hist_ax = ax[0, axes[g]]
    ecdf_ax = ax[1, axes[g]]

    # Histogram
    hist, bins = np.histogram(d['post_mean'], bins=15, density=True)
    hist_ax.step(bins[:-1], hist, color=colors['orange'], lw=1, label='inferred')
    hist_ax.fill_between(bins[:-1], hist, step='pre', color=colors['light_orange'], alpha=0.4)
    
    # ECDF
    x, y = np.sort(d['post_mean']), np.arange(0, len(d), 1) / len(d)
    ecdf_ax.step(x, y, color=colors['orange'], lw=0.75, label='inferred')

# ##############################################################################
# LEGENDS
# ##############################################################################
ax[0, 0].legend(fontsize=8, handlelength=0.5)
plt.tight_layout()
plt.savefig('../figs/ch7_figS3.pdf', bbox_inches='tight')
plt.savefig('../figs/ch7_figS3.png', bbox_inches='tight')


#%%
# ##############################################################################
# FIGURE S4: SENSITIVITY
# ##############################################################################
fig, ax = plt.subplots(1, 2, figsize=(4, 2.5), sharex=True, sharey=True)
phd.viz.despine(ax.ravel())

# Add labels
for i in range(2):
    ax[i].set_xlabel('posterior shrinkage $s$')
ax[0].set_ylabel('posterior $z$-score')

# Add titles
phd.viz.titlebox(ax[0], r'DNA binding energy $\Delta\varepsilon_{RA}$',
                color=colors['black'], bgcolor='white', size=8, pad=0.05,
                boxsize='12%')
phd.viz.titlebox(ax[1], r'standard deviation $\sigma$',
                color=colors['black'], bgcolor='white', size=8, pad=0.05,
                boxsize='12%')

# Add panel labels
fig.text(0.02, 0.9, '(A)', fontsize=8)
fig.text(0.53, 0.9, '(B)', fontsize=8)

# Adjust scaling
for i in range(2):
    ax[i].set_xlim([0, 1.05])
    ax[i].set_ylim([-5, 5])

# Assign axes
axes = {'ep_RA':0, 'sigma':1}
# ##############################################################################
# SBC DATA
# ##############################################################################
for g, d in sbc_samples.groupby(['param']):
    _ax = ax[axes[g]]
    _ax.plot(d['shrinkage'], d['z_score'], 'o', color=colors['orange'], ms=2,
    markeredgecolor='white', markeredgewidth=0.25, alpha=0.5)

plt.tight_layout()
plt.savefig('../figs/ch7_figS4.pdf', bbox_inches='tight')
plt.savefig('../figs/ch7_figS4.png', bbox_inches='tight')
#%%
# ##############################################################################
# FIGURE 3: RANK DISTRIBUTION
# ##############################################################################
fig, ax = plt.subplots(2, 2, figsize=(6, 4))
phd.viz.despine(ax.ravel())
for a in ax.ravel():
    a.set_xlabel('rank statistic')

# Formatting and scale
for i in range(2):
    ax[0, i].set_xlim([0, 800])
    ax[1, i].set_ylabel('cumulative distribution')
    ax[0, i].set_ylabel('counts')

# Add labels
phd.viz.titlebox(ax[0, 0], r'DNA binding energy $\Delta\varepsilon_{RA}$',
                 color=colors['black'], bgcolor='white', pad=0.05, 
                 boxsize='12%')
phd.viz.titlebox(ax[0, 1], r'standard deviation $\sigma$',
                 color=colors['black'], bgcolor='white', pad=0.05, 
                 boxsize='12%')

fig.text(0.0, 0.95, '(A)', fontsize=8)
fig.text(0.5, 0.95, '(B)', fontsize=8)
# Define the axes
axes = {'ep_RA': 0, 'sigma':1}
# ##############################################################################
# TRUE UNIFORM DISTRIBUTION
# ##############################################################################
n_sim = sbc_samples.sim_idx.max()
bins = 20
L = np.arange(0, n_sim, 1)
R = 800

# Bounds for histogram 99% 
low = scipy.stats.binom.ppf(0.005, R, 2 * bins/R)
high = scipy.stats.binom.ppf(0.995, R, 2* bins/R)

# Envelope of cdf 99%
y = scipy.stats.randint.cdf(L, 0, R)
std = np.sqrt(y * (1 - y) / n_sim)
low_perc = np.concatenate((scipy.stats.norm.ppf(0.005, y[:-1], std[:-1]), (1.0, )))
high_perc = np.concatenate((scipy.stats.norm.ppf(0.995, y[:-1], std[:-1]), (1.0, )))

# ##############################################################################
#  DATA DISTRIBUTION
# ##############################################################################
for g, d in sbc_samples.groupby('param'):
    hist_ax = ax[0, axes[g]]
    ecdf_ax = ax[1, axes[g]]

    # Bin the histogram
    _ = hist_ax.hist(d['rank'], bins=bins, color=colors['light_orange'], 
                edgecolor=colors['orange'])

    # Percentile bounds
    _ = hist_ax.fill_between(L, low, high, color=colors['light_purple'], alpha=0.4, zorder=100)

    # ECDF
    x, y = np.sort(d['rank']), np.arange(0, len(d), 1) / len(d)
    ecdf_ax.step(x, y, color=colors['orange'])

    # Percentile_bounds
    ecdf_ax.fill_between(L, low_perc, high_perc, color=colors['light_purple'], alpha=0.4)

plt.tight_layout()
plt.savefig('../figs/ch7_figS5.pdf', bbox_inches='tight')
plt.savefig('../figs/ch7_figS5.png', bbox_inches='tight')

# %%
