# -*- coding: utf-8 -*-cd code/figures
#%%
import numpy as np
import pandas as pd
import scipy.stats
import matplotlib.pyplot as plt
import phd.viz
colors, palette = phd.viz.phd_style()

# Load the sbc data
sbc_data = pd.read_csv('../../data/ch7_mutants_si/IND_sbc.csv')

# ##############################################################################
# FIGURE INSTANTIATION 
# ##############################################################################
fig, ax = plt.subplots(3, 2, figsize=(6, 6))
phd.viz.despine(ax.ravel())
# Assign axes
model_axes = {'KaKi_only': 0, 'KaKi_epAI':1}
param_axes = {'sigma':0, 'ep_AI':0, 'Ka':1, 'Ki':1}

# Add labels
for i in range(2):
    for j in range(2):
        ax[j, i].set_ylim([-5, 5])
        ax[j, i].set_xlim([-0.1, 1.1])
        ax[j, i].set_xlabel('shrinkage', fontsize=8)
        ax[j, i].set_ylabel('z-score', fontsize=8)

for i in range(2):
        ax[2, i].set_xlabel('rank statistic', fontsize=8)
        ax[2, i].set_ylabel('cumulative distribution', fontsize=8)
        ax[2, i].set_xlim([0, 800])
phd.viz.titlebox(ax[0, 0], '$K_A$ and $K_I$ modified', bgcolor='white', 
                  color=colors['black'], pad=0.05, boxsize='12%')
phd.viz.titlebox(ax[0, 1], '$K_A$, $K_I$, and $\Delta\varepsilon_{AI}$ modified', bgcolor='white', 
                  color=colors['black'], pad=0.05, boxsize='12%')

# Add panel labels
fig.text(0, 0.95, '(A)', fontsize=8)
fig.text(0.5, 0.95, '(B)', fontsize=8)

# ##############################################################################
# SHRINKAGE AND Z-SCORE
# ##############################################################################
legend = {'Ka': '$K_A$', 'Ki':'$K_I$', 'ep_AI':r'$\Delta\varepsilon_{AI}$'}
for g, d in sbc_data.groupby(['model']):
    ka = d[d['param']=='Ka']
    ki = d[d['param']=='Ki']
    sig = d[d['param']=='sigma']
#     idx = ka < ki
    print(g, np.sum(ka['shrinkage'] > 0.9) / len(ka))
    if g != 'KaKi_only':
        ep = d[d['param']=='ep_AI']
        ax[0, 1].plot(ep['shrinkage'], ep['z_score'], 'o', color=colors['purple'], 
                label=legend['ep_AI'], ms=3, zorder=1000, markeredgecolor='white',
                markeredgewidth=0.25, alpha=0.75)

    ax[0, model_axes[g]].plot(sig['shrinkage'], sig['z_score'], 'o', color='k', ms=3,
                              label = '$\sigma$', markeredgecolor='white',
                              markeredgewidth=0.25, alpha=0.75)
    
    ax[1, model_axes[g]].plot(ka['shrinkage'], ka['z_score'], 'o', color=colors['red'], 
            label=legend['Ka'], ms=3, alpha=0.75, markeredgecolor='white', markeredgewidth=0.25)
    ax[1, model_axes[g]].plot(ki['shrinkage'], ki['z_score'], 'o', color=colors['blue'], 
            label=legend['Ki'], ms=3, markeredgecolor='white', markeredgewidth=0.25, alpha=0.75)

# ##############################################################################
# TRUE UNIFORM DISTRIBUTION
# ##############################################################################
n_sim = sbc_data.sim_idx.max()
L = np.arange(0, n_sim, 1)
R = sbc_data.rank_ndraws.unique()

# Envelope of cdf 99%
y = scipy.stats.randint.cdf(L, 0, R)
std = np.sqrt(y * (1 - y) / n_sim)
low_perc = np.concatenate((scipy.stats.norm.ppf(0.005, y[:-1], std[:-1]), (1.0, )))
high_perc = np.concatenate((scipy.stats.norm.ppf(0.995, y[:-1], std[:-1]), (1.0, )))
for i in range(2):
    ax[2, i].fill_between(L, low_perc, high_perc, color=colors['light_purple'], alpha=0.4,
                        label='__nolegend__')

# ##############################################################################
# RANK DISTRIBUTION
# ##############################################################################
for g, d in sbc_data.groupby(['model']):
    _ax = ax[2, model_axes[g]]
    ka = d[d['param']=='Ka']
    ki = d[d['param']=='Ki']
    sig = d[d['param']=='sigma']
    sig_x = np.sort(sig['rank'])
    ka_x = np.sort(ka['rank'])
    ki_x = np.sort(ki['rank'])
    y = np.arange(0, len(ka), 1) / len(ka)
    if g != 'KaKi_only':
        ep = d[d['param']=='ep_AI']
        print(np.median(ep['shrinkage']))
        ep_x = np.sort(ep['rank'])
        _ax.step(ep_x, y, color=colors['blue'], label=r'$\Delta\varepsilon_{AI}$', lw=0.75)
    _ax.step(ka_x, y, color=colors['orange'], label='$K_A$', lw=0.75)
    _ax.step(ki_x, y, color=colors['green'], label='$K_I$', lw=0.75)
    _ax.step(sig_x, y, color='k', label='$\sigma$', lw=0.75)

# ##############################################################################
# LEGEND DISTRIBUTION
# ##############################################################################
for a in ax.ravel():
    leg = a.legend(title='parameter', fontsize=6,
                   labelspacing=0.01, loc='upper left')
    leg.get_title().set_fontsize(6)
plt.tight_layout()
plt.savefig('../figs/ch7_figS14.pdf')

# %%
