# -*- coding: utf-8 -*-
#%%
import numpy as np
import pandas as pd
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import seaborn as sns
import phd.viz
colors, palette = phd.viz.phd_style()

# Load the posterior check data. 
data = pd.read_csv('../../data/ch7_mutants_si/compiled_data.csv')
data = data[(data['mutant']=='Q294K') & (data['operator']=='O2')]
_kaki_data = pd.read_csv('../../data/ch7_mutants_si/IND_posterior_predictive_samples.csv')
kaki_data = _kaki_data[_kaki_data['model']=='KaKi_only']
kaki_epAI_data = _kaki_data[_kaki_data['model']=='KaKi_epAI']
kaki_samples = pd.read_csv('../../data/ch3_mutants/Chure2019_KaKi_only_samples.csv')
kaki_epAI_samples = pd.read_csv('../../data/ch3_mutants/Chure2019_KaKi_epAI_samples.csv')
kaki_samples = kaki_samples[(kaki_samples['mutant']=='Q294K') &
                             (kaki_samples['operator']=='O2')]
kaki_epAI_samples = kaki_epAI_samples[(kaki_epAI_samples['mutant']=='Q294K') &
                             (kaki_epAI_samples['operator']=='O2')]

cmap = sns.color_palette('magma', n_colors=len(kaki_samples) + 2)
kaki_samples.sort_values('lp__', inplace=True)
kaki_samples['color'] = cmap[:-2]

cmap = sns.color_palette('magma', n_colors=len(kaki_epAI_samples) + 2)
kaki_epAI_samples.sort_values('lp__', inplace=True)
kaki_epAI_samples['color'] = cmap[:-2]

# ##############################################################################
# FIGURE 1 -- KaKi ONLY INSTANTIATION
# ############################################################################## 
fig = plt.figure(figsize=(6, 3))
gs = gridspec.GridSpec(6, 11)
ax1 = fig.add_subplot(gs[0:2, 0:2])
ax2 = fig.add_subplot(gs[2:4, 0:2])
ax3 = fig.add_subplot(gs[2:4, 2:4])
ax4 = fig.add_subplot(gs[4:6, 0:2])
ax5 = fig.add_subplot(gs[4:6, 2:4])
ax6 = fig.add_subplot(gs[4:6, 4:6])
ax7 = fig.add_subplot(gs [:, 7:])
ax = [ax1, ax2, ax3, ax4, ax5, ax6, ax7]
phd.viz.despine(ax)

# Apply special formatting where needed
for a in ax:
    a.xaxis.set_tick_params(labelsize=6)
    a.yaxis.set_tick_params(labelsize=6)
for a in [ax1, ax2, ax3]:
    a.set_xticklabels([])
for a in [ax1, ax3, ax5, ax6]:
    a.set_yticks([])

# Adjust scaling
ax7.set_xscale('symlog', linthreshx=1E-2)
ax7.set_ylim([-0.4, 1.1])
# Add labels
ax4.set_xlabel('$K_A$ [µM]')
ax2.set_ylabel('$\sigma$ ')
ax4.set_ylabel('$K_I$ [µM]')
ax5.set_xlabel('$\sigma$')
ax6.set_xlabel('$K_I$ [µM]')
ax7.set_xlabel('IPTG [µM]')
ax7.set_ylabel('fold-change')

# Add panel labels
fig.text(0.02, 0.95, '(A)')
fig.text(0.5, 0.95, '(B)')
# ##############################################################################
# JOINT DISTRIBUTIONS
# ##############################################################################
ax2.scatter(kaki_samples['Ka'].values, kaki_samples['sigma'].values, marker='.',
            c=kaki_samples['color'].values,  s=0.5, rasterized=True)
ax4.scatter(kaki_samples['Ka'].values, kaki_samples['Ki'].values, marker='.',
            c=kaki_samples['color'].values, s=0.5, rasterized=True)
ax5.scatter(kaki_samples['sigma'].values, kaki_samples['Ki'].values, marker='.',
            c=kaki_samples['color'].values, s=0.5, rasterized=True)
            
# ##############################################################################
# MARGINAL DISTRIBUTIONS
# ##############################################################################
ax1.hist(kaki_samples['Ka'], bins=30, histtype='stepfilled', color=colors['purple'], 
         alpha=0.5)
ax3.hist(kaki_samples['sigma'], bins=30, histtype='stepfilled', color=colors['purple'], 
         alpha=0.5)
ax6.hist(kaki_samples['Ki'], bins=30, histtype='stepfilled', color=colors['purple'],
         alpha=0.5)

# ##############################################################################
# PERCENTILES
# ##############################################################################
percs = [99, 95, 80, 50, 20, 10, 5]
perc_cmap = sns.color_palette('magma', n_colors=(len(percs) + 2))
zorder = [11, 12, 13, 14, 15, 16, 17]
z = {p:z for p, z in zip(percs, zorder)}
c = {p:c for p, c in zip(percs, perc_cmap)}
grouped = kaki_data.groupby(['IPTGuM'])
df = pd.DataFrame([])
for g, d in grouped:
    for p in percs:     
        remainder = 100 - p
        low = remainder / 2
        upper = p + remainder / 2 
        _percs = np.percentile(d['fc_mu'], [low, upper])
        df = df.append({'percentile': p,
                         'IPTGuM': g,
                         'fc_low':_percs[0],
                         'fc_high': _percs[1]},
                      ignore_index=True)

for g, d in df.groupby('percentile'):
    ax7.fill_between(d['IPTGuM'], d['fc_low'], d['fc_high'], color=c[g], 
                    zorder=z[g], label=int(g), alpha=0.5)
  
# ##############################################################################
# EXPERIMENTAL MEASUREMENTS
# ##############################################################################
ax7.plot(data['IPTGuM'], data['fold_change'], 'o', color=colors['blue'], 
         markeredgecolor='w', markeredgewidth=0.5, ms=4, label='data', zorder=1000)
leg = ax7.legend(title='percentile', loc='upper left', fontsize=6, ncol=2)
leg.get_title().set_fontsize(6)
plt.subplots_adjust(hspace=0.2, wspace=0.2)
plt.savefig('../figs/ch7_figS15.pdf',  bbox_inches='tight')
plt.savefig('../figs/ch7_figS15.png',  bbox_inches='tight')

#%%

# ##############################################################################
# FIGURE 2 -- KaKi + ∆epAI INSTANTIATION
# ##############################################################################
fig = plt.figure(figsize=(6, 3))
gs = gridspec.GridSpec(4, 9)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[1, 0])
ax3 = fig.add_subplot(gs[1, 1])
ax4 = fig.add_subplot(gs[2, 0])
ax5 = fig.add_subplot(gs[2, 1])
ax6 = fig.add_subplot(gs[2, 2])
ax7 = fig.add_subplot(gs[3, 0])
ax8 = fig.add_subplot(gs[3, 1])
ax9 = fig.add_subplot(gs[3, 2])
ax10 = fig.add_subplot(gs[3, 3])
ax11 = fig.add_subplot(gs[:, 5:])
ax = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11]
phd.viz.despine(ax)
# Adjust scaling
ax11.set_ylim([0, 1.1])
ax11.set_xscale('symlog', linthreshx=1E-2)

# Apply special formatting where needed
for a in ax:
    a.xaxis.set_tick_params(labelsize=5)
    a.yaxis.set_tick_params(labelsize=5)
for a in [ax1, ax2, ax3, ax4, ax5, ax6]:
    a.set_xticklabels([])
for a in [ax1, ax3, ax5, ax6, ax8, ax9, ax10]:
    a.set_yticklabels([])
for a in [ax1, ax3, ax6, ax10]:
    a.set_yticks([])


# Add labels
ax2.set_ylabel('$\sigma$')
ax4.set_ylabel('$K_A$ [µM]')
ax4.set_ylabel('$K_I$ [µM]')
ax7.set_ylabel(r'$\Delta\varepsilon_{AI}$ [$k_BT$]')
ax7.set_xlabel('$K_A$ [µM]')
ax9.set_xlabel('$K_I$ [µM]')
ax8.set_xlabel(r'$\sigma$')
ax10.set_xlabel(r'$\Delta\varepsilon_{AI}$ [$k_BT$]')
ax11.set_xlabel('IPTG [µM]')
ax11.set_ylabel('fold-change')

# Add panel labels
fig.text(0.02, 0.95, '(A)', fontsize=8)
fig.text(0.5, 0.95, '(B)', fontsize=8)
# ##############################################################################
# JOINT DISTRIBUTIONS
# ##############################################################################
ax2.scatter(kaki_epAI_samples['Ka'].values, kaki_epAI_samples['sigma'].values, marker='.',
            c=kaki_epAI_samples['color'].values,  s=0.5, rasterized=True)
ax4.scatter(kaki_epAI_samples['Ka'].values, kaki_epAI_samples['Ki'].values, marker='.',
            c=kaki_epAI_samples['color'].values, s=0.5, rasterized=True)
ax5.scatter(kaki_epAI_samples['sigma'].values, kaki_epAI_samples['Ki'].values, marker='.',
            c=kaki_epAI_samples['color'].values, s=0.5, rasterized=True)
ax7.scatter(kaki_epAI_samples['Ka'].values, kaki_epAI_samples['ep_AI'].values, marker='.',
            c=kaki_epAI_samples['color'].values, s=0.5, rasterized=True)
ax8.scatter(kaki_epAI_samples['sigma'].values, kaki_epAI_samples['ep_AI'].values, marker='.',
            c=kaki_epAI_samples['color'].values, s=0.5, rasterized=True)
ax9.scatter(kaki_epAI_samples['Ki'].values, kaki_epAI_samples['ep_AI'].values, marker='.',
            c=kaki_epAI_samples['color'].values, s=0.5, rasterized=True)
         
# ##############################################################################
# MARGINAL DISTRIBUTIONS
# ##############################################################################
_ = ax1.hist(np.log10(kaki_epAI_samples['Ka']), bins=30, histtype='stepfilled', color=colors['purple'], 
         alpha=0.5)
_ = ax3.hist(kaki_epAI_samples['sigma'], bins=30, histtype='stepfilled', color=colors['purple'], 
         alpha=0.5)
_ = ax6.hist(kaki_epAI_samples['Ki'], bins=30, histtype='stepfilled', color=colors['purple'],
         alpha=0.5)
_ = ax10.hist(kaki_epAI_samples['ep_AI'], bins=30, histtype='stepfilled', color=colors['purple'],
         alpha=0.5)

# ##############################################################################
# PERCENTILES
# ##############################################################################
percs = [99, 95, 80, 50, 20, 10, 5]
perc_cmap = sns.color_palette('magma', n_colors=(len(percs) + 2))
zorder = [11, 12, 13, 14, 15, 16, 17]
z = {p:z for p, z in zip(percs, zorder)}
c = {p:c for p, c in zip(percs, perc_cmap)}
grouped = kaki_epAI_data.groupby(['IPTGuM'])
df = pd.DataFrame([])
for g, d in grouped:
    for p in percs:     
        remainder = 100 - p
        low = remainder / 2
        upper = p + remainder / 2 
        _percs = np.percentile(d['fc_mu'], [low, upper])
        df = df.append({'percentile': p,
                         'IPTGuM': g,
                         'fc_low':_percs[0],
                         'fc_high': _percs[1]},
                      ignore_index=True)

for g, d in df.groupby('percentile'):
    ax11.fill_between(d['IPTGuM'], d['fc_low'], d['fc_high'], color=c[g], 
                    zorder=z[g], label=int(g), alpha=0.5)
  
# ##############################################################################
# EXPERIMENTAL MEAUREMENTS
# ##############################################################################
ax11.plot(data['IPTGuM'], data['fold_change'], 'o', color=colors['blue'], 
         label='data', zorder=1000, ms=4, markeredgewidth=0.5, markeredgecolor='white')
leg = ax11.legend(title='percentile', loc='upper left', fontsize=6, ncol=2)
leg.get_title().set_fontsize(6)
plt.subplots_adjust(hspace=0.2, wspace=0.2)
plt.savefig('../figs/ch7_figS16.pdf', bbox_inches='tight')
plt.savefig('../figs/ch7_figS16.png', bbox_inches='tight')


# %%
