# -*- coding: utf-8 -*- 
#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import phd.thermo
import phd.viz  
import phd.stats
import seaborn as sns
constants = phd.thermo.load_constants()
pboc = phd.viz.color_selector('pboc')
colors = phd.viz.color_selector('mut')
phd.viz.phd_style()

# Load delta F data
data = pd.read_csv('../../data/ch3_mutants/Chure2019_empirical_F_statistics.csv')
data = data[data['class']=='IND']
KaKi_epAI_samples = pd.read_csv('../../data/ch3_mutants/Chure2019_KaKi_epAI_samples.csv')
c_range = np.logspace(-3, 4, 500)
c_range[0] = 0

# ##############################################################################
# FIGURE INSTANTIATION
# ##############################################################################
fig, ax = plt.subplots(3, 4, figsize=(6,4))
phd.viz.despine(ax.ravel())
for a in ax.ravel():
    a.set_xscale('symlog', linthreshx=1E-2)
    a.set_ylim([-6, 5])
    a.set_xticks([0, 1E-2, 1E0, 1E2, 1E4])
for i in range(3):
    ax[i, 0].set_ylabel('$\Delta F$ [$k_BT$]', fontsize=6)
    ax[-1,i].set_xlabel('IPTG [µM]', fontsize=6)

# Hide unnecessary spines
for i in range(1, 3):
    for j in range(1, 4):
        ax[0, j].set_yticks([])
        ax[0, j].spines['left'].set_visible(False)
        ax[i, j].set_yticks([])
        ax[i, j].spines['left'].set_visible(False)

for i in range(2):
    for j in range(4):
        ax[i, j].set_xticks([])
        ax[i, j].spines['bottom'].set_visible(False)

# Define the axes. 
muts = {'Q294K':0, 'Q294R':1, 'F164T':2, 'Q294V':3}
ops = {'O1':0, 'O2':1, 'O3':2}
sns_colors = sns.color_palette('magma', n_colors=4)
op_colors = {'O1':sns_colors[0], 'O2':sns_colors[1], 'O3':sns_colors[2]}

# Set titles
for m, a in muts.items():
    ax[0, a].set_title(m, fontsize=6,  y=1.06)
for o, a in ops.items():
    ax[a, 0].text(-0.53, 0.5, o, fontsize=6, rotation='vertical', 
                  transform=ax[a, 0].transAxes)

# ##############################################################################
# LEGEND INFO
# ##############################################################################
for o, c in op_colors.items():
    ax[0, 0].plot([], [], 'o', color=c, label=o, markeredgewidth=0.5, markeredgecolor='white', ms=5)
leg = ax[0, 0].legend(fontsize=6, bbox_to_anchor=(2.7, 1.8), title='operator',
                      ncol=3)
leg.get_title().set_fontsize(6)

# ##############################################################################
#  DELTA F DATA
# ##############################################################################
for g, d in data.groupby(['mutant', 'operator', 'IPTGuM']):
    dF = d[d['parameter']=='delta_bohr']
    fc_mu = d[d['parameter']=='fc_mu']['median'].values[0]
    fc_sigma = d[d['parameter']=='fc_sigma']['median'].values[0]
    if (fc_mu < fc_sigma) | (1 - fc_mu < fc_sigma):
        pass
    else:
        for o, c in op_colors.items():
            if o == g[1]:
                face = 'w'
                edge = op_colors[g[1]]
                zorder=1000
            else:
                face = op_colors[g[1]] 
                edge = 'w'
                zorder=100
            _ax = ax[ops[o], muts[g[0]]]
            _ax.plot(dF['IPTGuM'], dF['median'], 'o', ms=3, color=op_colors[g[1]],
            markerfacecolor=face, markeredgewidth=0.5, label='__nolegend__',
            zorder=zorder, markeredgecolor=edge)
            _ax.vlines(dF['IPTGuM'], dF['hpd_min'], dF['hpd_max'], lw=0.75,
            color=op_colors[g[1]])

# ##############################################################################
# PREDICTIONS
# ##############################################################################
ref = phd.thermo.MWC(ka=constants['Ka'], ki=constants['Ki'], ep_ai=constants['ep_AI'],
                    effector_conc=c_range).pact()
for g, d in data.groupby(['mutant', 'operator']):
    cred_region = np.zeros((2, len(c_range)))
    _samps = KaKi_epAI_samples[(KaKi_epAI_samples['mutant']==g[0]) & 
                                (KaKi_epAI_samples['operator']==g[1])]

    for i, c in enumerate(c_range):
        pact = phd.thermo.MWC(ka=_samps['Ka'].values[::5], ki=_samps['Ki'].values[::5], ep_ai=_samps['ep_AI'].values[::5],
                            effector_conc=c).pact()
        dbohr = -np.log(pact / ref[i])
        cred_region[:, i] = phd.stats.compute_hpd(dbohr, 0.95)
    _ax=ax[ops[g[1]], muts[g[0]]]
    _ax.fill_between(c_range, cred_region[0, :], cred_region[1, :], 
                    color=op_colors[g[1]], alpha=0.5)

plt.subplots_adjust(hspace=0.05, wspace=0.05)
plt.savefig('../figs/ch7_figS19.pdf', bbox_inches='tight')
plt.savefig('../figs/ch7_figS19.png', bbox_inches='tight')
