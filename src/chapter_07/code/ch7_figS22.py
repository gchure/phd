
#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import phd.viz
import phd.thermo
import phd.stats
constants = phd.thermo.load_constants()
colors, palette = phd.viz.phd_style()
palette = sns.color_palette('magma', n_colors=4)
op_colors = {'O1':palette[0], 'O2':palette[1], 'O3':palette[2]} 

# Load the data and sampling information
data = pd.read_csv('../../data/ch3_mutants/Chure2019_summarized_data.csv', comment='#')
IND_data = data[data['class']=='IND']
samples = pd.read_csv('../../data/ch7_mutants_si/Chure2019_global_KaKi_epAI_samples.csv', comment="#")
empirical_F = pd.read_csv('../../data/ch3_mutants/Chure2019_empirical_F_statistics.csv', comment='#')
empirical_F = empirical_F[empirical_F['class']=='IND']

c_range = np.logspace(-2, 4, 200)
c_range[0] = 0

# ##############################################################################
# FIGURE INSTANTIATION
# ##############################################################################
fig, ax = plt.subplots(2, 4, figsize=(6, 4))
phd.viz.despine(ax.ravel(0))
for a in ax.ravel():
    a.set_xscale('symlog', linthreshx=1E-2)
    a.set_xlabel('IPTG [µM]', fontsize=8)
    a.set_xlim([-.001, 1E4])

for i in range(4):
    ax[0, i].set_ylim([-0.01, 1.15])
    ax[1, i].set_ylim([-8, 8])
for i in range(3):
    ax[0, i+1].set_yticklabels([])
    ax[1, i+1].set_yticklabels([])
ax[0, 0].set_ylabel('fold-change', fontsize=8)
ax[1, 0].set_ylabel(r'$\Delta F$ [$k_BT$]', fontsize=8)

# Define the axes. 
axes = {'Q294K':0, 'Q294R':1, 'F164T':2, 'Q294V':3}
for m, a in axes.items():
        if m == 'Q294K':
                label = 'Q291K'
        elif m == 'Q294R':
                label = 'Q291R'
        elif m == 'Q294V':
                label = 'Q291V'
        elif m == 'F164T':
                label = 'F161T'
        else:
                label = m
        phd.viz.titlebox(ax[0, a], label, bgcolor='white', color=colors['black'],
                        pad=0.05, boxsize='12%')

# ##############################################################################
# INDUCTION PROFILES
# ##############################################################################
for g, d in samples.groupby(['mutant']):
    for o, c in op_colors.items():
        _ax = ax[0, axes[g]]
        cred_region = np.zeros((2, len(c_range)))
        for i, _c in enumerate(c_range):
            arch = phd.thermo.SimpleRepression(R=constants['RBS1027'], 
                    ep_r=constants[o], ka=d['Ka'], ki=d['Ki'], ep_ai=d['ep_AI'],
                    effector_conc=_c).fold_change()
            cred_region[:, i] = phd.stats.compute_hpd(arch, 0.95)
        _ax.fill_between(c_range, cred_region[0, :], cred_region[1, :], color=c, 
            alpha=0.5) 

# ##############################################################################
# DELTA F
# ##############################################################################
wt_pact = phd.thermo.MWC(ka=constants['Ka'], ki=constants['Ki'], effector_conc=c_range,
                        ep_ai=constants['ep_AI']).pact()
for g, d in samples.groupby(['mutant']):
    _ax = ax[1, axes[g]]
    cred_region = np.zeros((2, len(c_range)))
    for i, _c in enumerate(c_range):
        arch = phd.thermo.MWC(ka=d['Ka'], ki=d['Ki'], ep_ai=d['ep_AI'],
                    effector_conc=_c).pact()
        dBohr = -np.log(arch / wt_pact[i])
        cred_region[:, i] = phd.stats.compute_hpd(dBohr, 0.95)
    _ax.fill_between(c_range, cred_region[0, :], cred_region[1, :], color='slategray', 
            alpha=0.5) 

# ##############################################################################
# DATA
# ##############################################################################
for g, d in IND_data.groupby(['operator', 'mutant']):
    _ax = ax[0, axes[g[1]]]
    _ax.errorbar(d['IPTGuM'], d['mean'], d['sem'], lw=0.75, capsize=1, 
                color=op_colors[g[0]], label=g[0], linestyle='None',
                fmt='o', ms=4.5, markeredgewidth=0.5, markeredgecolor='white')

# ##############################################################################
# DELTA F DATA
# ##############################################################################
for g, d in empirical_F.groupby(['mutant', 'operator', 'IPTGuM']):
    mu = d[d['parameter']=='fc_mu']['median'].values[0]
    sig = d[d['parameter']=='fc_sigma']['median'].values[0]
    dF = d[d['parameter']=='delta_bohr']
    _ax = ax[1, axes[g[0]]]
    if (mu > sig) & (1 - mu > sig):
        _ax.plot(dF['IPTGuM'], dF['median'], 'o', color=op_colors[g[1]], ms=4.5, 
                 linestyle='None', markeredgecolor='w', markeredgewidth=0.5)
        _ax.vlines(dF['IPTGuM'], dF['hpd_min'], dF['hpd_max'], lw=0.75, 
                   color=op_colors[g[1]])

# ##############################################################################
# LEGEND AND SAVING
# ##############################################################################
leg = ax[0, 0].legend(title='operator', fontsize=6)
leg.get_title().set_fontsize(6)
plt.tight_layout()
plt.savefig('../figs/ch7_figS22.pdf', bbox_inches='tight')
plt.savefig('../figs/ch7_figS22.png', bbox_inches='tight')





# %%
