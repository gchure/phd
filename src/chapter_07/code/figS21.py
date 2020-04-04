
#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import phd.viz
import phd.thermo
constants = phd.thermo.load_constants()
colors, palette = phd.viz.phd_style()

# Define styling and axis mapping
rep_colors = {60: colors['brown'],  
              124: colors['green'], 
              260: colors['orange'], 
              1220: colors['purple']}
             
# Load the necessary data and statistics
data = pd.read_csv('../../data/ch3_mutants/Chure2019_summarized_data.csv', comment='#')
DNA_data = data[data['class']=='DNA']
epRA_stats =  pd.read_csv('../../data/ch7_mutants_si/Chure2019_global_DNA_binding_energy_summary.csv')
empirical_F = pd.read_csv('../../data/ch3_mutants/Chure2019_empirical_F_statistics.csv')
empirical_F = empirical_F[empirical_F['class']=='DNA']

# Define plotting constants
c_range = np.logspace(-2, 4, 200)
c_range[0] = 0

# ##############################################################################
# FIGURE INSTANTIATION
# ##############################################################################
fig, ax = plt.subplots(2, 3, figsize=(6, 4))
phd.viz.despine(ax.ravel())
for a in ax.ravel():
    a.xaxis.set_tick_params(labelsize=6)
    a.yaxis.set_tick_params(labelsize=6)
    a.set_xscale('symlog', linthreshx=1E-2)
    a.set_xlabel('IPTG [µM]', fontsize=8)
    a.set_xlim([-.001, 1E4])

for i in range(3):
    ax[0, i].set_ylim([-0.01, 1.15])
    ax[1, i].set_ylim([-8, 8])
for i in range(2):
    ax[0, i+1].set_yticklabels([])
    ax[1, i+1].set_yticklabels([])
ax[0, 0].set_ylabel('fold-change', fontsize=8)
ax[1, 0].set_ylabel(r'$\Delta F$ [$k_BT$]', fontsize=8)

# Define the axes. 
axes = {'Q21M':0, 'Q21A':1, 'Y20I':2}
for m, a in axes.items():
    phd.viz.titlebox(ax[0, a], m, bgcolor='white', color=colors['black'], 
                    pad=0.05, boxsize='12%')

# ##############################################################################
#  INDUCTION PROFILES
# ##############################################################################
for g, d in DNA_data.groupby(['repressors', 'mutant']):
    _ax = ax[0, axes[g[1]]]
    epRA = epRA_stats[(epRA_stats['parameter']=='ep_RA') &
                      (epRA_stats['mutant']==g[1])][['hpd_min', 'hpd_max']].values[0]
    ep, c = np.meshgrid(epRA, c_range)
    arch = phd.thermo.SimpleRepression(R=g[0], ep_r=ep, effector_conc=c, 
                        ka=constants['Ka'], ki=constants['Ki'], 
                        ep_ai=constants['ep_AI']).fold_change()
    _ax.fill_between(c_range, arch[:, 0], arch[:, 1], color=rep_colors[g[0]],
                    alpha=0.4, label='__nolegend__')

for i in range(3):
    ax[1, i].hlines(0, -0.001, 1E4, 'k', linestyle='--', label='__nolegend__')

# ##############################################################################
#  DELTA F
# ##############################################################################
for g, d in epRA_stats.groupby('mutant'):
    _ax = ax[1, axes[g]]
    epRA = d[(d['parameter']=='ep_RA')][['hpd_min', 'hpd_max']].values[0]
    _ax.fill_between(c_range, epRA[0] - constants['O2'], 
                     epRA[1] - constants['O2'], color='slategray',
                     label='__nolegend__')
# ##############################################################################
# DATA
# ##############################################################################
for g, d in DNA_data.groupby(['repressors', 'mutant']):
    _ax = ax[0, axes[g[1]]]
    _ax.errorbar(d['IPTGuM'], d['mean'], d['sem'], lw=0.75, capsize=1, 
                color=rep_colors[g[0]], label=int(g[0]), linestyle='None',
                fmt='o', ms=4.5, markeredgewidth=0.5, markeredgecolor='white')

# ##############################################################################
# DELTA F DATA
# ##############################################################################
for g, d in empirical_F.groupby(['mutant', 'repressors', 'IPTGuM']):
    mu = d[d['parameter']=='fc_mu']['median'].values[0]
    sig = d[d['parameter']=='fc_sigma']['median'].values[0]
    dF = d[d['parameter']=='delta_bohr']
    _ax = ax[1, axes[g[0]]]
    if (mu > sig) & (1 - mu > sig):
        _ax.plot(dF['IPTGuM'], dF['median'], 'o', color=rep_colors[g[1]], ms=4.5, 
                 linestyle='None', markeredgecolor='white', markeredgewidth=0.5)
        _ax.vlines(dF['IPTGuM'], dF['hpd_min'], dF['hpd_max'], lw=0.75, 
                   color=rep_colors[g[1]])

# ##############################################################################
# LEGEND AND SAVING
# ##############################################################################
leg = ax[0, 0].legend(title='rep. / cell', fontsize=6)
leg.get_title().set_fontsize(6)
plt.tight_layout()
plt.savefig('../figs/ch7_figS21.pdf', bbox_inches='tight')
plt.savefig('../figs/ch7_figS21.png', bbox_inches='tight')

# %%
