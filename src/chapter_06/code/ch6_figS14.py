#%%
import os
import glob
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import phd.viz
import phd.thermo
import phd.stats
import seaborn as sns
colors, palette = phd.viz.phd_style()
constants = phd.thermo.load_constants()

data = pd.read_csv('../../data/ch6_induction_si/oid_foldchange_measurements.csv',
                   comment='#')
data['repressors'] *= 2
data = data[(data['repressors'] > 0) & (data['fold_change'] <= 1)]

data = data.groupby(['operator', 
'repressors', 'IPTG_uM'])['fold_change'].agg(('mean', 'sem')).reset_index()

# Load the flat-chain
with open('../../data/ch6_induction_si/main_text_KaKi.pkl', 'rb') as file:
    unpickler = pickle.Unpickler(file)
    gauss_flatchain_O2 = unpickler.load()
    gauss_flatlnprobability_O2 = unpickler.load()

# map value of the parameters
ka_main = np.exp(-gauss_flatchain_O2[:, 0])
ki_main = np.exp(-gauss_flatchain_O2[:, 1])

with open('../../data/ch6_induction_si/SI_E_global.pkl', 'rb') as file:
    unpickler = pickle.Unpickler(file)
    flatchain = unpickler.load()

# Set the indexing for the MCMC dataframe.
index = ['log_ka', 'log_ki', 'sigma', 'R22', 'R60', 'R124',
         'R260', 'R1220', 'R1740', 'Oid', 'O1', 'O2', 'O3']

global_df = pd.DataFrame(flatchain, columns=index)
global_df['Ka'] = np.exp(-global_df['log_ka'])
global_df['Ki'] = np.exp(-global_df['log_ki'])
index = global_df.columns


#%%
#===============================================================================
# Plot the theory vs data for all 4 operators with the credible region
#===============================================================================
# Define the IPTG concentrations to evaluate
c_range = np.logspace(-2, 4, 200)

rep_colors = {22: colors['red'], 60:colors['brown'], 124:colors['green']}

# Initialize the plot to set the size
fig, ax = plt.subplots(1, 2, figsize=(5, 2))
phd.viz.despine(ax)

# Format the axes as needed and plot the data
title_energy = [-17.0, -17.7]
for i, a in enumerate(ax):
    a.set_xlabel('IPTG [µM]')
    a.set_ylabel('fold-change')
    a.set_xscale('log')
    phd.viz.titlebox(a, r'$\Delta\varepsilon_{RA} = %s\, k_BT$' %str(title_energy[i]),
                    bgcolor='white', color=colors['black'], size=6, pad=0.05,
                    boxsize='12%')

# Plot the predictions. 
for g, d in data.groupby(['repressors']):
    old_cred_region = np.zeros((2, len(c_range)))
    new_cred_region = np.zeros((2, len(c_range)))
    for i, c in enumerate(c_range):
        old_arch = phd.thermo.SimpleRepression(R=g, ep_r=-17,
                                               ka=ka_main[::5], ki=ki_main[::5],
                                               ep_ai=4.5, effector_conc=c).fold_change()
        new_arch = phd.thermo.SimpleRepression(R=global_df[f'R{g}'].values[::5] * 2, 
                                               ep_r=global_df['Oid'].values[::5],
                                               ka=global_df['Ka'].values[::5],
                                               ki=global_df['Ki'].values[::5],
                                               ep_ai=4.5, effector_conc=c).fold_change()
        old_cred_region[:, i] = phd.stats.compute_hpd(old_arch, 0.95)
        new_cred_region[:, i] = phd.stats.compute_hpd(new_arch, 0.95)
    ax[0].fill_between(c_range, old_cred_region[0, :], old_cred_region[1, :],
                        color=rep_colors[g], alpha=0.5, label='__nolegend__') 
    ax[1].fill_between(c_range, new_cred_region[0, :], new_cred_region[1, :],
                        color=rep_colors[g], alpha=0.5, label='__nolegend__') 

for a in ax:
    for g, d in data.groupby(['repressors']):
        a.errorbar(d['IPTG_uM'], d['mean'], d['sem'], fmt='o', ms=4.5, 
                    markeredgewidth=0.5, markeredgecolor='white',
                    color=rep_colors[g], label=g)
        
leg = ax[0].legend(title='repressors per cell', fontsize=6)
leg.get_title().set_fontsize(6)
plt.tight_layout()
fig.text(0.01, 0.9, '(A)', fontsize=8)
fig.text(0.51, 0.9, '(B)', fontsize=8)
plt.savefig('../figs/ch6_figS14.pdf', bbox_inches='tight')
plt.savefig('../figs/ch6_figS14.png', bbox_inches='tight')
# %%
