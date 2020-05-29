# -*- coding: utf-8 -*-
#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import phd.viz
import phd.thermo
import phd.stats
import seaborn as sns
import pickle
import seaborn as sns
colors, palette = phd.viz.phd_style()
constants = phd.thermo.load_constants()
viridis = sns.color_palette('viridis', n_colors=4)

# Load data and correct names
data = pd.read_csv('../../src/data/ch2_induction/RazoMejia_2018.csv', comment='#')
data = data[data['repressors']> 0].copy()
data['repressors'] *= 2
data.rename(columns={'IPTG_uM':'IPTGuM'}, inplace=True)

# Identifier as to plot all of the data or not. 
ALL_DATA = 1
FIT_STRAIN = ['O2', 260]

# Compute the summary data
summary = data.groupby(['operator', 'repressors', 
            'IPTGuM']).agg(('mean','sem')).reset_index()

# ##############################################################################
# FIGURE INSTANTIATION
# ##############################################################################
# Set up the figure. 
with open('../../src/data/ch2_induction/mcmc/main_text_KaKi.pkl', 'rb') as pkl:
    chain = pickle.load(pkl)
ka_chain = np.exp(-chain[:, 0])[::10]
ki_chain = np.exp(-chain[:, 1])[::10]

# ##############################################################################
# FIGURE INSTANTIATION
# ##############################################################################
#%%
fig, ax = plt.subplots(2, 3, figsize=(7, 4.5), dpi=100)
ax[0,0].set_visible(False)
ax = [ax[0, 1], ax[0, 2], ax[1, 0], ax[1, 1], ax[1, 2]]

for a in ax:
    a.set_xscale('log')
    a.set_xlabel('repressors per cell')
    phd.viz.despine(ax)
ax[0].set_yscale('log')
ax[3].set_yscale('log')

# Add specific y axis labels.
ax[0].set_ylabel('leakiness')
ax[1].set_ylabel('saturation')
ax[2].set_ylabel('dynamic range')
ax[3].set_ylabel('$[EC_{50}]$ [µM]')
ax[4].set_ylabel('effective Hill coefficient')

# Fix the axis limits
ax[0].set_ylim([1E-4, 1.25])
ax[1].set_ylim([1E-4, 1.02])
ax[2].set_ylim([1E-4, 1.02])
ax[3].set_ylim([1E-2, 500])
ax[4].set_ylim([1.15, 1.95])

# Define teh operator colors.
op_colors = {'O1':viridis[0], 'O2':viridis[1], 'O3':viridis[2]}

# Define the property axes 
axes = {'leakiness':ax[0], 'saturation':ax[1], 'dynamic_range':ax[2],
        'EC50':ax[3], 'effective_hill':ax[4]}
# # ##############################################################################
# #  THEORY CURVES
# # ##############################################################################
rep_range = np.logspace(0, 4, 200)

for i, o in enumerate(['O1', 'O2', 'O3']):
    # Compute the modes
    arch = phd.thermo.SimpleRepression(R=rep_range, ep_r=constants[o],
                               ka=constants['Ka'], ki=constants['Ki'],
                               effector_conc=0, 
                               ep_ai=constants['ep_AI']).compute_properties()
    for k, v in arch.items():
        axes[k].plot(rep_range, v, color=op_colors[o], label=o)

    # Credible regions
    _ = np.zeros((2, len(rep_range)))
    cred_regions = {'saturation': _.copy(),
                    'dynamic_range':_.copy(), 'EC50':_.copy(),
                    'effective_hill':_.copy()}
    for j, r in enumerate(rep_range):
        arch = phd.thermo.SimpleRepression(R=r,ep_r=constants[o], 
                                           ka=ka_chain, ki=ki_chain,
                                           ep_ai=constants['ep_AI'],
                                           effector_conc=0).compute_properties()
        for k, v in arch.items():
            if k != 'leakiness':
                cred_regions[k][:, j] = phd.stats.compute_hpd(v, 0.95)
    for k, v in cred_regions.items():
        axes[k].fill_between(rep_range, v[0, :], v[1, :], color=op_colors[o],
                            alpha=0.3, label='__nolegend__')

# # ##############################################################################
# # DATA
# # ##############################################################################
for g, d in data.groupby(['repressors', 'operator']):
        leakiness = d[d['IPTGuM']==0]['fold_change_A'].values
        sat =d[d['IPTGuM']==d['IPTGuM'].max()]['fold_change_A'].values
        try:
            dyn_rng = sat - leakiness  
        except:
            dyn_rng = sat - leakiness[:-1]
        ax[0].errorbar(g[0], np.mean(leakiness), 
                    np.std(leakiness)/len(leakiness), lw=0.75, capsize=1, 
                    color=op_colors[g[-1]], label='__nolegend__', linestyle='none',
                    ms=4, fmt='o', markerfacecolor=op_colors[g[1]], 
                    markeredgecolor=colors['grey'], markeredgewidth=0.5) 
        ax[1].errorbar(g[0], np.mean(sat), np.std(sat)/len(sat), lw=0.75, 
                    capsize=1, color=op_colors[g[-1]], label='__nolegend__',
                    linestyle='none', ms=4, fmt='o', markeredgecolor=colors['grey'],
                    markerfacecolor=op_colors[g[1]], markeredgewidth=0.5) 
        ax[2].errorbar(g[0], np.mean(dyn_rng), np.std(dyn_rng)/len(dyn_rng), 
                          lw=0.75, capsize=1, color=op_colors[g[1]], 
                          label='__nolegend__', linestyle='none', ms=4, fmt='o',
                          markerfacecolor=op_colors[g[1]], 
                          markeredgecolor=colors['grey'], markeredgewidth=0.5) 

        # Inferred points
        with open(f'../../src/data/ch2_induction/mcmc/SI_I_{g[-1]}_R{int(g[0])}.pkl', 'rb') as pkl:
            chain = pickle.load(pkl)
        ka_chain = np.exp(-chain[:,0])
        ki_chain = np.exp(-chain[:,1])
        ka_median = np.median(ka_chain)
        ki_median = np.median(ki_chain)

        # Compute the inferred EC50 and effective hill
        arch = phd.thermo.SimpleRepression(R=g[0], ep_r=constants[g[1]],
                                           ka=ka_median, ki=ki_median, 
                                           ep_ai=constants['ep_AI'],
                                           effector_conc=0).compute_properties()
        ax[3].plot(g[0], arch['EC50'], 's', color=op_colors[g[1]], ms=4,
        markerfacecolor=op_colors[g[1]], markeredgewidth=0.5, markeredgecolor=colors['grey'])
        ax[4].plot(g[0], arch['effective_hill'], 's', color=op_colors[g[1]], 
                   ms=4, markerfacecolor=op_colors[g[1]], 
                   markeredgewidth=0.5, markeredgecolor=colors['grey'])

        # Compute the credible regions:
        arch = phd.thermo.SimpleRepression(R=g[0], ep_r=constants[g[1]],
                                           ka=ka_chain, ki=ki_chain, 
                                           ep_ai=constants['ep_AI'],
                                           effector_conc=0).compute_properties()
        ec50_cred = phd.stats.compute_hpd(arch['EC50'], 0.95)
        hill_cred = phd.stats.compute_hpd(arch['effective_hill'], 0.95)
    
        # Plot the vlines. 
        ax[3].vlines(g[0], ec50_cred[0], ec50_cred[1], lw=0.75,
                    color=op_colors[g[1]])
        ax[4].vlines(g[0], hill_cred[0], hill_cred[1], lw=0.75,
                    color=op_colors[g[1]])

# # ##############################################################################
# # LEGEND AND SAVING
# # ##############################################################################
leg = ax[0].legend(title='operator')
leg.get_title().set_fontsize(6)
plt.tight_layout()
plt.savefig('../figs/properties_data.pdf', bbox_inches='tight')


#

# %%
