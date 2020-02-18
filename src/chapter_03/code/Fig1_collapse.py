# -*- coding: utf-8 -*-
#%%
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import phd.thermo
import phd.stats
import phd.viz
constants = phd.thermo.load_constants()
constants.update({'Oid':-17.0})
colors, palette = phd.viz.phd_style()

colors = {c:i for c, i in colors.items() if ('light' not in c) & ('pale' not in c)}

# Load the data.
old_gods = pd.read_csv('../../data/other/Garcia2011_Brewster2014.csv', comment='#')
new_gods = pd.read_csv('../../data/ch2_induction/RazoMejia_2018.csv', comment='#')
new_gods = new_gods[(new_gods['repressors'] > 0) & (new_gods['repressors'] <= 130) & (new_gods['repressors'] >11)]

# Prune non-physical fold-change. 
old_gods = old_gods[(old_gods['fold_change'] >= -0.1) & (old_gods['fold_change'] <= 1.2)]
new_gods = new_gods[(new_gods['fold_change_A'] >= -0.1) & (new_gods['fold_change_A'] <= 1.2)]

# Summarize the induction data. 
new_gods.loc[:, 'fold_change'] = new_gods['fold_change_A']
new_gods = new_gods.groupby(['repressors', 'operator', 'IPTG_uM']).agg(('mean', 'sem')).reset_index()


# Plot the collapse curve. 
bohr_range = np.linspace(-10, 10, 200)
collapse = (1 + np.exp(-bohr_range))**-1

#%%
# Instantiate the figure.
fig, ax = plt.subplots(1, 3, figsize=(6, 2 ))
for a in ax.ravel():
    a.xaxis.set_tick_params(labelsize=8)
    a.yaxis.set_tick_params(labelsize=8)
for a in ax:
    a.grid(False)
ax[2].set_xlabel('Bohr parameter [$k_BT$]', fontsize=8)
ax[2].set_ylabel('fold-change', fontsize=8)
ax[0].set_xscale('log')
ax[1].axis('off')
ax[0].set_ylabel('fold-change', fontsize=8)
ax[0].set_xlabel('IPTG [µM]', fontsize=8)
ax[0].set_xlim([1E-2, 1E4])
ax[0].set_yticks([0, 0.5, 1.0])
ax[2].set_yticks([0, 0.5, 1.0])
ax[0].set_xticks([1E-1, 1E1, 1E3])
ax[2].set_xticks([-8, 0, 8])

# Plot the master curve. 
_ = ax[2].plot(bohr_range, collapse, 'k-', label='__nolegend__', lw=1)

# Iterate through the brewster and garcia data and compute the bohr parameter.
_color = {'O1':colors['red'], 'O2':colors['purple'], 'O3':colors['blue'], 'Oid':'k'}
glyphs = {'garcia':'o', 'brewster':'D'}

# Plot the flow data. 
iter = 0
for g, d in new_gods.groupby(['operator',  'repressors']):
    bohr = phd.thermo.SimpleRepression(R=2 * g[1], ep_r=constants[g[0]], ka=constants['Ka'],
                                      ki=constants['Ki'], ep_ai=constants['ep_AI'],
                                      n_sites=constants['n_sites'], n_ns=constants['Nns'],
                                      effector_conc=d['IPTG_uM']).bohr_parameter()
    _ = ax[2].plot(bohr, d['fold_change_A']['mean'], 'o', markerfacecolor=palette[iter], ms=2.5, 
                  label='__nolegend__', alpha=0.75, zorder=2, markeredgecolor='k',
                  markeredgewidth=0.25)
    _ = ax[0].plot(d['IPTG_uM'], d['fold_change_A']['mean'], 'o', markerfacecolor=palette[iter], ms=2.5, label='__nolegend__',
                   markeredgecolor='k', markeredgewidth=0.25)
    iter += 1

# Plot the theory curves for each. 
iter = 0
for i, o in enumerate(_color.keys()):
    # Plot the induction profiles. 
    if o != 'Oid':
        c_range = np.logspace(-2, 4, 200)
        for r in new_gods['repressors'].unique():
            arch = phd.thermo.SimpleRepression(R=2 * r, ep_r=constants[o], ep_ai=constants['ep_AI'],
                                      ka=constants['Ka'], ki=constants['Ki'], effector_conc=c_range,
                                      n_sites=constants['n_sites'], n_ns=constants['Nns']).fold_change()
            _ = ax[0].plot(c_range, arch, '-', lw=1, color=palette[iter], label='__nolegend__')
            iter += 1 

# # Add the appropriate labels
plt.tight_layout()
plt.savefig('../figs/fig1_collapse_plot.svg', bbox_inches='tight')



# %%
