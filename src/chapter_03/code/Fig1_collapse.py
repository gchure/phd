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
ax[2].set_xlabel('Bohr parameter [$k_BT$]', fontsize=8)
ax[2].set_ylabel('fold-change', fontsize=8)
ax[0].set_xscale('log')
ax[0].set_yscale('log')
ax[1].set_xscale('log')
ax[1].set_ylabel('fold-change', fontsize=8)
ax[0].set_ylabel('leakiness', fontsize=8)
ax[0].set_xlabel('repressors per cell', fontsize=8)
ax[1].set_xlabel('IPTG [M]', fontsize=8)
ax[1].set_ylabel('fold-change', fontsize=8)
ax[0].set_xlim([1, 1E4])
ax[1].set_xlim([1E-8, 1E-2])
# Plot the master curve. 
_ = ax[2].plot(bohr_range, collapse, 'k-', label='__nolegend__', lw=1)

# Iterate through the brewster and garcia data and compute the bohr parameter.
_color = {'O1':colors['red'], 'O2':colors['purple'], 'O3':colors['blue'], 'Oid':'k'}
glyphs = {'garcia':'o', 'brewster':'D'}

for g, d in old_gods.groupby(['author', 'operator', 'repressor']):
    # Compute and plot bohrconstants.
    bohr = phd.thermo.SimpleRepression(R=g[2], ep_r=constants[g[1]], ka=constants['Ka'],
                                      ki=constants['Ki'], effector_conc=0, n_sites=constants['n_sites'], 
                                       n_ns=constants['Nns'], ep_ai=constants['ep_AI']).bohr_parameter()
    _ = ax[2].plot(bohr, d['fold_change'], glyphs[g[0]], color=_color[g[1]], ms=1.5, label='__nolegend__', alpha=0.75)
    
    # Plot the leakiness values.
    _ = ax[0].plot(g[-1], d['fold_change'], glyphs[g[0]], color=_color[g[1]], ms=1.5, label='__nolegend__', alpha=0.75)

# Plot the flow data. 
for g, d in new_gods.groupby(['operator', 'IPTG_uM', 'repressors']):
    bohr = phd.thermo.SimpleRepression(R=2 * g[2], ep_r=constants[g[0]], ka=constants['Ka'],
                                      ki=constants['Ki'], ep_ai=constants['ep_AI'],
                                      n_sites=constants['n_sites'], n_ns=constants['Nns'],
                                      effector_conc=g[1]).bohr_parameter()
    _ = ax[2].plot(bohr, d['fold_change_A']['mean'], 'o', color=_color[g[0]], ms=1.5, label='__nolegend__', alpha=0.75, zorder=0)
    
    if g[1] == 0:
        _ = ax[0].plot(2 * g[2], d['fold_change_A']['mean'], 'x', color=_color[g[0]], ms=1.5, label='__nolegend__')
    
    else:
        _ = ax[1].plot(g[1] / 1E6, d['fold_change_A']['mean'], 'x', color=_color[g[0]], ms=1.5, label='__nolegend__')
        
# Plot the theory curves for each. 
for i, o in enumerate(_color.keys()):
    # Plot the leakiness
    rep_range = np.logspace(0, 4, 200)
    arch = phd.thermo.SimpleRepression(R=rep_range, ep_r=constants[o], ep_ai=constants['ep_AI'],
                                      ka=constants['Ka'], ki=constants['Ki'], effector_conc=0,
                                      n_sites=constants['n_sites'], n_ns=constants['Nns']).fold_change()
    _ = ax[0].plot(rep_range, arch, '-', lw=1, color=_color[o], label=o)
    
    # Plot the induction profiles. 
    if o != 'Oid':
        c_range = np.logspace(-2, 4, 200)
        for r in new_gods['repressors'].unique():
            arch = phd.thermo.SimpleRepression(R=2 * r, ep_r=constants[o], ep_ai=constants['ep_AI'],
                                      ka=constants['Ka'], ki=constants['Ki'], effector_conc=c_range,
                                      n_sites=constants['n_sites'], n_ns=constants['Nns']).fold_change()
            _ = ax[1].plot(c_range / 1E6, arch, '-', lw=1, color=_color[o], label='__nolegend__')
    
# # Add the appropriate labels
ax[2].plot([], [], 'o', ms=2, color='slategray', label='Garcia & Phillips\n2011')
ax[2].plot([], [], 'D', ms=2, color='slategray', label='Brewster et al.\n2014')
ax[2].plot([], [], 'x', ms=2, color='slategray', label='Razo-Mejia et al.\n2018')

ax[-1].legend(fontsize=6, handletextpad=0.1, handlelength=1, loc='upper left')
ax[0].legend(fontsize=6, handletextpad=0.2,  handlelength=1, loc='lower left')
plt.tight_layout()
# plt.savefig('Fig1_collapse.svg')



# %%
