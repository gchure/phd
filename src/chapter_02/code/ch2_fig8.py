#%%
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import phd.viz 
import phd.thermo
import phd.stats
import pickle
colors, palette = phd.viz.phd_style()
constants = phd.thermo.load_constants()

#%% Load the experimental data and summarize
data = pd.read_csv('../../data/ch2_induction/RazoMejia_2018.csv',
                    comment='#')
data['repressors'] *= 2
data = data[data['repressors'] > 0]
grouped = data.groupby(['operator', 'repressors', 
                        'IPTG_uM']
                        )['fold_change_A'].agg(('mean', 'sem')).reset_index()

#%% Load the sampling chains for Ka Ki estimates. 
# Load the flatchains for the prediction measurements. 
with open('../../data/ch2_induction/mcmc/SI_I_O2_R260.pkl', 'rb') as file:
    unpickler = pickle.Unpickler(file)
    gauss_flatchain = unpickler.load()
    gauss_flatlnprobability = unpickler.load()
ka_fc = np.exp(-gauss_flatchain[:, 0])[::100]
ki_fc = np.exp(-gauss_flatchain[:, 1])[::100]
# %%
# Set up the figure canvas
fig, ax = plt.subplots(1, 1, figsize=(2.5, 2), dpi=150)
phd.viz.despine(ax)
ax.set_xlabel('free energy [$k_BT$]')
ax.set_ylabel('fold-change')

# Define the master curve. 
bohr_range = np.linspace(-10, 10, 200)
master = (1 + np.exp(-bohr_range))**-1

_ = ax.plot(bohr_range, master, 'k-', lw=0.5, label='__nolegend__')

# Define the glyphs and repressor colors
rep_colors = {22: colors['red'], 
              60: colors['brown'],  
              124: colors['green'], 
              260: colors['orange'], 
              1220: colors['purple'], 
              1740: colors['blue']} 
op_glyphs = {'O1':'s', 'O2':'o', 'O3':'D'}

# Iterate through each condition and plot 
for g, d in grouped.groupby(['repressors', 'operator', 'IPTG_uM']):
    bohr_best = phd.thermo.SimpleRepression(R=g[0], 
                                            ep_r=constants[g[1]],
                                            ka=constants['Ka'], 
                                            ki=constants['Ki'],
                                            ep_ai=constants['ep_AI'],
                                            effector_conc=g[-1]).bohr_parameter()
    cred_region = phd.thermo.SimpleRepression(R=g[0],
                                              ep_r=constants[g[1]],
                                              ka=ka_fc,
                                              ki=ki_fc,
                                              ep_ai=constants['ep_AI'],
                                              effector_conc=g[-1]).bohr_parameter()
    low, high = phd.stats.compute_hpd(cred_region, 0.95)

    # Plot
    _ = ax.hlines(d['mean'], low, high, color=rep_colors[g[0]], lw=0.5, label='__nolegend__')
    _ = ax.errorbar(bohr_best, d['mean'], d['sem'], lw=0.5, fmt=op_glyphs[g[1]],
                    color=rep_colors[g[0]], markeredgecolor='white', label='__nolegend__',
                    ms=3, markeredgewidth=0.25)
# Add the legend
for r, c in rep_colors.items():
    ax.plot([], [], '-', lw=2, color=c, label=r)
for o, g in op_glyphs.items():
    ax.plot([], [], g, color='slategrey', markeredgecolor='white', 
            markeredgewidth=.25, ms=3, label=o)
ax.legend(fontsize=6, handlelength=0.25)
plt.savefig('../figs/fig7_collapse.svg', bbox_inches='tight')
# %%
