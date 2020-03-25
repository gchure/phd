#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import phd.viz
import phd.thermo
constants = phd.thermo.load_constants()
colors, palette = phd.viz.phd_style()


# %%
# Load induction fold-change data. 
data = pd.read_csv('../../data/ch2_induction/RazoMejia_2018.csv', comment='#')
data['repressors'] *= 2
data = data[data['repressors'] > 0]

# Compute the summarized values.
data = data.groupby(['repressors', 'operator', 
                    'IPTG_uM'])['fold_change_A'].agg(('mean', 'sem')).reset_index()

#%%
rep_colors = {22: colors['red'], 
              60: colors['brown'],  
              124: colors['green'], 
              260: colors['orange'], 
              1220: colors['purple'], 
              1740: colors['blue']} 
op_glyphs = {'O1':'^', 'O2':'o', 'O3':'v'}
c_range = np.logspace(-2, 4, 200)
bohr_range = np.linspace(-9, 10, 200)
collapse = (1 + np.exp(-bohr_range))**-1
# %%
fig, ax = plt.subplots(1, 2, figsize=(4, 2))
phd.viz.despine(ax)
ax[0].set_xscale('log')

# Give proper labels to things
ax[0].set_xlabel('IPTG [µM]')
ax[0].set_ylabel('fold-change in gene expression')
ax[1].set_ylabel('fold-change in gene expression')
ax[1].set_xlabel('free energy [$k_BT$]')

# Plot the master collapse curve
ax[1].plot(bohr_range, collapse, 'k-', label='__nolegend__', lw=0.5)

# Iterate through and plot the data. 
for g, d in data.groupby(['repressors', 'operator']):
    # Plot the data
    _ = ax[0].errorbar(d['IPTG_uM'], d['mean'], d['sem'], 
                   fmt=op_glyphs[g[1]], ms=3, color=rep_colors[g[0]],
                    markeredgecolor='white', markeredgewidth=0.5,
                    linestyle='none', linewidth=0.65)

    # Compute and plot the bohr parameter
    bohr = phd.thermo.SimpleRepression(R=g[0], ep_r=constants[g[1]], 
                                        ka=constants['Ka'], ki=constants['Ki'],
                                        ep_ai=constants['ep_AI'],
                                        effector_conc=d['IPTG_uM'].values).bohr_parameter()
    _ = ax[1].errorbar(bohr, d['mean'], d['sem'], fmt=op_glyphs[g[1]],
                color=rep_colors[g[0]], markeredgecolor='white', 
                markeredgewidth=0.5, label='__nolegend__', ms=3)

    # Plot the 'fit'
    fit = phd.thermo.SimpleRepression(R=g[0], ep_r=constants[g[1]],
                                      ka=constants['Ka'], ki=constants['Ki'],
                                      ep_ai=constants['ep_AI'], 
                                      effector_conc=c_range).fold_change()
    _ = ax[0].plot(c_range, fit, color=rep_colors[g[0]], lw=0.5, label='__nolegend__')
plt.savefig('../figs/ch1_figX_plots.svg', bbox_inches='tight')
# %%
