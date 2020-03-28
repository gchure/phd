#%%
import os
import glob
import pickle
import re
import numpy as np
import pandas as pd
import phd.viz
import phd.stats
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn as sns
colors, palette = phd.viz.phd_style()

#===============================================================================
# Read the data
#===============================================================================
# Define working directory
datadir = '../../data/'
flow_data = pd.read_csv('../../data/ch2_induction/RazoMejia_2018.csv', comment='#')
oid_data = pd.read_csv('../../data/ch6_induction_si/oid_foldchange_measurements.csv', comment='#')
oid_data = oid_data[oid_data['fold_change'] <= 1]
flow_data.rename(columns={'fold_change_A':'fold_change'}, inplace=True)
df = pd.concat([flow_data, oid_data], sort=False)
df['repressors'] *= 2

# Now we remove the autofluorescence and delta values
df = df[(df.rbs != 'auto') & (df.rbs != 'delta')]

# COmpute the summary dataframe. 
summary = df.groupby(['repressors', 'operator', 
                    'IPTG_uM'])['fold_change'].agg(('mean', 'sem')).reset_index()

# Load the flat-chain
with open('../../data/ch6_induction_si/SI_E_global.pkl', 'rb') as file:
    unpickler = pickle.Unpickler(file)
    gauss_flatchain = unpickler.load()
    gauss_flatlnprobability = unpickler.load()

# Generate a Pandas Data Frame with the mcmc chain
columns = np.concatenate([['ea', 'ei', 'sigma'],
                          [df[df.repressors == r].rbs.unique()[0] for r in
                           np.sort(df.repressors.unique())],
                          [df[df.binding_energy == o].operator.unique()[0] for o in
                           np.sort(df.binding_energy.unique())]])

mcmc_df = pd.DataFrame(gauss_flatchain, columns=columns)
mcmc_df['ka'] = np.exp(-mcmc_df['ea'].values)
mcmc_df['ki'] = np.exp(-mcmc_df['ei'].values)
mcmc_df.rename(columns={'HG104':'R22', 'RBS1147':'R60', 'RBS446':'R124',
                         'RBS1027':'R260', 'RBS1':'R1220', 'RBS1L':'R1740'},
                         inplace=True)
# Define the modes (NOTE THIS IS TAKEN FROM THE INDUCTION SI DIRECTLY). 
modes = {'O1':-15.2, 'O2':-13.6, 'O3':-9.4, 'Oid':-17.7}

#%%
# Define the IPTG concentrations to evaluate
c_range = np.logspace(-2, 4, 200)
rep_colors = {22: colors['red'], 
              60: colors['brown'],  
              124: colors['green'], 
              260: colors['orange'], 
              1220: colors['purple'], 
              1740: colors['blue']}

# Define the operators and their respective energies
operators = ['O1', 'O2', 'O3', 'Oid']
energies = {'O1': -15.3, 'O2': -13.9, 'O3': -9.7, 'Oid': -17.0}

# Initialize the plot to set the size
fig, ax = plt.subplots(2, 2, figsize=(6, 4))
ax = ax.ravel()
phd.viz.despine(ax)

# Format the axes
for a in ax:
    a.set_xscale('log')
    a.set_xlabel('IPTG [µM]')
    a.set_ylabel('fold-change')
    a.set_xlim([1E-2, 1E4])
    a.set_ylim([-0.01, 1.2])

# Map axes to operators
op_ax = {'O1':ax[0], 'O2':ax[1], 'O3':ax[2], 'Oid':ax[3]}

# Loop through operators
for g, d in summary.groupby(['repressors', 'operator']):
    # Compute the credible region. 
    cred_region = np.zeros((2, len(c_range)))
    for i, c in enumerate(c_range):
        arch = phd.thermo.SimpleRepression(R=mcmc_df[f'R{g[0]}'].values[::5] * 2,
                                           ep_r=mcmc_df[g[1]].values[::5],
                                           ka=mcmc_df['ka'].values[::5],
                                           ki=mcmc_df['ki'].values[::5],
                                           ep_ai=4.5,
                                           effector_conc=c).fold_change()
        cred_region[:, i] = phd.stats.compute_hpd(arch, 0.95)
    op_ax[g[1]].fill_between(c_range, cred_region[0, :], cred_region[1, :],
                            color=rep_colors[g[0]], alpha=0.5, label=g[0])
    
    # Plot the data. 
    op_ax[g[1]].errorbar(d['IPTG_uM'], d['mean'], d['sem'], color=rep_colors[g[0]],
                        markeredgecolor='white', markeredgewidth=0.5, markersize=4,
                        label='__nolegend__', linestyle='none', linewidth=0.75,
                        fmt='o')
    
    # Add a title. 
    phd.viz.titlebox(op_ax[g[1]], f'{g[1]} ' + r'$\Delta\varepsilon_{RA} = $' + f'{modes[g[1]]} ' + r'$k_BT$',
                    color=colors['black'], bgcolor='white', size=8, boxsize='14%',
                    pad=0.05)

# Add a legend. 
leg = ax[0].legend(fontsize=6, title='repressors per cell')
leg.get_title().set_fontsize(6)
plt.tight_layout()
fig.text(0, 1, '(A)', fontsize=8)
fig.text(0.5, 1, '(B)', fontsize=8)
fig.text(0, 0.5, '(C)', fontsize=8)
fig.text(0.5, 0.5, '(D)', fontsize=8)
plt.savefig('../figs/ch6_figS12.pdf', bbox_inches='tight')
plt.savefig('../figs/ch6_figS12.png', bbox_inches='tight')

# %%
