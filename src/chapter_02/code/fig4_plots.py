#%%
import os
import glob
import pickle
import re
import numpy as np
import pandas as pd
import sys
import phd.viz
import phd.stats
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import matplotlib.colors as plc
import altair as alt
import seaborn as sns
import scipy.stats
colors, palette = phd.viz.phd_style()
DPI = 227
data = pd.read_csv('../../data/ch2_induction/RazoMejia_2018.csv', comment='#')

# Now we remove the autofluorescence and delta values
data = data[(data['rbs'] != 'auto') & 
            (data['rbs'] != 'delta') & 
            (data['operator']!='Oid')]

# Load the flat-chain
with open('../../data/ch2_induction/mcmc/SI_I_O2_R260.pkl', 'rb') as file:
    unpickler = pickle.Unpickler(file)
    gauss_flatchain = unpickler.load()
    gauss_flatlnprobability = unpickler.load()

# Map value of the parameters
max_idx = np.argmax(gauss_flatlnprobability, axis=0)
ea, ei, sigma = gauss_flatchain[max_idx]
lnprob =  gauss_flatchain[:, 2][::100]
ka_fc = np.exp(-gauss_flatchain[:, 0])[::100]
ki_fc = np.exp(-gauss_flatchain[:, 1])[::100]

#%%
# Define the IPTG concentrations to evaluate
c_range = np.logspace(-2, 4, 200)
c_range[0] = 0 
rep_range = np.logspace(0, 4, 200)

# Set the colors for the strains
rep_colors = {22: colors['light_red'], 
              60: colors['light_brown'],  
              124: colors['light_green'], 
              260: colors['light_orange'], 
              1220: colors['light_purple'], 
              1740: colors['light_blue']} 

# Define the operators and their respective energies
operators = ['O1', 'O2', 'O3']
energies = {'O1': -15.3, 'O2': -13.9, 'O3': -9.7}


#%%
# ##############################################################################
# SAMPLING JOINTPLOT
# ##############################################################################
sampling_df = pd.DataFrame(np.array([ka_fc, ki_fc, lnprob]).T, 
                          columns=['ka', 'ki', 'logprob'])
sampling_df.sort_values('logprob', inplace=True)

# Set up the relatively complex plot object. 
fig = plt.figure(figsize=(3, 2))
gs = gridspec.GridSpec(4, 4)
ax0 = fig.add_subplot(gs[1:, :3])
ax1 = fig.add_subplot(gs[0, :3])
ax2 = fig.add_subplot(gs[1:, 3])

for a in [ax1, ax2]:
    a.axis('off')

# Add appropriate axis labels. 
ax0.set_xlabel('$K_I$ [µM]')
ax0.set_ylabel('$K_A$ [µM]')


_ = ax0.plot(sampling_df['ki'], sampling_df['ka'], linestyle='none', marker='.', 
             color=colors['purple'], ms=1, alpha=0.75, markeredgewidth=0)
_ = ax1.hist(sampling_df['ki'].values, bins=25, density=True, 
            color=colors['light_purple'], edgecolor=colors['black'],
            lw=0.25)
_ = ax2.hist(sampling_df['ka'].values, bins=25, density=True, 
            color=colors['light_purple'], edgecolor=colors['black'],
            lw=0.25, orientation='horizontal')
plt.subplots_adjust(hspace=0.02, wspace=0.02)
plt.savefig('../figs/fig5_sampling_jointplot.svg', bbox_incehs='tight')

# ##############################################################################
# PREDICTED INDUCTION PROFILES
# ##############################################################################
# Set up the dataframe for the fold-change. 
fc_df = pd.DataFrame()

# Iterate through each operator, repressor, and IPTG concentration to calculate
# the fold-change
for op, op_en in energies.items():
    for r, _ in rep_colors.items():
        for i, c in enumerate(c_range):
            arch = phd.thermo.SimpleRepression(r, op_en, ka=ka_fc, ki=ki_fc, 
                                               ep_ai=4.5, effector_conc=c)
            
            fc_min, fc_max = phd.stats.compute_hpd(arch.fold_change(), 0.95)
            fc_df = fc_df.append({'fc_min': fc_min,
                                  'fc_max': fc_max,
                                  'repressors': r,
                                  'IPTGuM': c,
                                  'operator': op,
                                  'binding_energy': op_en}, ignore_index=True)

#  Iterate through each operator and compute the properties
prop_df = pd.DataFrame([])
for op, op_en in energies.items():
    for i, r in enumerate(rep_range):
        arch = phd.thermo.SimpleRepression(r, op_en, ka=ka_fc, ki=ki_fc, ep_ai=4.5,
                                           effector_conc=0).compute_properties()
        for prop, val in arch.items():
            if prop == 'leakiness':
                val_min = val 
                val_max = val 
            else:
                val_min, val_max = phd.stats.compute_hpd(val, 0.95)
            prop_df = prop_df.append({'val_min':val_min,
                              'val_max':val_max,
                              'property': prop,
                              'repressors': r,
                              'operator': op,
                              'binding_energy':op_en}, ignore_index=True)
        

#%%

# Set up the canvas for the operator plots. 
fig, ax = plt.subplots(1, 3, figsize=(6, 2))
axes={'O1':ax[0], 'O2':ax[1], 'O3':ax[2]}
for a in ax:
    a.set_xscale('symlog', linthreshx=1E-2)
    a.set_ylim([-0.01, 1.15])
    a.set_ylabel('fold-change')
    a.set_xlabel('IPTG [µM]')

for g, d in fc_df.groupby('operator'):
    for r, _d in d.groupby('repressors'):
        _ = axes[g].fill_between(_d['IPTGuM'], _d['fc_min'], _d['fc_max'],
                                 alpha=0.75, facecolor=rep_colors[r], label=int(r),
                                 edgecolor=colors['black'], lw=0.1)
    phd.viz.titlebox(axes[g], f'operator {g}', color=colors['black'], 
                    boxsize='12%')

# Isolate the data for the fit strain
fit_strain = data[(data['operator']=='O2') & 
                  (data['repressors']==130)].groupby(['IPTG_uM']
                  ).agg(('mean', 'sem')).reset_index()

# Plot the aggregated points of the fit strain
ax[1].errorbar(fit_strain['IPTG_uM'], fit_strain['fold_change_A']['mean'], 
               yerr=fit_strain['fold_change_A']['sem'], fmt='.', 
               markerfacecolor=colors['light_orange'], markeredgecolor=colors['black'],
               markeredgewidth=0.5,
               label='__nolegend__', lw=0.1, color=colors['orange'])

# Add the legend. 
leg = ax[0].legend(loc='upper left', title='rep. per cell', fontsize=6, handlelength=1.5)
leg.get_title().set_fontsize(6)
plt.tight_layout()
plt.savefig('../figs/fig5_induction_profiles.svg', bbox_inches='tight')

# ##############################################################################
#%% PROPERTIES
# ##############################################################################
fig, ax =  plt.subplots(2, 3, figsize=(6, 3.5))
ax[1, -1].axis('off')

# Assign the properties to axes. 
axes = {'leakiness':ax[0,0], 'saturation':ax[0, 1], 'dynamic_range':ax[0, 2],
        'EC50':ax[1,0], 'effective_hill':ax[1, 1]}

op_palette = sns.color_palette('viridis', n_colors=4)
op_colors = {'O1':op_palette[0], 'O2':op_palette[1], 'O3':op_palette[2]}

for g, d in prop_df.groupby(['property', 'operator']):
    _ax = axes[g[0]]
    if g[0] == 'leakiness':
        _ax.plot(d['repressors'], d['val_min'], 
                    color=op_colors[g[1]],
                    label=g[1], alpha=0.5)
    else:
        _ax.fill_between(d['repressors'], d['val_min'], d['val_max'], lw=0.1,
                    facecolor=op_colors[g[1]], edgecolor=colors['black'], 
                    label=g[1], alpha=0.5)

# Set the appropriate scaling. 
for a in ax.ravel():
    a.set_xscale('log')
    a.set_xlabel('repressors per cell')
for a in [ax[0, 0], ax[1, 0]]:
    a.set_yscale('log')

# Manually define the ylabels 
_ax = ax.ravel()
_ax[0].set_ylabel('leakiness')
_ax[1].set_ylabel('saturation')
_ax[2].set_ylabel('dynamic range')
_ax[3].set_ylabel('EC$_{50}$ [µM]')
_ax[4].set_ylabel('effective Hill coefficient')
leg = ax[0,0].legend(handlelength=1.5, title='operator', loc='lower left', fontsize=6)
leg.get_title().set_fontsize(6)
plt.tight_layout()
plt.savefig('../figs/fig5_properties.svg', bbox_inches='tight')
#%%
