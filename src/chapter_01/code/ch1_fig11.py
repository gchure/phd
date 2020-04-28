#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import phd.viz
import phd.thermo
colors, palette = phd.viz.phd_style()
constants = phd.thermo.load_constants()

# Load the data sets.
data = pd.read_csv('../../data/ch4_growth/analyzed_foldchange.csv', comment='#')
entropy = pd.read_csv('../../data/ch4_growth/pooled_entropic_parameter_summary.csv', comment='#')
old_data = pd.read_csv('../../data/ch4_growth/Garcia2011_Brewster2014_RazoMejia2018_Chure2019.csv', comment='#')
# old_data = old_data[old_data['IPTGuM']==0]

# Clean the data. 
data = data[(data['size']=='large') & (data['fold_change'] > 0) & (data['repressors'] > 0)]
data['repressors'] *= 1.16

# Summarize the data. 
summary = data.groupby(['atc_ngml', 'date', 
                        'carbon', 'temp'])[
                        ['repressors', 'fold_change']].mean().reset_index()
summary = summary.groupby(['atc_ngml', 'carbon', 'temp'])[
                            ['repressors', 'fold_change']
                            ].agg(('mean', 'sem')).reset_index()


#%%
# Set up the figure canvas
fig, ax = plt.subplots(2, 2, figsize=(6, 4.4))
phd.viz.despine(ax.ravel())

# Format the axes.
for i in range(2):
    ax[i, 0].set_xscale('log')
    ax[i, 0].set_yscale('log')
    ax[i, 1].set_yscale('log')
    ax[i, 0].set_xlabel('repressors per cell')
    ax[i, 1].set_xlabel('free energy [$k_BT$]')
    ax[i, 0].set_ylabel('fold-change')
    ax[i, 1].set_ylabel('fold-change')
    ax[i, 1].set_xlim([-5, 2])
    ax[i, 1].set_ylim([5E-3, 1.5])

# Compute and plot the collapse curve. 
bohr_range = np.linspace(-10, 10, 200)
collapse = (1 + np.exp(-bohr_range))**-1
for i in range(2):
    ax[i, 1].plot(bohr_range, collapse, 'k-', label='scaling function',
                zorder=10)

    # Plot the collapsed old data for reference
    ax[i, 1].plot(old_data['bohr_parameter'], old_data['mean'], 'o', ms=4, 
                color='lightgrey', zorder=9, label='previous data')



# Plot the data 
condition_colors = {'glucose':colors['purple'], 'glycerol':colors['green'],
                    'acetate':colors['brown'], 37:colors['purple'],
                     32:colors['blue'], 42:colors['red']}

# Plot the prediction for the carbon sources. 
rep_range = np.logspace(0, 3, 300)
fc = phd.thermo.SimpleRepression(rep_range, ep_r=constants['O2'],
                                ka=constants['Ka'], ki=constants['Ki'],
                                ep_ai=constants['ep_AI'], effector_conc=0).fold_change()
ax[0, 0].plot(rep_range, fc, 'k-', lw=1, label='prediction')

# Plot the various predictions for the temperature data.
temps = [32, 42]
for t in temps:

    # Compute the simple rescaling
    scale_factor = (37 + 273.15) / (t + 273.15)
    fc_simple = phd.thermo.SimpleRepression(rep_range, ep_r=scale_factor * constants['O2'],
                                    ka=constants['Ka'], ki=constants['Ki'],
                                    ep_ai=scale_factor * constants['ep_AI'],
                                    effector_conc=0).fold_change()

    # compute the proper rescaling
    params = entropy[entropy['temp']==t]
    ep_r = params[params['parameter']=='epRA_star']['median'].values[0]
    ep_ai = params[params['parameter']=='epAI_star']['median'].values[0]
    fc_entropy = phd.thermo.SimpleRepression(rep_range, ep_r=ep_r, ka=constants['Ka'],
                                    ki=constants['Ki'], ep_ai=ep_ai,
                                    effector_conc=0).fold_change()

    # Plot the two cases. 
    ax[1, 0].plot(rep_range, fc_simple, ':', color=condition_colors[t],
            lw=0.75, label='__nolegend__')
    ax[1, 0].plot(rep_range, fc_entropy, '-', color=condition_colors[t],
            lw=0.75, label='__nolegend__')

# Plot the data sources
for g, d in summary.groupby(['carbon', 'temp']):
    if (g[1] == 37):
        _index = 0 
        c = condition_colors[g[0]]

        # Compute the bohr parameter
        bohr = phd.thermo.SimpleRepression(d['repressors']['mean'], ep_r=constants['O2'],
                                           ka=constants['Ka'], ki=constants['Ki'],
                                           ep_ai=constants['ep_AI'], 
                                           effector_conc=0).bohr_parameter()
    else:
        _index = 1    
        c = condition_colors[g[1]]

        # Compute the bohr parameter
        params = entropy[entropy['temp']==g[1]]
        ep_r = params[params['parameter']=='epRA_star']['median'].values[0]
        ep_ai = params[params['parameter']=='epAI_star']['median'].values[0]
        bohr = phd.thermo.SimpleRepression(d['repressors']['mean'], ep_r=ep_r,
                                            ka=constants['Ka'], ki=constants['Ki'],
                                            ep_ai=ep_ai, effector_conc=0).bohr_parameter()

    # Plot the experimental measurements
    ax[_index, 0].errorbar(d['repressors']['mean'], d['fold_change']['mean'],
                           xerr=d['repressors']['sem'], yerr=d['fold_change']['sem'],
                           fmt='o', ms=4, color=c, markeredgecolor=colors['grey'],
                           markeredgewidth=0.5, label=f'{g[0]}, {g[1]}° C')

    # Plot the collapse data
    ax[_index, 1].errorbar(bohr, d['fold_change']['mean'], d['fold_change']['sem'],
                            fmt='o', color=c, markeredgecolor=colors['grey'],
                            ms=4.5, zorder=1000, label=f'{g[0]}, {g[1]}° C')

# Add legend as necessary
ax[1, 0].plot([], [], 'k:', label='simple rescaling')
ax[1, 0].plot([], [], 'k-', label='entropic penalty')
for a in ax.ravel():
    a.legend(fontsize=6)
ax[1, 1].legend(fontsize=6, loc='lower right')
for i in range(2):
    phd.viz.titlebox(ax[0, i], 'carbon quality variation', fontsize=6,
                    color=colors['black'], pad=0.05, boxsize='12%')
    phd.viz.titlebox(ax[1, i], 'temperature variation', fontsize=6,
                    color=colors['black'], pad=0.05, boxsize='12%')
plt.tight_layout()
fig.text(0, 0.95, '(A)', fontsize=8)
fig.text(0, 0.45, '(B)', fontsize=8)
plt.savefig('../figs/ch1_fig11.pdf', bbox_inches='tight')
plt.savefig('../figs/ch1_fig11.png', bbox_inches='tight')
# %%


# %%
