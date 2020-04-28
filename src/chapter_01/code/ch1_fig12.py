#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import phd.viz
import seaborn as sns
colors, palette = phd.viz.phd_style()

# Load the data. 
samples = pd.read_csv('../../data/ch5_mscl/complete_mcmc_traces.csv', comment='#')

# Define the channels over which to evaluate the probabilities
channel_range = np.logspace(0, 3.1, 300)

# Define the percentiles, labels, and colors
_colors = sns.color_palette('Purples_r', n_colors=10)
percs = {1:[45.5, 55.5, '1%', _colors[0], 100],
         5:[42.5, 57.5, '5%', _colors[1], 99],
        10:[45, 55, '10%', _colors[2], 98],
        25:[37.5, 62.5, '25%', _colors[3], 97],
        50:[25, 75, '50%', _colors[4], 96],
        75:[12.5, 87.5, '75%', _colors[5], 95],
        95:[2.5, 97.5, '95%', _colors[6], 94],
        99:[0.5, 99.5, '99%', _colors[7], 93]}

# Compute the percentiles
perc_df = pd.DataFrame([])
for i, r in enumerate([('fast', 1), ('slow', 0)]):
    beta_0 = samples[f'beta_0__{r[1]}'].values
    beta_1 = samples[f'beta_1__{r[1]}'].values
    for k in channel_range:
        prob = (1 + k ** -beta_1 * np.exp(-beta_0))**-1
        for p, l in percs.items():
            bottom, top, lab, c, z = l
            lower, upper = np.percentile(prob, l[:2])
            perc_df = perc_df.append({'shock_rate':r[0],
                                      'percentile':p,
                                      'label':lab,
                                      'color':c,
                                      'zorder':z,
                                      'lower_bound':lower,
                                      'upper_bound':upper,
                                      'channels':k}, ignore_index=True)
            
        
#%% Set up the figure canvas. 
fig, ax = plt.subplots(1, 2, figsize=(4.5, 2.5))
phd.viz.despine(ax.ravel())

# Format the axes
phd.viz.titlebox(ax[0], 'slow shock $<$ 1.0 Hz', fontsize=6, color=colors['black'],
                pad = 0.05, boxsize='12%')
phd.viz.titlebox(ax[1], 'fast shock $\geq$ 1.0 Hz', fontsize=6, color=colors['black'],
                pad = 0.05, boxsize='12%')

for i in range(2):
    ax[i].set_xlabel('channels per cell')
    ax[i].set_ylabel('survival probability')
    ax[i].set_ylim([0, 1])

shock_ax = {'slow':0, 'fast':1}
for g, d in perc_df.groupby('shock_rate'):
    for _g, _d in d.groupby(['percentile', 'color', 'zorder', 'label']):
        ax[shock_ax[g]].fill_between(_d['channels'], _d['lower_bound'], 
                                    _d['upper_bound'], label=_g[3], 
                                    color=_g[1], zorder=_g[2])


plt.tight_layout()
leg = ax[0].legend(loc='lower right', title='percentile', fontsize=6)
leg.get_title().set_fontsize(6)
plt.savefig('../figs/ch1_fig12_plots.svg', bbox_inches='tight')
# %%





# %%
34