#%%
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import phd.viz
import phd.stats
colors, color_list = phd.viz.phd_style()

# Load an example data set. 
fluct_data = pd.read_csv('../../data/ch4_growth/analyzed_fluctuations.csv', comment='#')
fluct_data = fluct_data[fluct_data['date']==20181002]
alpha_mean = np.round(fluct_data['alpha_mean'].unique(), -1)[0]
alpha_std = np.round(fluct_data['alpha_std'].unique())[0]

#%%
# Compute bins of 50 events each. 
BIN_SIZE = 75
bins = phd.stats.bin_by_events(fluct_data, average=['summed', 'fluct'], bin_size=BIN_SIZE)

# Compute the theory line
summed_range = np.logspace(1, 5.5, 200)

fig, ax = plt.subplots(1, 1, figsize=(3, 2))
phd.viz.despine(ax)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim([1E2, 10**5.5])
ax.set_xlabel(r'$I_1 + I_2$', fontsize=8)
ax.set_ylabel(r'$(I_1 - I_2)^2$', fontsize=8)
ax.plot(fluct_data['summed'], fluct_data['fluct'], '.', color=colors['black'], ms=0.75, label='division event')
ax.plot(summed_range, alpha_mean * summed_range, lw=0.75, color=colors['red'], 
        label=r'$\alpha$ ' + f'= {int(np.round(alpha_mean, decimals=-1))} ± {int(np.round(alpha_std, decimals=-1))} [a.u. / LacI]')
ax.plot(bins['summed'], bins['fluct'], marker='o', color=colors['red'], 
        markeredgecolor=colors['grey'], ms=4, label=f'{BIN_SIZE} events / bin',
        markeredgewidth=0.75, linestyle='none' )
ax.fill_between(summed_range, (alpha_mean + alpha_std) * summed_range, 
                (alpha_mean - alpha_std) * summed_range, color=colors['red'], alpha=0.25)
ax.legend(loc='lower right', fontsize=6)
phd.viz.titlebox(ax, 'glucose, 37 °C', size=6, color=colors['black'], bgcolor='white')
plt.savefig('../figs/Fig2_fluct_example.svg', bbox_inches='tight')


# %%
