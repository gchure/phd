#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import phd.viz
import seaborn as sns
colors, palette = phd.viz.phd_style()

# Load the data sets
data = pd.read_csv('../../data/ch1_introduction/compiled_absolute_measurements.csv')

#%%
fig, ax = plt.subplots(1, 2, figsize=(6, 3))
phd.viz.despine(ax)

# Format the axes
ax[0].set_ylabel('mass of proteome sector [fg / cell]')
ax[0].set_xlabel('growth rate [hr$^{-1}$]')
ax[1].set_xlabel('information storage and\nprocessing mass fraction')
ax[1].set_ylabel('metabolism mass fraction')

# Define the author glyphs
author_glyphs = {'li_2014':'^', 'peebo_2015':'v', 'schmidt_2016':'o', 'valgepea_2013':'D'}
cog_colors = {'cellular processes and signaling':colors['purple'],
              'metabolism': colors['blue'],
              'information storage and processing':colors['green'],
              'not assigned':colors['light_grey'],
              'poorly characterized':colors['black']}

summarized = data.groupby(['dataset_name', 'dataset', 'cog_class', 
                          'growth_rate_hr', 'condition'])['fg_per_cell'].agg('sum').reset_index()
for g, d in summarized.groupby(['cog_class', 'dataset', 'dataset_name']):
    ax[0].plot(d['growth_rate_hr'], d['fg_per_cell'], color=cog_colors[g[0].lower()],
                linestyle='none', marker=author_glyphs[g[1]], markersize=4.5,
                markeredgecolor=colors['grey'], markeredgewidth=0.5,
                alpha=0.75, label='__nolegend__')

# Plot the correlation of info storage and processing vs metabolism
sorted_growth = np.sort(summarized['growth_rate_hr'].unique())
growth_colors = {g:c for g, c in zip(sorted_growth, 
               sns.color_palette('viridis_r', n_colors=len(sorted_growth) + 5))}
for g, d in summarized.groupby(['dataset', 'dataset_name', 'condition', 'growth_rate_hr']):
    total_mass = d['fg_per_cell'].sum()
    info_frac = d[d['cog_class']=='information storage and processing']['fg_per_cell'] / total_mass
    metab_frac = d[d['cog_class']=='metabolism']['fg_per_cell'] / total_mass
    ax[1].plot(info_frac, metab_frac, linestyle='none', marker=author_glyphs[g[0]],
               color=growth_colors[g[-1]], markeredgecolor=colors['grey'],
               markeredgewidth=0.5, label='__nolegend__')

# Plot thing for legend entry
for g, d in summarized.groupby(['dataset', 'dataset_name']):
    ax[0].plot([], [], linestyle='none', marker=author_glyphs[g[0]], markerfacecolor=colors['grey'],
                   markeredgecolor=colors['black'], label=g[-1], ms=5)
for n, c in cog_colors.items():
    ax[0].plot([], [], '-', lw=2, color=c, label=n)

for g, c in growth_colors.items():
    ax[1].plot([], [], '-', lw=2, color=c, label=np.round(g, decimals=2))

ax[0].legend(fontsize=6, handlelength=1)

plt.tight_layout()
plt.savefig('../figs/ch1_fig13_plots.svg', bbox_inches='tight')

# %%
