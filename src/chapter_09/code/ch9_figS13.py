
# %%
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import phd.viz
colors, palette = phd.viz.phd_style()
fig5b_data = pd.read_csv('../../data/ch9_mscl_si/poolman_fig5b.csv')
figS6b_data = pd.read_csv('../../data/ch9_mscl_si/van_den_Berg_2016_figS6b.csv')

#%%
# Define the error reported in van den Berg et al 2016 SI
y_err = 0.3

# Instantiate the figure.
fig, ax = plt.subplots(1, 2, figsize=(6, 3))
phd.viz.despine(ax)

phd.viz.titlebox(ax[0], 'van den Berg et al. 2016 | Figure S6b', color=colors['black'],
                 bgcolor='white', pad=0.05, boxsize='10%')

ax[0].vlines([0, 1, 2], [0, 0, 0], figS6b_data['mean'].values, lw=30, 
            color=colors['light_purple'], alpha=0.5)
ax[0].vlines([0, 1, 2], figS6b_data['min'].values, figS6b_data['Max'].values, lw=1, 
            color=colors['purple'])

ax[0].set_xlim([-0.5, 2.5])
ax[0].set_xticks([0, 1, 2])
ax[0].set_ylabel('percent survival')
ax[0].set_xticklabels(figS6b_data['sample'].values, rotation=0)

# Add labels and format axes
ax[1].set_xlabel('MscL channels per cell')
ax[1].set_ylabel('percent survival')
phd.viz.titlebox(ax[1], 'van den Berg et al. 2016 | Figure 5b', color=colors['black'],
                 bgcolor='white', pad=0.05, boxsize='10%')

# Add data and errorbars
_ = ax[1].errorbar(fig5b_data['x'], fig5b_data['y'], y_err * fig5b_data['y'], fmt='o', 
                color=colors['purple'], ms=4.5, lw=0.75, markeredgecolor='white',
                markeredgewidth=0.5)
_ = ax[1].hlines(fig5b_data['y'], fig5b_data['xmin'], fig5b_data['xmax'], 
                    color=colors['purple'], lw=0.75)
plt.tight_layout()
fig.text(0.01, 0.96, '(A)', fontsize=8)
fig.text(0.5, 0.96, '(B)', fontsize=8)
plt.savefig('../figs/ch9_figS13.pdf', bbox_inches='tight')
plt.savefig('../figs/ch9_figS13.png', bbox_inches='tight')
# %%
