#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import phd.viz
import phd.bayes
import scipy.stats
import scipy.special
colors, palette = phd.viz.phd_style()

# %%
# Load data, restrict, and define new constants
fluct_data = pd.read_csv('../../data/ch8_growth_si/analyzed_fluctuations.csv', comment='#')
fluct_data = fluct_data[fluct_data['date']==20181002]
alpha_range = np.linspace(550, 950, 500)

#%% 
# Instantiate the figure
fig, ax = plt.subplots(1, 3, figsize=(6, 2))
phd.viz.despine(ax)

# Format the axes
ax[0].set_xlabel('$\log_{10}$ cell intensity')
ax[0].set_ylabel('number of observations')
ax[1].set_xlabel('$I_1 / (I_1 + I_2)$')
ax[1].set_ylabel('number of observations')
ax[2].set_xlabel('cell volume [fL]')
ax[2].set_ylabel('fractional intensity')
ax[2].set_xlim(0, 3.5)
avg_partitioning  = np.mean(fluct_data['I_1'] / fluct_data['summed'])

# Plot the fluctuations
ax[0].hist(np.log10(fluct_data['I_1']), bins=50, color=colors['purple'],
            edgecolor='white', label='$I_1$')
ax[0].hist(np.log10(fluct_data['I_2']), bins=50, color=colors['orange'],
            edgecolor='white', alpha=0.5, label='$I_2$')

ax[1].hist(fluct_data['I_1'] / fluct_data['summed'], bins=50, color='grey', 
             edgecolor=colors['black'], alpha=0.5, lw=0.25)
ax[1].vlines(0.5, 0, ax[1].get_ylim()[1], color=colors['orange'], label='even partitioning\n(0.50)')
ax[1].vlines(avg_partitioning, 0, ax[1].get_ylim()[1], 
        color=colors['dark_purple'], label=f'average value\n({np.round(avg_partitioning, decimals=3)})')
        
ax[2].plot(fluct_data['volume_1_birth'], fluct_data['I_1'] / fluct_data['summed'],
             'k.', ms=0.75, alpha=0.75)
ax[2].plot(fluct_data['volume_2_birth'], fluct_data['I_2'] / fluct_data['summed'],
              'k.', ms=0.75,  alpha=0.75)

ax[0].legend()
ax[1].legend(handlelength=0.5, fontsize=5)
# ax[2].legend(handlelength=0.5, fontsize=5)
plt.tight_layout()
fig.text(0, 0.95, '(A)', fontsize=8)
fig.text(0.34, 0.95, '(B)', fontsize=8)
fig.text(0.66, 0.95, '(C)', fontsize=8)
plt.savefig('../figs/ch8_figS4.pdf', bbox_inches='tight')
plt.savefig('../figs/ch8_figS4.png', bbox_inches='tight')



# %%
