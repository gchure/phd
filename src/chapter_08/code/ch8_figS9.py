#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import phd.viz
import phd.thermo
import phd.stats
colors, color_list = phd.viz.phd_style()

# Load the data sets and restrict
samples = pd.read_csv('../../data/ch8_growth_si/DNA_binding_energy_samples.csv', 
                      comment='#')
samples = samples[samples['temp'] == 37]
summary = pd.read_csv('../../data/ch8_growth_si/DNA_binding_energy_summary.csv', 
                      comment='#')
summary = summary[summary['temp']==37]

# Define the bins for the histograms
ep_bins = np.linspace(-15.5, -13.5, 75)
sig_bins = np.linspace(0.2, 0.6, 75)

# %%
# Set up the figure canvas
fig, ax = plt.subplots(2,2, figsize=(4, 4), dpi=100)
phd.viz.despine(ax.ravel())

# Format the axes
ax[0, 1].set_visible(False)
ax[0, 0].yaxis.grid(False)
ax[0, 0].set_xticks([])
ax[0, 0].spines['bottom'].set_visible(False)
ax[0, 0].set_yticks([])
ax[1, 1].yaxis.grid(False)
ax[1, 1].set_yticks([])
ax[0, 0].set_xticklabels([])
ax[1, 0].set_xlim([ep_bins[0], ep_bins[-1]])
# ax[1, 0].set_xticks([-14.75, -14.25, ])
ax[1, 0].set_ylim([sig_bins[0], sig_bins[-1]])

# Set the titles of the marginals for clarity
mwc.viz.titlebox(ax[0, 0], '$\epsilon$', boxsize="15%", pad=0.03, color=colors['black'])
mwc.viz.titlebox(ax[1, 1], '$\sigma$', boxsize="15%", pad=0.03, color=colors['black'])

# Add the axis labels
ax[1, 0].set_ylabel('$\sigma$')
ax[1, 0].set_xlabel('$\epsilon$ [$k_BT$]')
ax[1, 1].set_xlabel('$\sigma$')

carb_colors = {'glucose':colors['purple'], 'glycerol':colors['green'], 
          'acetate':colors['brown']}
iter = 1
for g, d in samples.groupby('carbon'):
    # Isoalte the parameters
    ep = d[d['parameter']=='epRA']['value'].values
    sig = d[d['parameter']=='sigma']['value'].values

    # Compute the posterior distributions
    ep_hist, _ = np.histogram(ep, bins=ep_bins, density=True)
    sig_hist, _ = np.histogram(sig, bins=sig_bins, density=True)

    # Plot the distributions
    ax[0, 0].step(ep_bins[:-1], ep_hist, color=colors['grey'], lw=0.75, 
                  alpha=0.75, zorder=iter + 1)
    ax[0, 0].fill_between(ep_bins[:-1], ep_hist, color=carb_colors[g], 
                          alpha=0.75, zorder=iter, step='pre')
    
    ax[1, 1].step(sig_bins[:-1], sig_hist, color=colors['grey'], lw=0.75, 
                alpha=0.75, zorder=iter + 1)
    ax[1, 1].fill_between(sig_bins[:-1], sig_hist, color=carb_colors[g], 
                          alpha=0.75, zorder=iter, step='pre')
    ax[1, 0].plot(ep, sig, ',', color=carb_colors[g], alpha=0.15)
    
    # Add blanks for a legend
    ax[0, 1].plot([], [], '-', lw=2, color=carb_colors[g], label=f'{g}, 37° C')
    iter += 1


# Add a legend
ax[0, 1].legend(loc='center', fontsize=8)
plt.subplots_adjust(wspace=0.05, hspace=0.05)
plt.savefig('../figs/ch8_figS9.pdf', bbox_inches='tight') 
plt.savefig('../figs/ch8_figS9.png', bbox_inches='tight') 
            
# %%
