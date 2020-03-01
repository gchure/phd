#%% 
import glob
import numpy as np
import pandas as pd
import scipy.stats
import phd.viz
import phd.stats
import phd.flow
import fcsparser
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn as sns
from matplotlib import ticker
colors, palette = phd.viz.phd_style()

# Purpose is to generate figures of distributions from flow cytometry for
# a flowchart diagram.
flow_glob = glob.glob('../../data/ch2_induction/example_flow/201608*O2*RBS1027*.fcs')
meta, cell_cloud = fcsparser.parse(flow_glob[3])

# Fit a two-dimensional gaussian to the data.
mu, cov = phd.flow.fit_2D_gaussian(cell_cloud, log=True)

# Compute the statistic for each of the pair of log scattering data
interval_array = phd.flow.gauss_interval(cell_cloud, mu, cov, log=True)
alpha = 0.40

#%%
# Find which data points fall inside the interval
idx = interval_array <= scipy.stats.chi2.ppf(alpha, 2)
disc_idx = interval_array >= scipy.stats.chi2.ppf(alpha, 2)

# Select those data which lay within the 40th percentile.
selected_cells = cell_cloud[idx]
discarded_cells = cell_cloud[disc_idx]

# Now just generate the plot.
fig, ax = plt.subplots(1, 1, figsize=(1.5, 1), dpi=150)
phd.viz.despine(ax)
plt.plot(discarded_cells['FSC-A'], discarded_cells['SSC-A'], marker=',', color=colors['black'],
                   rasterized=True, label='discarded cells', linestyle='none')
plt.plot(selected_cells['FSC-A'], selected_cells['SSC-A'], marker=',', color=colors['purple'],
               rasterized=True, label='selected cells', linestyle='none')


# Fix formatting and restrict bounds.
plt.xlabel('forward scatter [a.u.]', fontsize=8)
plt.ylabel('side scatter [a.u.]', fontsize=8)
plt.xscale('log')
plt.yscale('log')
plt.tick_params(labelsize=6)

# Restrict bounds for aesthetic reasons.
plt.xlim([1E3, 1E5])
plt.ylim([1E3, 3E5])

# Save the figure.
plt.savefig('../figs/fig4_flow_cloud.pdf', bbox_inches='tight')

#%%
fig, ax = plt.subplots(2,1, figsize=(3,2), sharex=True)
phd.viz.despine(ax)
# Now generate the example distributions.
colors_RBS1027 = sns.color_palette('Oranges_r', n_colors=8)
colors_delta = sns.color_palette('Purples_r', n_colors=8)
IPTG_range = [0, 25, 50, 100, 500, 5000]
ax[0].plot([], [], 'v', markersize=5, markeredgecolor='white',
            markerfacecolor=colors_RBS1027[-1], label='mean', linestyle='none')
ax[1].plot([], [], 'v', markersize=5, markeredgecolor='white', markeredgewidth=0.5,
            markerfacecolor=colors_delta[-1], label='mean', linestyle='none')
leg = ax[0].legend(title=r'repressors / cell = 260')
leg.get_title().set_fontsize(6)
leg = ax[1].legend(title=r'repressors / cell = 0')
leg.get_title().set_fontsize(6)

for i, val in enumerate(IPTG_range):

    glob_RBS1027 = glob.glob(f'../../data/ch2_induction/example_flow/201608*O2*RBS1027*_{val}uM*IPTG.fcs')
    _, data_RBS1027  = fcsparser.parse(glob_RBS1027[0])
    glob_delta = glob.glob(f'../../data/ch2_induction/example_flow/201608*O2*delta*_{val}uM*IPTG.fcs')
    _, data_delta = fcsparser.parse(glob_delta[0])

    # Fit a two-dimensional gaussian to the data.
    gate_RBS1027= phd.flow.gaussian_gate(data_RBS1027, alpha)
    gate_delta = phd.flow.gaussian_gate(data_delta, alpha)
    ax[0].hist(gate_RBS1027['FITC-A'], color=colors_RBS1027[i],alpha=0.5, bins=100,
            histtype='stepfilled', density=True, edgecolor=colors['dark_red'])
    mean_RBS1027 = np.mean(gate_RBS1027['FITC-A'])
    ax[0].plot(mean_RBS1027, 1.3E-4, 'v', markeredgecolor=colors['dark_red'], markeredgewidth=0.5,
             markerfacecolor=colors_RBS1027[i], markersize=6)
    ax[1].hist(gate_delta['FITC-A'], color=colors_delta[i], alpha=0.5, bins=100,
            histtype='stepfilled', density=True, edgecolor=colors['dark_purple'])
    mean_delta = np.mean(gate_delta['FITC-A'])
    ax[1].plot(mean_delta, 7.5E-5, 'v', markeredgecolor=colors['dark_purple'], markeredgewidth=0.5,
             markerfacecolor=colors_delta[i], markersize=6)

ax[0].yaxis.get_major_formatter().set_powerlimits((0, -1))
ax[1].yaxis.get_major_formatter().set_powerlimits((0, -1))
ax[0].xaxis.get_major_formatter().set_powerlimits((0, -1))
ax[1].xaxis.get_major_formatter().set_powerlimits((0, -1))

for a in ax:
        a.set_xlim([-0.1E5, 0.8E5])
fig.text(0, 0.45, 'frequency', fontsize=8, rotation='vertical')
ax[1].set_xlabel('total cell intensity (a.u.)', fontsize=8)
plt.savefig('../figs/fig4_flow_distributions.svg')
plt.show()


#%%

# Now plot a full titration from this set.
data = pd.read_csv('../../data/ch2_induction/RazoMejia_2018.csv', comment='#')

# Isolate to the proper date and run number 
data = data[(data['date']==20160813) & (data['operator']=='O2') & 
            (data['repressors']==130) & (data['username']=='mrazomej')]

fig, ax = plt.subplots(1, 1, figsize=(2.5, 1.75), dpi=150)
phd.viz.despine(ax)
plt.plot([], [], 'o', markersize=4, markeredgecolor=colors['orange'], markerfacecolor='white',
        label='experimental data', markeredgewidth=0.5)
leg = plt.legend(loc='upper left', title="repressors / cell = 260\n " + r"$\Delta\varepsilon_{RA} = -13.9\,k_BT$")
leg.get_title().set_fontsize(6)
plt.plot(data['IPTG_uM'], data['fold_change_A'], marker='o', linestyle=':', color=colors['orange'],
          markersize=4, markeredgecolor=colors['orange'], markeredgewidth=0.5, markerfacecolor='white',
          linewidth=0.5)
plt.xscale('log')
plt.ylabel('fold-change', fontsize=6)
plt.xlabel('IPTG [µM]', fontsize=6)
plt.xlim([1E-2, 1E4])
plt.ylim([-0.01, 1.1])
plt.tick_params(labelsize=6)
plt.savefig('../figs/fig4_titration.svg')

# %%
#