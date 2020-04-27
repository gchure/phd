#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.interpolate
colors, palette = phd.viz.phd_style()

# Load  the data
data = pd.read_csv('../../data/ch1_introduction/Monod1946_Fig5.csv')


# %%
fig, ax = plt.subplots(1, 1, figsize=(3, 2))
phd.viz.despine(ax)
ax.set_xlabel('hours')
ax.set_ylabel('optical density')

strain_colors = {'L+':colors['red'], 'L-':colors['black']}
labels = {'L+':'L$+$ (can digest lactose)', 'L-':'L$-$ (cannot digest lactose)'}
for g, d in data.groupby(['strain']):
    # Set up the time range
    time_range = np.linspace(d['hours'].min(), d['hours'].max(), 200)
    spline = scipy.interpolate.UnivariateSpline(d['hours'], d['optical_density'])   
    # spline.set_smoothing_factor(0.001)
    ax.plot(time_range, spline(time_range), '-', color=strain_colors[g], lw=0.75,
            label='__nolegend__')
    ax.plot(d['hours'], d['optical_density'], 'o', color=strain_colors[g],
            markeredgecolor=colors['grey'], markersize=4,
            label=labels[g], markeredgewidth=0.5)
phd.viz.titlebox(ax, 'Monod 1947', size=6, pad=0.05, boxsize='12%', color=colors['black'])
ax.legend()
plt.savefig('../figs/ch1_fig6.pdf', bbox_inches='tight')
plt.savefig('../figs/ch1_fig6.png', bbox_inches='tight')

# %%
