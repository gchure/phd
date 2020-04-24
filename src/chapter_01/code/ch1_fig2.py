#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import phd.viz
import scipy.interpolate
colors, palette = phd.viz.phd_style()

# Load the three data set. 
diauxie = pd.read_csv('../../data/ch1_introduction/Monod1941_Fig1_Fig2.csv')
triauxie = pd.read_csv('../../data/ch1_introduction/Monod1947_Fig11.csv')

#%%
fig, ax = plt.subplots(1, 2, figsize=(6, 2))
phd.viz.despine(ax.ravel())
for a in ax:
    a.set_xlabel('hours')
    a.set_ylabel('optical density')
    a.spines['bottom'].set_position(('outward', 2))
    a.spines['left'].set_position(('outward', 2))
ax[0].set_ylim([-0.12, 95])
ax[1].set_ylim([-0.12, 80])
ax[0].set_xlim([-0.12, 10])
ax[1].set_xlim([-0.12, 6])


# Diauxie plots. 
sugar_colors = {'glucose':colors['blue'], 'arabinose':colors['green']}
for g, d in diauxie.groupby('secondary_sugar'):
    time = np.linspace(0, d['hours'].max()+0.1, 200)
    spline = scipy.interpolate.UnivariateSpline(d['hours'], d['optical_density'])    
    ax[0].plot(time, spline(time), '-', color=sugar_colors[g], lw=0.75, label='__nolegend__')
    ax[0].plot(d['hours'], d['optical_density'], 'o', color=sugar_colors[g],
                markeredgecolor=colors['grey'], markeredgewidth=0.5, label=f'saccharose + {g}',
            ms=5, linewidth=0.5)

# Triauxie plot
time = np.linspace(0, triauxie['time_hrs'].max() + 0.1, 200)
spline = scipy.interpolate.UnivariateSpline(triauxie['time_hrs'], triauxie['optical_density'])
spline.set_smoothing_factor(0.5)
ax[1].plot(time, spline(time), '-', lw=0.75, color=colors['purple'])
ax[1].plot(triauxie['time_hrs'], triauxie['optical_density'], 'o', color=colors['light_purple'],
       markeredgecolor='white', markeredgewidth=0.5, ms=5)

ax[0].vlines(6.2, 0, 100, color='white', lw=18, alpha=0.75)
ax[0].text(3.2, 60, 'sucrose\n+ glucose', color=colors['blue'], fontsize=7)
ax[0].text(7.4, 38, 'sucrose\n+ arabinose', color=colors['green'], fontsize=7)
ax[1].text(0.5, 20, 'glucose', color=colors['purple'], fontsize=7)
ax[1].text(2, 40, 'sorbitol', color=colors['purple'], fontsize=7)
ax[1].text(3.25, 65, 'glycerol', color=colors['purple'], fontsize=7)
# Add titles
phd.viz.titlebox(ax[0], 'Monod 1941', color=colors['black'], bgcolor='white', size=7)
phd.viz.titlebox(ax[1], 'Monod 1947', color=colors['black'], bgcolor='white', size=7)
plt.subplots_adjust(wspace=0.4)

# Add panel labels. 
fig.text(0.05, 0.84, '(A)', fontsize=8)
fig.text(0.5, 0.84, '(B)', fontsize=8)
plt.savefig('../figs/ch1_fig2.png', bbox_inches='tight', dpi=300)
plt.savefig('../figs/ch1_fig2.pdf', bbox_inches='tight')


# %%
