#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import phd.viz
colors, palette = phd.viz.phd_style()

# Load the three data set. 
diauxie = pd.read_csv('../../data/ch1_introduction/Monod1941_Fig1_Fig2.csv')
con_diauxie = pd.read_csv('../../data/ch1_introduction/Monod1947_Fig6.csv')
muts = pd.read_csv('../../data/ch1_introduction/Monod1946_Fig5.csv')

# %% Set up the figure canvas. 
fig, ax = plt.subplots(1, 3, figsize=(6, 2))
for a in ax:
    a.set_xlabel('hours')
    a.set_ylabel('optical density')
ax[0].set_ylim([0, 130])
ax[1].set_ylim([-0.1, 80])
ax[0].set_xlim([-0.1, 9])
ax[1].set_xlim([-0.1, 13])
ax[2].set_xlim([-0.1, 16])
ax[2].set_ylim([5, 130])


# Diauxie plots. 
sugar_colors = {'glucose':colors['purple'], 'arabinose':colors['orange']}
for g, d in diauxie.groupby('secondary_sugar'):
    ax[0].plot(d['hours'], d['optical_density'], '--o', color=sugar_colors[g],
    markeredgecolor=colors['black'], markeredgewidth=0.25, label=f'saccharose + {g}',
    ms=3)

# Add indicator of diauxic shift
ax[0].plot(6, 36, '^', color='white', markeredgecolor=colors['orange'], 
            markeredgewidth=0.5, ms=4, label='diauxic shift')

# Controlled diauxie plots
ratio_colors = {1/3: colors['blue'], 1: colors['purple'], 3: colors['orange']}
for g, d in con_diauxie.groupby(['glucose_sorbitol']):
    if g < 1:
        color = colors['blue']
    else:
        color = ratio_colors[g]
    ax[1].plot(d['hours'], d['optical_density'], '--o', color=color,
        markeredgecolor=colors['black'], markeredgewidth=0.25, 
        label=np.round(g, decimals=2), ms=3)

ax[1].plot(2.4, 15, '^', color='white', ms=4, markeredgewidth=0.5, markeredgecolor=colors['blue'],
        label='__nolegend__')
ax[1].plot(7, 25, '^', color='white', ms=4, markeredgewidth=0.5, markeredgecolor=colors['purple'],
        label='__nolegend__')
ax[1].plot(10.5, 34, '^', color='white', ms=4, markeredgewidth=0.5, markeredgecolor=colors['orange'],
        label='__nolegend__')


# Mutation Plots
mut_colors = {'L+':colors['blue'], 'L-':colors['green']}
for g, d in muts.groupby(['strain']):
    if g == 'L+':
        label = 'lactose positive'
    else:
        label = 'lactose negative'

    ax[2].plot(d['hours'], d['optical_density'], '--o', color=mut_colors[g],
            markeredgecolor=colors['black'], markeredgewidth=0.25, ms=3,
            label=label)

# Add diauxie labels
ax[2].plot(6, 40, '^', color='white',markeredgecolor=colors['blue'], markeredgewidth=0.5,
            label='__nolegend__', ms=4)
ax[2].plot(11, 38, '^', color='white',markeredgecolor=colors['green'], markeredgewidth=0.5,
            label='__nolegend__', ms=4)
legs = [0, 0, 0]

legs = [0, 0, 0]
legs[0] = ax[0].legend(loc='upper left', fontsize=6)
legs[1] = ax[1].legend(loc='upper left', title='glucose / xylose', fontsize=6, ncol=3,
                    columnspacing=0.5)
legs[2] = ax[2].legend(loc='upper left', fontsize=6)

for l in legs:
    l.get_title().set_fontsize(6)

# Add titles
phd.viz.titlebox(ax[0], 'Monod 1941', color=colors['black'], bgcolor='white', size=6)
phd.viz.titlebox(ax[1], 'Monod 1947', color=colors['black'], bgcolor='white', size=6)
phd.viz.titlebox(ax[2], 'Monod & Audureau 1946', color=colors['black'], bgcolor='white', size=6)
plt.subplots_adjust(wspace=0.35)

# Add panel labels. 
fig.text(0.05, 0.84, '(A)', fontsize=8)
fig.text(0.35, 0.84, '(B)', fontsize=8)
fig.text(0.622, 0.84, '(C)', fontsize=8)
plt.savefig('../figs/ch1_fig2.png', bbox_inches='tight', dpi=300)
plt.savefig('../figs/ch1_fig2.pdf', bbox_inches='tight')


# %%
