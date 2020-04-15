#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import phd.viz
colors, palette = phd.viz.phd_style()


#%% Load the fluctuations and the snap infos
flucts = pd.read_csv('../../data/ch4_growth/analyzed_fluctuations.csv', comment='#')
snaps = pd.read_csv('../../data/ch4_growth/analyzed_foldchange.csv', comment="#")
snaps = snaps[snaps['atc_ngml']==10]
condition_colors = {'glucose':'purple', 'glycerol':'green', 'acetate':'brown',
                    42:'red', 32:'blue'}

MIN_THRESH = 2.5
MAX_THRESH = 3.5
#%%
bins = np.linspace(0, 10, 50)
fig, ax = plt.subplots(2, 5, figsize=(6, 2.8))
phd.viz.despine(ax.ravel())
for i in range(5):
    ax[1, i].set_xlabel('length [µm]')

ax[0, 0].set_ylabel('$\propto$ probability')
ax[1, 0].set_ylabel('cumulative distribution')
iter = 0
for g, d in snaps.groupby(['carbon', 'temp']):
    # Get the corresponding fluctuations
    _flucts = flucts[(flucts['carbon']==g[0]) & (flucts['temp']==g[1])]

    if g[1] == 37:
        tc = 'white'
        fc = colors[f'{condition_colors[g[0]]}']
    else:
        tc = 'white'
        fc = colors[f'{condition_colors[g[1]]}']
    phd.viz.titlebox(ax[0,iter], text=f'{g[0]}, {g[1]}° C',
                     size=6, color=colors['black'], bgcolor='white', boxsize="15%")

    # plot the histogram of the data lengths. 
    _ = ax[0, iter].hist(d['length_um'], bins=bins, histtype='stepfilled',
                        edgecolor=tc, facecolor=fc, alpha=0.75, density=True)
    _ = ax[0, iter].hist(_flucts[['length_1_birth', 'length_2_birth']].values.flatten(), 
                        bins=bins, histtype='stepfilled',
                        edgecolor=colors['black'], facecolor=colors['black'], 
                        alpha=0.25, density=True)
    # Plot the ecdfs. 
    lengths = _flucts[['length_1_birth', 'length_2_birth']].values.flatten()
    xb, yb = np.zeros(len(lengths) + 2), np.zeros(len(lengths) + 2)
    xb[-1] = 100
    yb[-1] = 1
    _x, _y = np.sort(lengths), np.arange(1, len(lengths)+1, 1) / len(lengths)
    xb[1:-1] = _x
    yb[1:-1] = _y


    # Plot the ecdf. 
    x, y = np.zeros(len(d) + 2), np.zeros(len(d) + 2)
    x[-1] = 100
    y[-1] = 1
    _x, _y = np.sort(d['length_um']), np.arange(1, len(d)+1, 1) / len(d)
    x[1:-1] = _x
    y[1:-1] = _y

    _ = ax[1, iter].step(x, y, color=fc)
    _ = ax[1, iter].step(xb, yb, color=colors['black'])
    for k in range(2):
        _ = ax[k, iter].vlines(MIN_THRESH, 0, 2, color=colors['red'], lw=0.75)
        _ = ax[k, iter].vlines(MAX_THRESH, 0, 2, color=colors['red'], lw=0.75)
    iter += 1

for a in ax.ravel():
    a.set_xlim([0, 6])
for i in range(5):
    ax[0, i].spines['bottom'].set_visible(False)
    ax[0, i].set_xticks([])
    ax[0, i].set_xticklabels([])
    ax[0, i].set_ylim([0, 1.1])
    ax[1, i].set_ylim([0, 1])
    if i > 0:
        for j in range(2):
            ax[j, i].set_yticklabels([])
            ax[j,i].set_yticks([])
            ax[j, i].spines['left'].set_visible(False)

plt.subplots_adjust(wspace=0.1, hspace=0.1)
plt.savefig('../figs/ch8_figS7.pdf', bbox_inches='tight')
plt.savefig('../figs/ch8_figS7.png', bbox_inches='tight')


# %%
