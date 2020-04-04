# -*- coding: utf-8 -*-
# %%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import phd.thermo
import phd.stats
import seaborn as sns
import phd.viz
colors, palette = phd.viz.phd_style()

# Load the data sets
fug_data = pd.read_csv('../../data/ch7_mutants_si/fugacity_data.csv')
ssr_data = pd.read_csv('../../data/ch7_mutants_si/RazoMejia2018_Garcia2011_Brewster2014_tidy.csv')
ssr_data = ssr_data[ssr_data['operator'] != 'Oid']
fug_data = fug_data[fug_data['operator'] != 'Oid']

# Load the summary information
razo_samples = pd.read_csv('../../data/ch7_mutants_si/razomejia_epai_summary.csv')
matt_samples = pd.read_csv('../../data/ch7_mutants_si/matthews_epai_summary.csv')
daber_samples = pd.read_csv('../../data/ch7_mutants_si/daber_epai_summary.csv')

# Summarize the ssr data
ssr_summ = ssr_data.groupby(['operator', 'repressors', 'IPTGuM', 'author'
                            ])['fold_change'].agg(('mean', 'sem')).reset_index()


#%% 
# Instantiate the figure canvas
fig, ax = plt.subplots(2, 3, figsize=(7, 5))
phd.viz.despine(ax.ravel())
# Formatting for induction profiles
ops = ['O1', 'O2', 'O3']
for i in range(3):
    ax[0, i].set_xlabel('IPTG [µM]')
    ax[0, i].set_ylabel('fold-change')
    ax[0, i].set_xscale('symlog', linthreshx=1E-2)
    ax[0, i].set_xticks([0, 1E-1, 1E1, 1E3])
    ax[0, i].set_ylim([-0.0, 1.1])
    ax[0, i].set_xlim([-0.001, 1E4])
    phd.viz.titlebox(ax[0, i], f'operator {ops[i]}', bgcolor='white', 
                    color=colors['black'], pad=0.05, boxsize='12%')

for i in range(2):
    ax[1, i+1].set_xlabel('repressors per cell')
    ax[1, i+1].set_ylabel('fold-change')
    ax[1, i+1].set_yscale('log')
    ax[1, i+1].set_xscale('log')


# Blank one axis out for a legend
ax[1, 0].axis(False)
ax[1, 1].set_xlim([1, 1E3])

# Define styling and axis mapping
rep_colors = {22: colors['red'], 
              60: colors['brown'],  
              124: colors['green'], 
              260: colors['orange'], 
              1220: colors['purple'], 
              1740: colors['blue']}

_viridis = sns.color_palette('viridis', n_colors=5)
op_colors = {'O1':_viridis[0], 'O2':_viridis[1], 'O3':_viridis[2]}
axes = {'O1':ax[0,0], 'O2':ax[0, 1], 'O3':ax[0, 2]}
deep = sns.color_palette('deep', n_colors=3)
N_colors = {64:deep[0], 52:deep[1]}
op_mark = {'O1': 'o', 'Oid': 'o'}

# Add legend information
ax[1, 0].plot([], [], 'k-', label=r'$\Delta\varepsilon_{AI} = 4.5\, k_BT$' + '\n Razo-Mejia\n et al. 2018')
ax[1, 0].plot([], [], 'k--', label=r'$\Delta\varepsilon_{AI} = -1.75\, k_BT$' + "\n Daber \n et al. 2011")
ax[1, 0].plot([], [], 'k:', label=r'$\Delta\varepsilon_{AI} = 0.35\, k_BT$' + "\n O'Gorman\n et al. 1980")
ax[1, 0].legend(loc='center', fontsize=6)

for r, c in rep_colors.items():
    ax[0, 0].plot([], [], linestyle='none', label=int(r), marker='o', color=c, ms=3, 
                 markeredgecolor='white', markeredgewidth=0.5)
leg = ax[0, 0].legend(loc='upper left', title='rep. / cell', fontsize=6)
leg.get_title().set_fontsize(6)

for o, c in op_colors.items():
    ax[1, 1].plot([], [], 'o', label=o, color=c, ms=3, markeredgecolor='w', markeredgewidth=0.5)
ax[1, 1].legend(loc='lower left', fontsize=7)

for n, c in N_colors.items():
    ax[1, 2].plot([], [], 'o', color=c, ms=3, label=int(n), markeredgecolor='w', markeredgewidth=0.5)
leg = ax[1, 2].legend(title='num. promoters', fontsize=6)
leg.get_title().set_fontsize(6)

# Plot the summarized data
for g, d in ssr_summ.groupby(['author', 'operator', 'repressors']):
    if g[0] == 'razo-mejia':
        a = axes[g[1]]
        a.errorbar(d['IPTGuM'], d['mean'], d['sem'], capsize=1, ms=4, fmt='o',
               label=int(g[2]), linestyle='none', color=rep_colors[int(g[2])],
               lw=0.75, markeredgecolor='w', markeredgewidth=0.5)

    else:
       ax[1, 1].plot(d['repressors'], d['mean'], color=op_colors[g[1]], marker='o',
                    ms=4, markeredgecolor='w', markeredgewidth=0.5)

# Plot the fugacity data
for g, d in fug_data.groupby(['N', 'operator']):
    ax[1, 2].plot(d['repressor'], d['fold_change'], marker=op_mark[g[1]], 
                  color=N_colors[g[0]], linestyle='none', ms=4,
                  markeredgecolor='w', markeredgewidth=0.5)


# Plot the theory curves
c_range = np.logspace(-2, 4, 200)
c_range[0] = 0
rep_range = np.logspace(0, 3.5, 200)

# Parse out the parameter values
linestyles = ['-', '--', ':']
linewidths = [0.75, 0.75, 0.75]
for i, stats in enumerate([razo_samples, daber_samples, matt_samples]):
    # Retrieve the parameters
    ka = stats[stats['parameter']=='ka']['median'].values[0]
    ki = stats[stats['parameter']=='ki']['median'].values[0]
    ep_ai = stats['ep_ai'].unique()

    for o in ['O1', 'O2','O3']:
        ep_r = stats[(stats['parameter']=='ep_r') & (stats['operator']==o)]['median'].values[0]
        for r, c in rep_colors.items(): 
            fc = phd.thermo.SimpleRepression(r, ep_r, ka=ka, ki=ki, ep_ai=ep_ai,
                                             effector_conc=c_range).fold_change()
            axes[o].plot(c_range, fc, linestyle=linestyles[i], color=c, lw=linewidths[i])

        leak = phd.thermo.SimpleRepression(rep_range, ep_r, ka=ka, ki=ki, ep_ai=ep_ai,
                                           effector_conc=0).fold_change()
        ax[1,1].plot(rep_range, leak, linestyle=linestyles[i], linewidth=linewidths[i],
                    color=op_colors[o])
        if o == 'O1':
            o1 = ep_r 
    for n in [64, 52]:
        x = np.exp(-o1)
        pact = 1 / (1 + np.exp(-ep_ai))
        r = pact * rep_range
        B = r * (1 + x) - n * x - 4.6E6
        A = x * (r - n - 4.6E6)
        numer = -B - np.sqrt(B * B - 4 * A * r)
        denom = 2 * A
        _fug = 1 / (1 + (numer/denom) * x)
        ax[1, 2].plot(rep_range, _fug, linestyle=linestyles[i], linewidth=linewidths[i],
                    color=N_colors[n])


ax[1, 2].set_xlim([1, 800])
ax[1, 2].set_ylim([1E-3, 3])
plt.subplots_adjust(hspace=0.3, wspace=0.38)
fig.text(0.05, 0.9, '(A)', fontsize=8)
fig.text(0.35, 0.9, '(B)', fontsize=8)
fig.text(0.63, 0.9, '(C)', fontsize=8)
fig.text(0.34, 0.47, '(D)', fontsize=8)
fig.text(0.63, 0.47, '(E)', fontsize=8)
plt.savefig('../figs/ch7_figS20.pdf', bbox_inches='tight')
plt.savefig('../figs/ch7_figS20.png', bbox_inches='tight')
#%%
