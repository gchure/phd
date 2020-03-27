#%% 
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import phd.viz
import phd.thermo
constants = phd.thermo.load_constants()
colors, palette = phd.viz.phd_style()

# Load the various datasets
newgods = pd.read_csv('../../data/ch2_induction/RazoMejia_2018.csv', comment='#')
newgods = newgods[(newgods['IPTG_uM']==0) & (newgods['repressors'] > 0)]
newgods['repressors'] *= 2
newgods = newgods.groupby(['operator', 
            'repressors'])['fold_change_A'].agg(('mean', 'sem')).reset_index()
oldgods = pd.read_csv('../../data/other/Garcia2011_Brewster2014.csv', comment='#')
oldgods = oldgods[oldgods['operator']!='Oid']

# Define the colors. 
viridis = sns.color_palette('viridis', n_colors=4)
op_colors = {'O1':viridis[0], 'O2':viridis[1], 'O3':viridis[2]}

# Set up constants for plotting
rep_range = np.logspace(0, 4, 200)
r, ep = np.meshgrid(rep_range, np.array([constants['O1'], constants['O2'], constants['O3']]))

# Compute the fold-change curves
fc = phd.thermo.SimpleRepression(R=r, ep_r=ep, ka=constants['Ka'], 
                                ki=constants['Ki'], ep_ai=constants['ep_AI'],
                                effector_conc=0).fold_change()

# Set up the figure canvas.
fig, ax = plt.subplots(1, 1, figsize=(4, 3))
phd.viz.despine(ax)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('repressors per cell', fontsize=8)
ax.set_ylabel('fold-change', fontsize=8)
ax.set_xlim([1, 5E3])

# plot the theory curves
iter = 0
for o, c in op_colors.items():
    ax.plot(rep_range, fc[iter], color=c, label=o)
    iter += 1

# Plot the data.
for g, d in oldgods.groupby(['author', 'operator']):
    if g[0] == 'brewster':
        glyph = '^'
    elif g[0] == 'garcia':
        glyph = 'X'
    ax.plot(d['repressor'], d['fold_change'], color=op_colors[g[1]], 
            marker=glyph, markeredgecolor='white', markersize=5,
            label='__nolegend__', linestyle='none', markeredgewidth=0.75)

for g, d in newgods.groupby(['operator']):
    ax.errorbar(d['repressors'], d['mean'], d['sem'], fmt='o', ms=4.5, linestyle='none',
    markeredgecolor='white', markeredgewidth=0.75, label='__nolegend__', color=op_colors[g])
    
# Add the legend info
ax.plot([], [], 'kX', label='Garcia and Phillips, 2011\n(Miller Assay)',
        markeredgecolor='white', markersize=5, markeredgewidth=0.75)
ax.plot([], [], 'k^', label='Brewster et al., 2014\n(Microscopy)',
        markeredgecolor='white', markersize=5, markeredgewidth=0.75)
ax.plot([], [], 'ko', label='flow cytometry', markeredgecolor='white',
        markeredgewidth=0.75, markersize=4.5)

ax.legend(fontsize=6)
plt.savefig('../figs/ch6_figS9.pdf', bbox_inches='tight')
plt.savefig('../figs/ch6_figS9.png', bbox_inches='tight')
# %%



# %%
