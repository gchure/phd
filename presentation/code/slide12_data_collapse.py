#%%
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import phd.thermo
import phd.viz
colors, palette = phd.viz.phd_style()

# Load the compiled data
data = pd.read_csv('../../src/data/other/Garcia2011_Brewster2014_RazoMejia2018_Chure2019.csv')
data = data[data['author']=='razo-mejia']

# %%

#define the modifiers
DATA = True

fig, ax = plt.subplots(1, 1, figsize=(3, 3))
phd.viz.despine(ax)
ax.set_xlabel('free energy [$k_BT$]')
ax.set_ylabel('fold-change')
ax.set_xlim([-10, 10])
ax.set_ylim([-0.05, 1.18])

# Define the theoretical curves. 
bohr_range = np.linspace(-10, 10, 200)
scaling_fn = (1 + np.exp(-bohr_range))**-1

ax.plot(bohr_range, scaling_fn, 'k-', label='scaling function', lw=0.75)
rep_colors = {22: colors['red'], 
              60: colors['brown'],  
              124: colors['green'], 
              260: colors['orange'], 
              1220: colors['purple'], 
              1740: colors['blue']} 
glyphs = {'O1':'v', 'O2':'o', 'O3':'X'}

if DATA:
    for g, d in data.groupby(['repressors', 'operator']):
        ax.errorbar(d['bohr_parameter'], d['mean'], d['sem'], ms=4.5, color=rep_colors[g[0]],
                    linestyle='none', marker=glyphs[g[1]], label='__nolegend__', 
                    markeredgecolor=colors['grey'], markeredgewidth=0.5)

    for r, c in rep_colors.items():
        ax.plot([], [], '-', color=c, label=f'{int(r)} repressors', lw=2)

    for g, m in glyphs.items():
        ax.plot([], [], 'k', ms=4.5, markeredgecolor=colors['grey'], markeredgewidth=0.75, 
                linestyle='none', marker=m, label=f'operator {g}')

    ax.legend(fontsize=6, handlelength=1.5)

    name = 'data'
else:
    name='nodata'
plt.savefig(f'../figs/slide12_induction_collapse_{name}.pdf', bbox_inches='tight')
# %%
