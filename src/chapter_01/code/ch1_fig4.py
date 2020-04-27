#%% 
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import phd.viz
import phd.thermo
constants = phd.thermo.load_constants()
colors, palette = phd.viz.phd_style()

# %%
# Load the experimental data. 
data = pd.read_csv('../../data/ch2_induction/RazoMejia_2018.csv', comment='#')
data = data[data['repressors'] > 0]
data['repressors'] *= 2
data_summarized = data.groupby(['operator', 'IPTG_uM', 'repressors'])['fold_change_A'].mean().reset_index()

#%%
c_range = np.logspace(-2, 4, 200)
rep_colors = {22:colors['red'], 60:colors['brown'], 124:colors['green'],
              260:colors['orange'], 1220:colors['purple'], 1740:colors['blue']}

fig, ax = plt.subplots(3, 1, sharey=True, figsize=(3, 4))
phd.viz.despine(ax.ravel())
for i, a in enumerate(ax):
    a.set_xscale('log')
    a.set_xticks([1E-1, 1E1, 1E3])
    a.set_xlim([1E-2, 1E4])
    if i < 2:
        a.set_xticks([])
        a.spines['bottom'].set_visible(False)
    a.set_ylabel('fold-change', fontsize=6)
ax[-1].set_xlabel('inducer [µM]', fontsize=6)
# Define the operator axes
op_ax = {'O1':ax[0], 'O2':ax[1], 'O3':ax[2]}

for o, a in op_ax.items():
    phd.viz.titlebox(a, f'operator {o}', size=6, color=colors['black'],
                    boxsize='15%')


for g, d in data_summarized.groupby(['operator', 'repressors']):
    theo = phd.thermo.SimpleRepression(R=g[1], ep_r=constants[g[0]],
                                       ka=constants['Ka'], ki=constants['Ki'],
                                       ep_ai=constants['ep_AI'], 
                                       effector_conc=c_range).fold_change()
    if (g[0] == 'O2') & (g[1] == 260):
        fill = 'white'
        edge = rep_colors[g[1]]
        zorder = 1000
    else:
        fill = rep_colors[g[1]]
        edge = colors['grey']
        zorder = 11
    op_ax[g[0]].plot(c_range, theo, '-', lw=0.75, color=rep_colors[g[1]], 
                    label=int(g[1]), zorder=10)
    op_ax[g[0]].plot(d['IPTG_uM'], d['fold_change_A'], 'o', markerfacecolor=fill,
                    markeredgecolor=edge, markeredgewidth=0.5, ms=2.5, label='__nolegend__', 
                    zorder=zorder)

leg.get_title().set_fontsize(6)
plt.savefig('../figs/ch1_fig4_induction_plot.svg')
# %%

