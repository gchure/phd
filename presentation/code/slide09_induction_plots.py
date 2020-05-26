#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle
import phd.viz
import phd.thermo
colors, palette = phd.viz.phd_style()
constants = phd.thermo.load_constants()

# Define the identifiers
SLIDE = 9
FIT_DATA = True
FIT_CURVE = True
ALL_PREDICTIONS = True
ALL_DATA = True



# Load the induction paper dataset
data = pd.read_csv('../../src/data/ch2_induction/RazoMejia_2018.csv', comment='#')
data['repressors'] *= 2
data = data[data['repressors'] > 0]
data = data.groupby(['repressors', 'operator', 
                     'IPTG_uM'])['fold_change_A'].agg(('mean', 'sem')).reset_index()

# Load the sampling information
with open('../../src/data/ch2_induction/mcmc/main_text_KaKi.pkl', 'rb') as file:
    unpickler = pickle.Unpickler(file)
    gauss_flatchain = unpickler.load()

ka_fc = np.exp(-gauss_flatchain[:, 0])[::100]
ki_fc = np.exp(-gauss_flatchain[:, 1])[::100]

# Set up the figure canvas. 
fig, ax = plt.subplots(1, 3, figsize=(7, 2), sharex=True, sharey=True)
for a in ax:
    phd.viz.despine(a)
    a.set_xscale('log')
    a.set_xlim([1E-2, 1E4])
    a.set_ylim([-0.05, 1.15])

# Define the colors
rep_colors = {22: colors['red'], 
              60: colors['brown'],  
              124: colors['green'], 
              260: colors['orange'], 
              1220: colors['purple'], 
              1740: colors['blue']}

# Define the axes
axes = {'O1':ax[0], 'O2':ax[1], 'O3':ax[2]}

# Add labels as necessary
for i, a in enumerate(ax):
    phd.viz.titlebox(a, f'operator O{i+1}', color=colors['black'], size=6)
    a.set_xlabel('IPTG [µM]')
ax[0].set_ylabel('fold-change')
# plot the data and curves if desired
for g, d in data.groupby(['operator', 'repressors']):
    if (g[0] == 'O2') & (g[1] == 260):
        face = colors['grey']
        edge = rep_colors[g[1]]
        _ = axes[g[0]].errorbar(d['IPTG_uM'], d['mean'], d['sem'], fmt='o',
                        linestyle='none', color=rep_colors[g[1]], 
                        markeredgecolor=edge, markerfacecolor=face,
                        markeredgewidth=0.75, ms=3, label=int(g[1]))
    else:
        face = rep_colors[g[1]]
        edge = colors['grey']

        if ALL_DATA:
            _ = axes[g[0]].errorbar(d['IPTG_uM'], d['mean'], d['sem'], fmt='o',
                        linestyle='none', color=rep_colors[g[1]], 
                        markeredgecolor=edge, markerfacecolor=face,
                        markeredgewidth=0.75, ms=4, label='__nolegend__')


# Plot curves as necessary
c_range = np.logspace(-2, 4, 200)
if (FIT_CURVE != False) | (ALL_PREDICTIONS != False):
    for g, d in data.groupby(['operator', 'repressors']):
        cred_region = np.zeros((2, len(c_range)))
        for i, c in enumerate(c_range):
            arch = phd.thermo.SimpleRepression(g[1], ep_r=constants[g[0]],
                                               ka=ka_fc, ki=ki_fc,
                                               ep_ai=constants['ep_AI'],
                                               effector_conc=c).fold_change()
            cred_region[:, i] = phd.stats.compute_hpd(arch, 0.95)

        if (g[0] == 'O2') & (g[1] == 260) & (FIT_CURVE != False):
            axes[g[0]].fill_between(c_range, cred_region[0, :], cred_region[1, :],
                            color=rep_colors[g[1]], alpha=0.5, label=int(g[1]))
        elif ALL_PREDICTIONS == True:
            axes[g[0]].fill_between(c_range, cred_region[0, :], cred_region[1, :],
                            color=rep_colors[g[1]], alpha=0.5, label=int(g[1]))

leg = ax[1].legend(title='repressors')
leg.get_title().set_fontsize(6)

name = ''
if (FIT_DATA == True):
    name += 'fitdata_'
if (FIT_CURVE == True):
    name += 'fitcurve_'
if (ALL_PREDICTIONS==True):
    name += 'allcurve_'
if (ALL_DATA == True):
    name += 'alldata_'
plt.savefig(f'../figs/slide_{SLIDE}_induction_{name}.pdf', bbox_inches='tight')
# %%

# %%


