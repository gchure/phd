
#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import phd.viz
import phd.stats
colors, palette = phd.viz.phd_style()

# Define the repressor copy number range. 
rep_range = np.logspace(0, 4, 200)

# Define the target points for the explanation
example_reps = np.array([25, 50, 100, 200, 400])
example_facecolors = {25:colors['purple'], 50: colors['orange'],
                    100:'white', 200:colors['red'], 400:colors['blue']}

ref_rep = 100

# Define thermodynamic constants
ep_ai = 4.5 # in kT
ep_RA = -13.9 # in kT
Nns = 4.6E6 # in bp
sigma = 0.2

# Define the fold-change and the deltaF theory curves
pact = (1 + np.exp(-ep_ai))**-1
fc = (1 + pact * (rep_range / Nns) * np.exp(-ep_RA))**-1 
fc_min = (1 + pact * (rep_range / Nns) * np.exp(-(ep_RA + sigma)))**-1 
fc_max = (1 + pact * (rep_range / Nns) * np.exp(-(ep_RA - sigma)))**-1 
example_fc = (1 + pact * (example_reps/ Nns) * np.exp(-(ep_RA - sigma)))**-1 
example_delF = -np.log(example_reps / ref_rep)
F_ref = -np.log(pact) - np.log(ref_rep/Nns) + ep_RA
deltaF = -np.log(rep_range / ref_rep) # in kT

# Load the actual data. 
fc_data = pd.read_csv('../../data/ch4_growth/analyzed_foldchange.csv', comment='#')
fc_data = fc_data[(fc_data['temp']==37) & (fc_data['size']=='large') & (fc_data['fold_change'] > 0) & (fc_data['repressors'] > 0)]
inferred_F = pd.read_csv('../../data/ch4_growth/inferred_empirical_F.csv', comment='#')

# Isolate the fc data to the relevant measurements 
inferred_F = inferred_F[inferred_F['temp']==37].copy()

# Compute the summary statistics
rep_summary = fc_data.groupby(['date', 'run_number', 
                               'atc_ngml', 'carbon']).mean().reset_index()
summary = rep_summary.groupby(['atc_ngml', 'carbon']).agg(('mean', 'sem')).reset_index()

#%% Set up the figure axis
fig, ax = plt.subplots(2, 4, figsize=(7, 3.5), dpi=100)
phd.viz.despine(ax.ravel())

# Define the axes and colors
ax_map = {'glucose':0, 'glycerol':1, 'acetate':2}
facecolors = {'glucose':colors['purple'], 'glycerol':colors['green'],
              'acetate':colors['brown']}

for g, d in summary.groupby(['carbon']):
    phd.viz.titlebox(ax[0, ax_map[g]], g, color=colors['black'], size=6, pad=0.05, 
                    boxsize='15%')
    ax[0, ax_map[g]].errorbar(d['repressors']['mean'], d['fold_change']['mean'], 
                      xerr=d['repressors']['sem'], yerr=d['fold_change']['sem'],
                      fmt='o', ms=4, color=facecolors[g], 
                      markeredgecolor=colors['grey'], markeredgewidth=0.75)

for g, d in inferred_F.groupby(['carbon']):
        _reps = summary[(summary['carbon']==g)]
        F = d[d['parameter']=='empirical_F']['median']
        F_max = d[d['parameter']=='empirical_F']['hpd_max']
        F_min = d[d['parameter']=='empirical_F']['hpd_min']
        delF = F - F_ref   
        delF_max = F_max - F_ref   
        delF_min = F_min - F_ref   
        ax[1, ax_map[g]].errorbar(_reps['repressors']['mean'][:-1] / ref_rep, delF, xerr=_reps['repressors']['sem'][:-1] / ref_rep, 
        fmt='o', markeredgecolor=colors['grey'], color=facecolors[g], linestyle='none', linewidth=1,
        ms=4, markeredgewidth=0.75)
        ax[1, ax_map[g]].vlines(_reps['repressors']['mean'][:-1] / ref_rep, delF_max, delF_min, lw=0.75,
                        color=facecolors[g])

for i in range(4):
    ax[0, i].fill_between(rep_range, fc_min,  fc_max, color='k', alpha=0.25)
    ax[0, i].set_xscale('log')
    # ax[1, i].set_xscale('log')
    ax[1, i].plot(rep_range / ref_rep, deltaF, 'k-')
    ax[0, i].set_yscale('log')
    ax[0, i].set_ylabel('fold-change', fontsize=8)
    ax[1, i].set_ylabel('$\Delta F$ [$k_BT$]', fontsize=8)

for i, a in enumerate(ax.ravel()):
    if i <= 3:
        a.set_xlabel('repressors per cell', fontsize=8)
    else: 
        a.set_xlabel('$R_C / R_{ref}$', fontsize=8)



# Plot the explanatory points. 
ax[0, -1].grid(False)
ax[1, -1].grid(False)
phd.viz.titlebox(ax[0, -1], 'example', size=6, color=colors['black'], pad=0.05,
                 boxsize='12%')
ax[0, -1].fill_between(rep_range, example_fc[2], 1.05, color=colors['red'],
            alpha=0.5)
ax[0, -1].fill_between(rep_range, example_fc[2], -.05, color=colors['green'],
            alpha=0.5)
ax[1, -1].fill_betweenx(np.linspace(-5, 5, 200), 1, 4.5, 
                    color=colors['red'], alpha=0.5)
ax[1, -1].fill_betweenx(np.linspace(-5, 5, 100), 0, 0.99, 
                    color=colors['green'], alpha=0.5)
for i, r in enumerate(example_reps):
    if r == ref_rep:
        edge = colors['black']
    else:
        edge = colors['grey']
    ax[0, -1].plot(r, example_fc[i], 'o', markerfacecolor=example_facecolors[r],
            markeredgecolor=edge, markeredgewidth=0.75)
    ax[1, -1].plot(r / ref_rep, example_delF[i], 'o', markerfacecolor=example_facecolors[r],
            markeredgecolor=edge, markeredgewidth=0.75)

ax[1, -1].set_xscale('log')
for i in range(2):
    for j in range(0, 3):
        ax[0, j].set_xlim([2, 300])
        ax[1, j].set_xlim([.01, 2.5])
        ax[1, j].set_xscale('log')
        ax[1, j].set_xticks([0.1, 0.5, 1, 1.5, 2])
        ax[0, j].set_ylim([1E-3, 1.1])
for i in range(4):
    ax[1, i].set_xticks([0.1, 1, 10])

ax[1, -1].set_xlim([0.1, 5])

plt.tight_layout()
plt.savefig('../figs/Fig4_delF_carbon_quality_plots.svg', bbox_inches='tight') 



#%%
