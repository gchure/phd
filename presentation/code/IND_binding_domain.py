#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import phd.viz
import phd.thermo
import phd.stats
colors, palette = phd.viz.phd_style()
constants = phd.thermo.load_constants()
data = True
free_energy = True
fit = True


# %%
# Load the data sources. 
fc_data = pd.read_csv('../../src/data/ch3_mutants/Chure2019_summarized_data.csv', comment='#')
fc_data = fc_data[fc_data['class']=='IND']
emp_F = pd.read_csv('../../src/data/ch3_mutants/Chure2019_empirical_F_statistics.csv', comment='#')
emp_F = emp_F[emp_F['class']=='IND']
KaKi_epAI = pd.read_csv('../../src/data/ch3_mutants/Chure2019_KaKi_epAI_samples.csv', comment='#')
KaKi_epAI = KaKi_epAI[KaKi_epAI['operator']=='O2']

# %%
# Define the repressor colors
base_colors = ['red', 'orange', 'brown']
ops = ['O1', 'O2', 'O3']
op_edge_colors = {o:colors[f'{c}'] for o, c in zip(ops, base_colors)}


#%% Instantiate the figure canvas. 
fig, ax = plt.subplots(2, 4, figsize=(9, 5.2), dpi=100)
for a in ax.ravel():
    phd.viz.despine(a)

for i in range(4):
    ax[1, i].hlines(0, 0, 1E4, color=colors['black'], alpha=0.75)

# Assign axes and translate names to proper notation. 
axes = {'F164T':0, 'Q294R':1, 'Q294K':2, 'Q294V':3}
names = {'F164T':'F161T', 'Q294R':'Q291R', 'Q294K':'Q291K', 'Q294V':'Q291V'}
for i, a in enumerate(ax.ravel()):
    a.set_xscale('symlog', linthreshx=1E-2) 
    a.set_xlabel('IPTG [µM]', fontsize=10)
    a.set_xlim([-0.001, 1E4])
    a.set_xticks([0, 1E-2, 1E0, 1E2, 1E4])
    a.tick_params(labelsize=9)
    if i <= 3:
        a.set_ylim([-0.03, 1.15])
    else:
        a.set_ylim([-6, 5])
ax[0, 0].set_ylabel('fold-change', fontsize=10)
ax[1, 0].set_ylabel('free energy shift [$k_BT$]', fontsize=10)

for i in range(1, 4):
    ax[0, i].set_yticklabels([])
    ax[1, i].set_yticklabels([])

    # Assign the title boxes. 
for k, v in axes.items():
    for i in range(2):
        _ax = ax[i, v]
        phd.viz.titlebox(_ax, names[k], size=10, color=colors['black'],
                    bgcolor=colors['grey'], pad=0.05, boxsize="15%")

if data == True:
    for g, d in fc_data.groupby(['mutant', 'operator']):
        _ax = ax[0, axes[g[0]]]
        if g[1] == 'O2':
            face = colors['grey']
            edge = op_edge_colors[g[1]]
        else:
            face = op_edge_colors[g[1]]
            edge = colors['grey']
        if fit == False:
            linestyle = 'none'
        else:
            linestyle = 'none'
        _ax.errorbar(d['IPTGuM'], d['mean'], d['sem'], color=op_edge_colors[g[1]],
                    fmt='o', ms=4, linestyle=linestyle, lw=1, markeredgewidth=0.5,
                    markerfacecolor=face, label=f'{constants[g[1]]}',
                    markeredgecolor=edge)

if free_energy == True:
    for g, d in emp_F.groupby(['mutant', 'operator', 'IPTGuM']):
        _ax = ax[1, axes[g[0]]]
        mu = d[d['parameter']=='fc_mu']['median'].values[0]
        sig = d[d['parameter']=='fc_sigma']['median'].values[0]
        if (mu > sig) & ((1 - mu) > sig):
           median, low, high = d[d['parameter']=='delta_bohr'][
                                ['median', 'hpd_min', 'hpd_max']].values[0]
           if g[1] == 'O2':
                zorder = 100
                face = colors['grey']
                edge = op_edge_colors[g[1]]

           else:
                zorder = 1
                edge = colors['grey']
                face = op_edge_colors[g[1]]
           _ax.vlines(d['IPTGuM'], low, high, color=op_edge_colors[g[1]],
                   lw=1, zorder=zorder) 
           _ax.plot(g[-1], median, marker='o', ms=4, zorder=zorder,
                color=op_edge_colors[g[1]], markerfacecolor=face,
                markeredgewidth=0.5, markeredgecolor=edge,  label='__nolegend__')

if fit == True:
    c_range = np.logspace(-4, 4, 200)
    c_range[0] = 0
    for g, d in KaKi_epAI.groupby(['mutant']):
        # Isolate the samples.
        ka = d['Ka'].values
        ki = d['Ki'].values
        epAI = d['ep_AI'].values

        # compute the credible regions for the fold-change
        for i, o in enumerate(['O1', 'O2', 'O3']):
            cred_region = np.zeros((2, len(c_range)))
            for j, c in enumerate(c_range):
                arch = phd.thermo.SimpleRepression(R=constants['RBS1027'], 
                        ep_r=constants[o], ka=ka, ki=ki, ep_ai=epAI,
                        effector_conc=c).fold_change()
                cred_region[:, j] = phd.stats.compute_hpd(arch, 0.95)
            ax[0, axes[g]].fill_between(c_range, cred_region[0, :], 
                    cred_region[1, :], color=op_edge_colors[o], alpha=0.45,
                    label='__nolegend__')

        # Compute the delta F
        wt_pact = phd.thermo.MWC(ka=constants['Ka'], ki=constants['Ki'], 
                                ep_ai=constants['ep_AI'], effector_conc=c_range).pact()
        cred_region = np.zeros((2, len(c_range)))
        for j, c in enumerate(c_range):
            mut_pact = phd.thermo.MWC(ka=ka, ki=ki, ep_ai=epAI, effector_conc=c).pact()
            delF = -np.log(mut_pact / wt_pact[j]) 
            cred_region[:, j] = phd.stats.compute_hpd(delF, 0.95)
        ax[1, axes[g]].fill_between(c_range, cred_region[0, :], 
                    cred_region[1, :], color=colors['orange'],
                    alpha=0.5, label='__nolegend__')

if data == True:
    leg = ax[0, 1].legend(title=r'$\Delta\varepsilon_{RA}$ [$k_BT$]', fontsize=9, 
                         handlelength=0)
    leg.get_title().set_fontsize(9)
plt.subplots_adjust(wspace=0.1, hspace=0.3)

# Determine the plot title and save
if data == True:
    data_name = 'data'
else:
    data_name = 'nodata'
if fit == True:
    fit_name = 'fit'
else:
    fit_name = 'nofit'
if free_energy == True:
    en_name = 'delF'
else:
    en_name = 'nodelF'
plt.savefig(f'../figs/IND_muts_{data_name}_{fit_name}_{en_name}.pdf', 
            bbox_inches='tight', facecolor=None)
#

# %%


# %%
