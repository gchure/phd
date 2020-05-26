#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import phd.viz
import phd.thermo
colors, palette = phd.viz.phd_style()
constants = phd.thermo.load_constants()
# %%
# Define what should be plotted
data = False
free_energy = False
fit = False

# Define the repressor colors
base_colors = ['blue', 'green', 'orange', 'purple']
reps = [60, 124, 260, 1220]
rep_edge_colors = {r:colors[f'{c}'] for r, c in zip(reps, base_colors)}



# %%
# Load the fold-change, empirical F, and inference data
fc_data = pd.read_csv('../../src/data/ch3_mutants/Chure2019_summarized_data.csv', comment='#')
fc_data = fc_data[fc_data['class']=='DNA']
F_data = pd.read_csv('../../src/data/ch3_mutants/Chure2019_empirical_F_statistics.csv', comment='#')
F_data = F_data[F_data['class']=='DNA']
epRA = pd.read_csv('../../src/data/ch3_mutants/Chure2019_DNA_binding_energy_summary.csv', comment='#')
epRA = epRA[epRA['repressors']==260]
# %%
fig, ax = plt.subplots(2, 3, figsize=(7, 5), dpi=100)
for a in ax.ravel():
    phd.viz.despine(a)

for i in range(3):
    ax[1, i].hlines(0, 0, 1E4, color=colors['black'], lw=0.75, alpha=0.75)
# Define the mutant axes
axes = {'Y20I':1, 'Q21M':0, 'Q21A':2}

# Define the renamed mutants.
rename = {'Y20I':'Y17I', 'Q21M':'Q18M', 'Q21A':'Q18A'}

# Format the axes
for i, a in enumerate(ax.ravel()):
    a.set_xscale('symlog',linthreshy=1E-2)
    a.tick_params(labelsize=9)
    a.set_xlabel('IPTG [µM]', fontsize=10)
    a.set_xlim([-0.2, 1E4])
    if i == 0:
        a.set_ylabel('fold-change', fontsize=10)
    elif i==3:
        a.set_ylabel('free energy shift [$k_BT$]', fontsize=10)
    else:
        a.set_yticklabels([])
    if i < 3:
        a.set_ylim([-0.05, 1.15])
    else:
        a.set_ylim([-8, 8])

for k, v in axes.items():
    for i in range(2):
        phd.viz.titlebox(ax[i, axes[k]], rename[k], size=10, color=colors['black'],
                    boxsize="12.5%", bgcolor=colors['grey'], pad=0.05)

# Plot the fold-change data 
if data == True:
    for g, d in fc_data.groupby(['mutant', 'repressors']):
        if g[1] == 260:
            face = colors['grey']
            edge = rep_edge_colors[g[1]]
        else:
            face = rep_edge_colors[g[1]]
            edge = colors['grey']
        if fit == False:
            linestyle='none'
        else:
            linestyle= 'none'
        ax[0, axes[g[0]]].errorbar(d['IPTGuM'], d['mean'], d['sem'], fmt='o',
                                   ms=5, color=rep_edge_colors[g[1]], markerfacecolor=face,
                                   markeredgewidth=0.5, markeredgecolor=edge, linestyle=linestyle, label=int(g[1]),
                                   lw=0.75) 

# Plot the empirical F
if free_energy == True:
    for g, d in F_data.groupby(['mutant', 'repressors', 'IPTGuM']):
        median, low, high = d[d['parameter']=='delta_bohr'][
                            ['median', 'hpd_min', 'hpd_max']].values[0]
        mu = d[d['parameter']=='fc_mu']['median'].values[0]
        sig = d[d['parameter']=='fc_sigma']['median'].values[0]
        if g[1] == 260:
            face = colors['grey']
            edge = rep_edge_colors[g[1]]
            zorder=1000
        else:
            face = rep_edge_colors[g[1]]
            edge = colors['grey']
            zorder = 1
        if fit == False:
            linestyle='--'
        else:
            linestyle= 'none'

        if (mu > sig) & ((1 - mu) > sig):
            ax[1, axes[g[0]]].vlines(g[-1], low, high, color=rep_edge_colors[g[1]],
                                    lw=1, zorder=zorder)
            ax[1, axes[g[0]]].plot(g[-1], median, color=rep_edge_colors[g[1]],
                                  marker='o', ms=4.5, linestyle=linestyle,
                                  markerfacecolor=face, markeredgewidth=0.75, 
                                  zorder=zorder, markeredgecolor=edge)

# Plot the fit. 
c_range = np.logspace(-2, 4, 200)
c_range[0] = 0
if fit == True:
    for g, d in epRA.groupby(['mutant']):
        low, high = d[d['parameter']=='ep_RA'][['hpd_min', 'hpd_max']].values[0]

        # Plot the empirical F first
        delF_low = constants['O2'] - low 
        delF_high = constants['O2'] - high

        # Plot the delF cred region 
        ax[1, axes[g]].fill_between(c_range, -delF_low, -delF_high, color=colors['orange'],
                                    alpha=0.75, label='__nolegend__')
        for i, r in enumerate(rep_edge_colors.keys()):
            # Compute the theoretical fold-change. 
            c, ep = np.meshgrid(c_range, np.array([low, high]))
            arch = phd.thermo.SimpleRepression(R=r, ep_r=ep, ka=constants['Ka'],
                                            ki=constants['Ki'], ep_ai=constants['ep_AI'],
                                            effector_conc=c_range).fold_change()
            ax[0, axes[g]].fill_between(c_range, arch[0, :], arch[1, :], 
                            color=rep_edge_colors[r], alpha=0.45, label='__nolegend__') 

if data == True:
    leg = ax[0, 0].legend(title='repressors', fontsize=9, handlelength=0)
    leg.get_title().set_fontsize(9)
plt.subplots_adjust(wspace=0.1, hspace=0.4)

if fit == True:
    fit_name = 'fit'
else:
    fit_name = 'nofit'
if free_energy == True:
    data_name = 'ind'
else:
    data_name = 'ind_empF'
if data == True:
    plot_id = 'data'
else:
    plot_id = 'nodata'
plt.savefig(f'../figs/DNA_muts_{fit_name}_{data_name}_{plot_id}.pdf', bbox_inches='tight',
            facecolor=None)


# %%
