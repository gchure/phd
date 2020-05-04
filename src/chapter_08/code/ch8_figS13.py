#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import phd.viz
import phd.stats
colors, color_list = phd.viz.phd_style()

# Define the repressor copy number range. 
rep_range = np.logspace(0, 4, 200)
ref_rep = 100

# Define dynamic constants
ep_ai = 4.5 # in kT
ep_RA = -13.9 # in kT
T_ref = 37 + 273.15 # in K
Nns = 4.6E6 # in bp
sigma = 0.2 # Error in estimation of binding energy

pact_ref = (1 + np.exp(-ep_ai))**-1
bohr_ref = -np.log(pact_ref) - np.log(ref_rep/4.6E6) + ep_RA

# Load the actual data. 
fc_data = pd.read_csv('../../data/ch4_growth/analyzed_foldchange.csv', comment="#")
fc_data = fc_data[fc_data['size']=='large']
fc_data['repressors'] = fc_data['repressors'].round()
inferred_F = pd.read_csv('../../data/ch8_growth_si/inferred_empirical_F.csv', comment='#')
fc_data = fc_data[(fc_data['fold_change'] >= 0) & (fc_data['temp']!=37) ]
independent_params = pd.read_csv('../../data/ch8_growth_si/pooled_entropic_parameter_samples.csv', comment='#')
dependent_params = pd.read_csv('../../data/ch8_growth_si/pooled_entropic_parameter_temp_dependence_samples.csv', comment='#')

# Isolate the fc data to the relevant measurements 
fc_data = fc_data[(fc_data['carbon']=='glucose') & 
                  (fc_data['volume_death'] >= 3.75) & 
                  (fc_data['temp']!=37)] 
inferred_F = inferred_F[inferred_F['temp']!=37].copy()

# Compute the summary statistics
rep_summary = fc_data.groupby(['date', 'run_number', 
                               'atc_ngml', 'temp']).mean().reset_index()
summary = rep_summary.groupby(['atc_ngml', 'temp']).agg(('mean', 'sem')).reset_index()


# Compute the reference bits
ref_pact = (1 + np.exp(-ep_ai))**-1
ref_fc = (1 + ref_pact * rep_range * np.exp(-ep_RA) / Nns)**-1

# Define the temperature colors. 
colors_fill = {32:colors['light_blue'],  42:colors['light_red']}
colors_edge = {32:colors['dark_blue'], 42:colors['dark_red']}

#%%
# Set up the figure canvas
fig, ax = plt.subplots(2, 2, figsize=(4.9, 4.5), dpi=150)
phd.viz.despine(ax.ravel())

# # Add the titles
for i in range(2):
    mwc.viz.titlebox(ax[i, 0], '32° C', color=colors['black'],
            bgcolor='white', size=8, boxsize="12%")
    
    mwc.viz.titlebox(ax[i, 1], '42° C', color=colors['black'],
            bgcolor='white', size=8, boxsize="12%")
            
for i in range(2):
    ax[0, i].set_xscale('log')
    ax[0, i].set_yscale('log')
    ax[0, i].set_xlim([1, 700])
    ax[1, i].set_xscale('log')
    ax[1, i].set_xlim([1E-2, 5])

# Add labels.
for i in range(2):
    ax[0, i].set_xlabel('repressors per cell', fontsize=8)
    ax[0, i].set_ylabel('fold-change', fontsize=8)
    ax[1, i].set_xlabel('$R_{T} / R_{ref}$', fontsize=8)
    ax[1, i].set_ylabel('$\Delta F$ [$k_BT$]')

# Define the axes for temperatures.
temp_axes = {32:0, 42:1}

ax[0, 0].set_xlim([1, 500]) 
ax[0, 1].set_xlim([10, 800]) 
ax[0, 0].set_ylim([2E-4, 1])
ax[0, 1].set_ylim([1E-2, 1])

ax[1, 1].set_xlim([1E-2, 5 ])
ax[1, 0].set_ylim([-3, 3])
ax[1, 1].set_ylim([-2, 4])


params = [dependent_params, independent_params]
labels = ['T independent', 'T dependent']
_colors = [colors['light_grey'], colors['dark_purple']]
for i, p in enumerate(params):
    for g, d in p.groupby(['temp']):
        fc_ax = ax[0, temp_axes[g]]
        delF_ax = ax[1, temp_axes[g]]

        # Extract the relevant parameters. 
        epRA_star = d[(d['parameter']=='epRA_star')]['value'].values
        epAI_star = d[(d['parameter']=='epAI_star')]['value'].values
        pact = (1 + np.exp(-epAI_star))**-1

        # Compute the credible regions for the delta F and fold-change
        fc_cred = np.zeros((2, len(rep_range)))
        delF_cred = np.zeros((2, len(rep_range)))
        for j, r in enumerate(rep_range):
            fc = (1 + pact * r * np.exp(-epRA_star) / 4.6E6)**-1
            bohr = -np.log(pact) - np.log(r / 4.6E6) + epRA_star
            dBohr = bohr - bohr_ref
            fc_cred[:, j] = mwc.stats.compute_hpd(fc, 0.95)
            delF_cred[:, j] = mwc.stats.compute_hpd(dBohr, 0.95)

        if i == 0:
            hatch=""
            facecolor=_colors[i]
            edgecolor=_colors[i]
            alpha=0.5
        else:
            hatch='/////'
            facecolor='none'
            edgecolor=_colors[i]
            alpha=1
        fc_ax.fill_between(rep_range, fc_cred[0, :], fc_cred[1, :], 
                        facecolor=facecolor, alpha=alpha, label = labels[i],
                        hatch=hatch, edgecolor=edgecolor, lw=0.5)
        delF_ax.fill_between(rep_range / ref_rep, delF_cred[0, :], delF_cred[1, :], 
                        facecolor=facecolor,  lw=0.5, alpha=alpha, label=labels[i],
                        hatch=hatch, edgecolor=edgecolor)

# Plot the data
for g, d in summary.groupby(['temp']):
    ax[0, temp_axes[g]].errorbar(d['repressors']['mean'], d['fold_change']['mean'], 
            xerr=d['repressors']['sem'], yerr=d['fold_change']['sem'], fmt='.', 
            linestyle='none', markersize=8, markerfacecolor=colors_fill[g],
            markeredgecolor=colors_edge[g], markeredgewidth=0.75, 
            color=colors_edge[g])


# Plot the empirical delta F vs the repressor copy number
for g, d in inferred_F.groupby(['temp']):
    _reps = summary[summary['temp']==g]['repressors']
    F_median = d[d['parameter']=='empirical_F']['median'].values
    F_min = d[d['parameter']=='empirical_F']['hpd_min'].values
    F_max = d[d['parameter']=='empirical_F']['hpd_max'].values
    F_theo = -np.log(ref_pact) - np.log(ref_rep / Nns) + ep_RA
    delF_median = F_median - F_theo
    delF_max = F_max - F_theo
    delF_min = F_min - F_theo

    ax[1, temp_axes[g]].vlines(_reps['mean'].values[:len(delF_min)] / ref_rep, delF_min, delF_max, linewidth=0.75,
                            color=colors_edge[g])
    ax[1, temp_axes[g]].errorbar(_reps['mean'].values[:len(delF_min)] / ref_rep, delF_median, xerr=_reps['sem'].values[:len(delF_min)] / ref_rep,
                linestyle='none', color=colors_edge[g], markerfacecolor=colors_fill[g],
                markeredgecolor=colors_edge[g], markeredgewidth=0.75, ms=8, fmt='.')

plt.subplots_adjust(wspace=0.3, hspace=0.4)
ax[0, 0].legend(fontsize=6, handlelength=1, ncol=1)
plt.savefig('../figs/ch8_figS13.pdf',  bbox_inches='tight')
plt.savefig('../figs/ch8_figS13.png',  bbox_inches='tight')

# %%
