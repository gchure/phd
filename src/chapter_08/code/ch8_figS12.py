#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import phd.viz
import phd.thermo
import phd.stats
colors, palette = phd.viz.phd_style()

# %%
# Load the data sets and restrict to the carbon sources
data = pd.read_csv('../../data/ch4_growth/analyzed_foldchange.csv', comment='#')
stats = pd.read_csv('../../data/ch8_growth_si/entropic_parameter_samples.csv', comment='#')
pooled_stats = pd.read_csv('../../data/ch8_growth_si/pooled_entropic_parameter_samples.csv', comment='#')
data = data[(data['carbon']=='glucose') & (data['volume_death'] >= 3.75)] 
data = data[data['temp'] !=37]
summary = data.groupby(['temp', 'date', 'run_number', 'atc_ngml']).mean().reset_index()
summary = summary.groupby(['temp', 'atc_ngml']).agg(('mean', 'sem')).reset_index()
stats = stats[(stats['temp']!=37)]

# Define the constants for plotting
rep_range = np.logspace(-1, 3, 100)
# %%
# Set up the figure canvas
fig, ax = plt.subplots(2, 2, figsize=(4, 4), dpi=100)
phd.viz.despine(ax.ravel())
for a in ax.ravel():
    a.set_xscale('log')
    a.set_yscale('log')
    a.set_xlim([0.5, 800])

ax[0, 0].spines['bottom'].set_visible(False)
ax[0, 1].spines['bottom'].set_visible(False)
ax[0, 1].set_yticks([])
ax[0, 1].set_xticks([])
ax[0, 1].spines['left'].set_visible(False)
ax[1, 1].spines['left'].set_visible(False)
ax[1, 1].set_yticks([])
ax[0, 0].set_xticks([])



for i in range(2):
    ax[-1, i].set_xlabel('repressors per cell')
    ax[i, 0].set_ylabel('fold-change')

titles = ['32° C', '42° C']
title_colors = [colors['dark_blue'], colors['dark_red']]
bgcolors = [colors['pale_blue'], colors['pale_red']]
face_colors = [colors['light_blue'], colors['light_red']]
for i in range(2):
    if i > 0:
        # apply offset transform to all y ticklabels.
        dx = -13 / fig.dpi
        dy = 0
        offset = matplotlib.transforms.ScaledTranslation(dx, dy, fig.dpi_scale_trans)
        for label in ax[i, 0].yaxis.get_majorticklabels():
            label.set_transform(label.get_transform() + offset)


# Plot the predictions
for i, fit_temp in enumerate(titles):
    fit_temp = int(fit_temp.split('°')[0]) 

    # Get the binding energy values for the prediction strain 
    delta_SR = stats[(stats['parameter']=='delta_SR') & 
                    (stats['temp']==fit_temp)]['value'].values
    delta_SAI = stats[(stats['parameter']=='delta_SAI') & 
                    (stats['temp']==fit_temp)]['value'].values
    for j, pred_temp in enumerate(titles): 
        pred_temp = int(pred_temp.split('°')[0])
        
        # Compute the adjusted energies. 
        rel_T = (37 + 273.15) / (pred_temp + 273.15)
        if i == j:
            epRA_star = stats[(stats['parameter']=='epRA_star') & 
                              (stats['temp']==fit_temp)]['value'].values
            epAI_star = stats[(stats['parameter']=='epAI_star') & 
                              (stats['temp']==fit_temp)]['value'].values
        else:
            epRA_star =  delta_SR * ((37 + 273.15) - (pred_temp + 273.15)) -13.9
            epAI_star = delta_SAI * ((37 + 273.15) - (pred_temp + 273.15)) + 4.5

        p_act = (1 + np.exp(-epAI_star))**-1

        # Compute the credible region
        cred_region = np.zeros((2, len(rep_range)))
        for k, r in enumerate(rep_range):
            theo =  (1 + p_act * (r / 4.6E6) * np.exp(-epRA_star))**-1 
            cred_region[:, k] = phd.stats.compute_hpd(theo, 0.95)

        ax[i, j].fill_between(rep_range, cred_region[0, :], cred_region[1, :], 
                            color=title_colors[i], alpha=0.25)
        # Plot the simple rescaling prediction.
        relT = (37 + 273.15) / (pred_temp + 273.15)
        p_act = (1 + np.exp(-relT * 4.5))**-1
        theo =  (1 + p_act * (rep_range / 4.6E6) * np.exp(-(relT * -13.9)))**-1 
        ax[i, j].plot(rep_range, theo, 'k--') 

# Plot the pooled prediction
for i, t1 in enumerate([32, 42]):
    for j, t2 in enumerate([32, 42]):
        _samples = pooled_stats[pooled_stats['temp']==t2]
        epRA_star = _samples[_samples['parameter']=='epRA_star']['value'].values
        epAI_star = _samples[_samples['parameter']=='epAI_star']['value'].values
        cred_region = np.zeros((2, len(rep_range)))
        p_act = (1 + np.exp(-epAI_star))**-1
        for k, r in enumerate(rep_range):
            theo = (1 + p_act * (r / 4.6E6) * np.exp(-epRA_star))**-1
            cred_region[:, k] = phd.stats.compute_hpd(theo, 0.95)

        ax[i, j].fill_between(rep_range, cred_region[0, :], cred_region[1, :], color='grey', alpha=0.5) 


# Plot the data
for i, temp in enumerate(titles):
    temp = int(temp.split('°')[0])
    for j in range(2):
        if i == j:
            fill  = 'white'
            edge = face_colors[i]
        else:
            fill = face_colors[i]
            edge = colors['grey']

        # Isolate the data. 
        d = summary[summary['temp']==temp]
        ax[j, i].errorbar(d['repressors']['mean'], d['fold_change']['mean'],
                        xerr=d['repressors']['sem'], yerr=d['fold_change']['sem'],
                        color=face_colors[i], markeredgecolor=edge, fmt='o', ms=5, 
                        markerfacecolor=fill, markeredgewidth=0.75, linestyle='none', lw=0.75)
plt.subplots_adjust(wspace=0.05, hspace=0.05)
plt.savefig('../figs/figS12_plots.svg', bbox_inches='tight')

# %%
