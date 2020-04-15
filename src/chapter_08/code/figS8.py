#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import phd.viz
import phd.stats
import phd.thermo
import matplotlib.gridspec as gridspec
colors, palette = phd.viz.phd_style()

# %%
# load the various data sets
old_gods = pd.read_csv('../../data/other/Garcia2011_Brewster2014.csv', comment='#')
new_gods = pd.read_csv('../../data/ch2_induction/RazoMejia_2018.csv', comment='#')
data = pd.read_csv('../../data/ch4_growth/analyzed_foldchange.csv', comment='#')

# Prune the data as necessary
old_gods.rename(columns={'repressor':'repressors'}, inplace=True)
old_gods = old_gods[old_gods['operator']=='O2']

new_gods = new_gods[(new_gods['repressors'] > 0) & (new_gods['IPTG_uM']==0) & 
                    (new_gods['operator']=='O2')]
new_gods.rename(columns={'fold_change_A':'fold_change'}, inplace=True)
new_gods['repressors'] *= 2
data = data[(data['carbon']=='glucose') & (data['temp']==37) & 
            (data['strain']=='dilution')]
large = data[data['size']=='large']
# Group the data for display. 
new_gods = new_gods.groupby('repressors').agg(('mean', 'sem')).reset_index()
data = data.groupby(['date', 'run_number', 'atc_ngml']).mean().reset_index()
data = data.groupby(['atc_ngml']).agg(('mean', 'sem')).reset_index()
large = large.groupby(['date', 'run_number', 'atc_ngml']).mean().reset_index()
large = large.groupby(['atc_ngml']).agg(('mean', 'sem')).reset_index()


# Load the parameter estimate summary
stats = pd.read_csv('../../data/ch8_growth_si/systematic_error_DNA_binding_energy_summary.csv', comment="#")

# %%
# Set teh expreimental parameters
rep_range = np.logspace(-1, 4, 200)
epRA = -13.9
sig = 0.2
epRA_low = epRA - sig
epRA_high = epRA + sig
theo_min = phd.thermo.SimpleRepression(R=rep_range, ep_r=epRA_low, effector_conc=0,
                                      ka=0.53, ki=139, ep_ai=4.5).fold_change()
theo_max = phd.thermo.SimpleRepression(R=rep_range, ep_r=epRA_high, effector_conc=0,
                                      ka=0.53, ki=139, ep_ai=4.5).fold_change()
# %%

# Set up the figure canvas
fig = plt.figure(figsize=(6, 4), dpi=100)
gs = gridspec.GridSpec(1, 5)
ax0 = fig.add_subplot(gs[:3])
ax1 = fig.add_subplot(gs[3:])
ax = [ax1, ax0]
phd.viz.despine(ax)

# Set the appropriate scaling
ax[1].set_xscale('log')
ax[1].set_yscale('log')
ax[1].set_xlim([1, 2000])
ax[0].set_ylim([-0.5, 6.5])
ax[0].set_xlim([-16, -13])


#Add the appropriate labels
phd.viz.titlebox(ax[0], 'DNA binding energy estimation', color=colors['black'], 
                 size=8, boxsize="5%")
phd.viz.titlebox(ax[1], 'theoretical prediction', color=colors['black'], 
                 size=8, boxsize="5%")

# Define the positioning of the data sets. 
position = {'correction': 1, 'no_correction':0, 'large_only':2, 'all_divided': 3,
            'garcia':4, 'brewster':5, 'razo-mejia':6}

# Define the colors. 
fill_colors = ['purple', 'purple', 'black', 'red', 'blue', 
                'green', 'orange']

fill_colors = {k:colors[v] for k, v in zip(position.keys(), fill_colors)}


# Add axis labels. 
ax[0].set_xlabel(r'DNA binding energy [$k_BT$]', fontsize=8)
ax[0].set_yticks([0, 1, 2, 3, 4, 5, 6])
ax[0].set_yticklabels(['this work\n(no correction)', 
                       'this work\n(correction applied\nto small cells)',
                       'this work\n(large cells only,\nno correction applied)',
                       'this work\n(correction applied\nto all cells)',
                       'Garcia & Phillips\n2011',
                       'Brewster et al.\n 2014',
                       'Razo-Mejia et al.\n 2018'])
ax[1].set_xlabel('repressors per cell')
ax[1].set_ylabel('fold-change')

# Plot the "gold standard" binding energy. 
ax[0].fill_betweenx([-1, 10], epRA_low, epRA_high, color=colors['black'], alpha=0.25, zorder=1)

# Plot the inferred DNA binding energies. 
for g, d in stats.groupby('source'):
    if g == 'no_correction':
        alpha = 0.25
    else:
        alpha = 1
    epRA_median, low, high = d[d['parameter']=='epRA'][
        ['median', 'hpd_min', 'hpd_max']].values[0]
    ax[0].hlines(position[g], low, high, lw=1, color=fill_colors[g],
                alpha=alpha, zorder=99)
    ax[0].plot(epRA_median, position[g], lw=1, marker='o', 
               markerfacecolor=fill_colors[g], markeredgecolor='white',
               markeredgewidth=0.5, ms=4.5, alpha=alpha, zorder=100)

# Plot the theory
iter = 2
ax[1].fill_between(rep_range, theo_min, theo_max, color=colors['black'], 
                  alpha=0.25, label='__nolegend__',
                  zorder=iter)

# Plot the fold-change data sets. 
# Brewster and Garcia
for g, d in old_gods.groupby(['author']):
    if g == 'garcia':
        label = 'Garcia & Phillips 2011'
    else:
        label = 'Brewster et al. 2014'
    ax[1].plot(d['repressors'], d['fold_change'], 'o', ms=4.5, 
           markerfacecolor=fill_colors[g], markeredgecolor='white',
           markeredgewidth=0.5, label=label, zorder=100, alpha=0.75)
    iter +=1

# Razo-Mejia
ax[1].errorbar(new_gods['repressors'], new_gods['fold_change']['mean'], 
               new_gods['fold_change']['sem'], lw=0.75, linestyle='none',
               capsize=1, color=fill_colors['razo-mejia'], fmt='o', ms=4.5,
               label='Razo-Mejia et al. 2018',
               zorder=iter,  alpha=0.75, markeredgecolor='white', markeredgewidth=0.5)
iter += 1
# Without correction
ax[1].errorbar(data['raw_repressors']['mean'], data['fold_change']['mean'],
            xerr=data['raw_repressors']['sem'], yerr=data['fold_change']['sem'],
            fmt='o', ms=4.5, lw=0.75, capsize=1, linestyle='none', 
            markerfacecolor=fill_colors['no_correction'],
            color=fill_colors['no_correction'], markeredgewidth=0.5,
            label='without correction', alpha=0.25, zorder=iter, markeredgecolor='white')

iter += 1
# With correction
ax[1].errorbar(data['repressors']['mean'], data['fold_change']['mean'],
            xerr=data['repressors']['sem'], yerr=data['fold_change']['sem'],
            fmt='o', ms=4.5, lw=0.75, capsize=1, linestyle='none', 
            markerfacecolor=fill_colors['correction'],
            color=fill_colors['correction'], markeredgewidth=0.75,
            label='correction applied\nto small cells', zorder=iter,
            markeredgecolor='white')

iter += 1
# large only
ax[1].errorbar(large['repressors']['mean'], large['fold_change']['mean'],
            xerr=large['repressors']['sem'], yerr=large['fold_change']['sem'],
            fmt='o', ms=4.5, lw=0.75, capsize=1, linestyle='none', 
            color=fill_colors['large_only'], markeredgewidth=0.5,
            markeredgecolor='white',
            label='large cells only,\nno correction applied', zorder=iter)
# all divided
ax[1].errorbar(data['repressors']['mean'] * 2, data['fold_change']['mean'],
            xerr=data['repressors']['sem']*2, yerr=data['fold_change']['sem'],
            fmt='o', ms=4.5, lw=0.75, capsize=1, linestyle='none', 
            color=fill_colors['all_divided'], markeredgewidth=0.5,
            label='correction applied\nto all cells', zorder=iter,
            markeredgecolor='white')

# Add legends
ax[1].legend()

# Add panel labels
fig.text(0.05, 0.98, '(A)', fontsize=8)
fig.text(0.55, 0.98, '(B)', fontsize=8)
ax[0].spines['right'].set_visible(True)
ax[0].spines['left'].set_visible(False)
ax[0].yaxis.tick_right()
plt.tight_layout()
plt.savefig('../figs/ch8_figS8.pdf', bbox_inches='tight')
plt.savefig('../figs/ch8_figS8.png', bbox_inches='tight')


# %%
