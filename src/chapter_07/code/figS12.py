#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import phd.thermo
import phd.viz
constants = phd.thermo.load_constants()
colors = phd.viz.color_selector('pboc')
_ = phd.viz.phd_style()

# Load the empirical F data
data = pd.read_csv('../../data/ch3_mutants/Chure2019_empirical_F_statistics.csv')
data = data[data['class']=='DNA']
epRA_stats = pd.read_csv('../../data/ch3_mutants/Chure2019_DNA_binding_energy_summary.csv')

# ##############################################################################
#  FIGURE INSTANTIATION
# ##############################################################################
fig, ax = plt.subplots(4, 3, figsize=(6, 4))
phd.viz.despine(ax.ravel())
for a in ax.ravel():
    a.set_xscale('symlog', linthreshx=2E-2)
    a.set_ylim([-8, 8])

for i in range(3):
    for j in range(3):
        ax[i, j].set_xticks([])
        ax[i, j].spines['bottom'].set_visible(False)
        if (i > 0) and (j > 0):
            ax[i, j].set_yticks([])
            ax[i, j].spines['left'].set_visible(False)

for i in range(1, 3):
    ax[0, i].set_yticks([])
    ax[0, i].spines['left'].set_visible(False)  
    ax[-1, i].set_yticks([])
    ax[-1, i].spines['left'].set_visible(False)

for i in range(3):
    ax[-1, i].set_xlabel('IPTG [µM]')

for i in range(4):
    ax[i, 0].set_ylabel('$\Delta F$ [$k_BT$]')

# Define the axes. 
cols = {'Q21A':1, 'Q21M':0, 'Y20I':2}
rows = {60:0, 124:1, 260:2, 1220:3}

# Add titles and row labels. 
for m, a in cols.items():
    ax[0, a].set_title(m, fontsize=6, bbox={'edgecolor':'k', 'facecolor':'white'})

for r, a in rows.items():
    ax[a, 0].text(-0.38, 0.40, f'R = {int(r)}', fontsize=6, rotation='vertical',
                backgroundcolor=colors['pale_yellow'],
                transform=ax[a, 0].transAxes, 
                bbox={'edgecolor':'k', 'facecolor':'white'})

# Define the repressor colors
rep_colors = {60:colors['purple'], 124:colors['blue'], 
              260:colors['red'], 1220:colors['green']}

# ##############################################################################
# PLOT LEGEND ENTRIES
# ##############################################################################
for r, c in rep_colors.items():
    ax[0, 0].plot([], [], 'o', ms=4, color=c, label=int(r), markeredgecolor='white',
                  markeredgewidth=0.5)
leg = ax[0,  0].legend(title='repressors per cell', fontsize=6, ncol=4,
                    columnspacing=0.001, handletextpad=0.1)
leg.get_title().set_fontsize(6)

# ##############################################################################
# PLOT DATA
# ##############################################################################
for g, d in data.groupby(['mutant', 'repressors', 'IPTGuM']):
    dbohr = d[d['parameter']=='delta_bohr']
    fc_mu = d[d['parameter']=='fc_mu']['median'].values[0]
    fc_sig = d[d['parameter']=='fc_sigma']['median'].values[0]
    if (fc_mu > fc_sig) & (1 - fc_mu > fc_sig):
        for r, a in rows.items():
            _ax = ax[a, cols[g[0]]]
            if g[1] == r:
                face = 'w'
                edge = rep_colors[g[1]]
                zorder = 1000
            else:
                face = rep_colors[g[1]]
                edge = 'w'
                zorder = 100
            _ax.plot(dbohr['IPTGuM'], dbohr['median'], 'o', color=face,
                     markeredgecolor=edge, markeredgewidth=0.75, label='__nolegend__',
                     zorder=zorder, ms=3.5)
            _ax.vlines(dbohr['IPTGuM'], dbohr['hpd_min'], dbohr['hpd_max'], lw=0.75,
                color=rep_colors[g[1]], label='__nolegend__', zorder=1)


# ##############################################################################
#  PREDICTIONS
# ##############################################################################
for c, ca in cols.items():
    for r, ra in rows.items():
        _ax = ax[ra, ca]
        epRA = epRA_stats[(epRA_stats['parameter']=='ep_RA') & 
        (epRA_stats['mutant']==c) & (epRA_stats['repressors']==r)]
        
        _ax.fill_between([0, 1E4], epRA['hpd_min'] - constants['O2'], 
                        epRA['hpd_max'] - constants['O2'], alpha=0.75, 
                        color=rep_colors[r])
plt.subplots_adjust(hspace=0.05, wspace=0.03)
plt.savefig('../figs/ch7_figS12.pdf', bbox_inches='tight')

# %%
