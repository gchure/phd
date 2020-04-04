#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import phd.thermo
import phd.stats
import phd.viz
_ = phd.viz.phd_style()
color = phd.viz.color_selector('mut')
constants = phd.thermo.load_constants()

#  Load the data
data = pd.read_csv('../../data/ch3_mutants/Chure2019_summarized_data.csv', comment='#')
data = data[(data['class'] == 'DNA') & (data['operator']=='O2')]
samples = pd.read_csv('../../data/ch7_mutants_si/DNA_binding_energy_samples.csv')
stats = pd.read_csv('../../data/ch7_mutants_si/DNA_binding_energy_summary.csv')

# Determine the unique repressor copy numbers
reps = np.sort(data['repressors'].unique())
c_range = np.logspace(-3, 4, 200)
c_range[0] = 0
fig, ax = plt.subplots(len(reps), len(reps), figsize=(6,6)) 
phd.viz.despine(ax.ravel())
# Format the axes
for a in ax.ravel():
    a.set_xscale('symlog', linthreshx=1E-2)
    a.set_xlim([-0.001, 1E4])

# Add appropriate labels
for i in range(len(reps)):
    ax[i, 0].set_ylabel('fold-change')
    ax[-1, i].set_xlabel('IPTG [µM]')
    ax[i, 0].text(-0.47, 0.38, '$R = $' + str(int(reps[i])), fontsize=8,
                 transform=ax[i,0].transAxes, rotation='vertical')
    ax[0, i].set_title('$R = $' + str(int(reps[i])), fontsize=8) 
for i in range(4):
    ax[-1, i].set_xticks([0, 1E-1, 1E1, 1E3])
    
# Add predictor titles
fig.text(0, 0.47, 'fit strain', fontsize=8, bbox={'boxstyle':'square', 'edgecolor':'k', 'facecolor':'white'}, rotation='vertical')
fig.text(0.435, 0.92, 'comparison strain', fontsize=8, bbox={'boxstyle':'square', 'edgecolor':'k', 'facecolor':'white'})
    
# Plot the data. 
for g, d in data.groupby(['mutant']):
    g = g.upper()
    for i, _ in enumerate(reps):
        for j, _ in enumerate(reps):
            _d = d[d['repressors'] == reps[j]]
            if i == j:
                face = 'w'
                edge = color[g]
            else:
                face = color[g]
                edge = 'white' 
            _ = ax[i, j].errorbar(_d['IPTGuM'], _d['mean'], _d['sem'], markerfacecolor=face,
                                 markeredgecolor=edge, color=edge, lw=0.15, linestyle='none', fmt='o',
                                 ms=4, label=f'{g[0]}{str(int(g[1:-1]) - 3)}{g[-1]}', markeredgewidth=0.5)
           
            # Plot the best-fit lines. 
            for k, m in enumerate(data['mutant'].unique()):
                _d = data[(data['mutant']=='m') & (data['repressors'] == reps[j])]
                
                # Get the binding energies. 
                epRA = stats[(stats['parameter']=='ep_RA') & (stats['mutant']==m) & (stats['repressors']==reps[i])][['median', 'hpd_min', 'hpd_max']].values[0]

                
                # Compute the fold-change
                epRA_mesh, c_mesh = np.meshgrid(epRA, c_range)
                fc = phd.thermo.SimpleRepression(R=reps[j], ep_r=epRA_mesh, ka=constants['Ka'], 
                                                 ki=constants['Ki'], ep_ai=constants['ep_AI'],
                                                effector_conc=c_mesh, n_sites=constants['n_sites'],
                                                n_ns=constants['Nns']).fold_change()
                
                # Plot the fit. 
                _ = ax[i, j].fill_between(c_range, fc[:, 1], fc[:, 2], color=color[m], alpha=0.2) 

            
_  = ax[0, 3].legend(fontsize=6, bbox_to_anchor=(1.04, 0.95))
for i in range(3):
    for j in range(4):
        ax[i, j].set_xticks([])
        ax[i, j].spines['bottom'].set_visible(False)
        if (j > 0):
            ax[i, j].set_yticks([])
            ax[i, j].spines['left'].set_visible(False)

for i in range(1, 4):
    ax[-1, i].set_yticks([])
    ax[-1, i].spines['left'].set_visible(False)
plt.subplots_adjust(wspace=0.05, hspace=0.05)
plt.savefig('../figs/ch7_figS11.pdf', bbox_inches='tight')
plt.savefig('../figs/ch7_figS11.png', bbox_inches='tight')


# %%
