
#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import phd.thermo
import phd.stats
import phd.viz
color = phd.viz.color_selector('mut')
constants = phd.thermo.load_constants()
_ = phd.viz.phd_style()

#  Load the data
data = pd.read_csv('../../data/ch3_mutants/Chure2019_summarized_data.csv', comment='#')
data = data[data['class'] == 'IND']
KaKi_only_samples = pd.read_csv('../../data/ch3_mutants/Chure2019_KaKi_only_samples.csv', comment='#')
KaKi_epAI_samples = pd.read_csv('../../data/ch3_mutants/Chure2019_KaKi_epAI_samples.csv', comment='#')


# Determine the unique repressor copy numbers
ops = np.sort(data['operator'].unique())
c_range = np.logspace(-2, 4, 200)
c_range[0] = 0
MODEL = 'KaKi_only'

# ##############################################################################
#  FIGURE WITH KAKI FIT ONLY
# ##############################################################################
fig, ax = plt.subplots(len(ops),len(ops), figsize=(6,5)) 
phd.viz.despine(ax.ravel())
# Format the axes
for a in ax.ravel():
    a.set_xscale('symlog', linthreshx=1E-2)
    a.set_ylim([0, 1.1])

for i in range(1, 3):
    ax[0, i].set_xticks([])
    ax[0, i].set_yticks([])
    ax[0, i].spines['bottom'].set_visible(False)
    ax[0, i].spines['left'].set_visible(False)

for i in range(2):
    ax[0, i].set_xticks([])
    ax[0, i].spines['bottom'].set_visible(False)
    ax[-1, i+1].set_yticks([])
    ax[-1, i+1].spines['left'].set_visible(False)
ax[1, 1].set_yticks([])
ax[1, 2].set_yticks([])
ax[1, 1].spines['left'].set_visible(False)
ax[1, 2].spines['left'].set_visible(False)
ax[1, 1].set_xticks([])
ax[1, 2].set_xticks([])
ax[1, 1].spines['bottom'].set_visible(False)
ax[1, 2].spines['bottom'].set_visible(False)
ax[1, 0].set_xticks([])
ax[1, 0].spines['bottom'].set_visible(False)



# Add appropriate labels
for i in range(len(ops)):
    ax[i, 0].set_ylabel('fold-change', fontsize=8)
    ax[-1, i].set_xlabel('IPTG [µM]', fontsize=8)
    ax[i, 0].text(-0.5, 0.55, ops[i], fontsize=8,
                 transform=ax[i,0].transAxes, rotation='vertical')
    ax[0, i].set_title(ops[i], fontsize=8)
for i in range(3):
    ax[-1, i].set_xticks([0, 1E-2, 1E0, 1E2, 1E4])
    
# Add predictor titles
fig.text(-0.04, 0.5, 'fit strain', fontsize=8, rotation='vertical', 
        bbox={'edgecolor':'k', 'facecolor':'w'})
fig.text(0.435, 0.94, 'comparison strain', fontsize=8,  rotation='horizontal', 
        bbox={'edgecolor':'k', 'facecolor':'w'})
    
# Plot the data. 
for g, d in data.groupby(['mutant']):
    g = g.upper()
    for i, _ in enumerate(ops):
        for j, _ in enumerate(ops):
            _d = d[d['operator'] == ops[j]]
            if i == j:
                face = 'w'
                edge = color[g]
            else:
                face = color[g]
                edge = 'w'
            _ = ax[i, j].errorbar(_d['IPTGuM'], _d['mean'], _d['sem'], markerfacecolor=face,
                                 markeredgecolor=edge, color=color[g], lw=0.75, linestyle='none', fmt='o',
                                 ms=4, label=g, markeredgewidth=0.5)
           
            # Plot the best-fit lines. 
            for k, m in enumerate(data['mutant'].unique()):
                _d = data[(data['mutant']==m) & (data['operator'] == ops[i])]
                
                # Get the binding energies. 
                if MODEL == 'KaKi_only':
                    _samps = KaKi_only_samples[(KaKi_only_samples['mutant']==m) &\
                         (KaKi_only_samples['operator']==ops[i])]
                    _samps['ep_AI'] = 4.5
                else:
                    _samps = KaKi_epAI_samples[(KaKi_epAI_samples['mutant']==m) &\
                        (KaKi_epAI_samples['operator']==ops[i])]
                Ka = _samps['Ka'][::10]
                Ki = _samps['Ki'][::10]
                epAI = _samps['ep_AI'][::10]
                cred_region = np.zeros((2, len(c_range)))
                for z, c in enumerate(c_range): 
                    # Compute the fold-change   
                    fc = phd.thermo.SimpleRepression(R=constants['RBS1027'], 
                                            ep_r=constants[ops[j]], ka=Ka, 
                                                 ki=Ki, ep_ai=epAI,
                                                effector_conc=c, n_sites=constants['n_sites'],
                                                n_ns=constants['Nns']).fold_change()
                    cred_region[:, z] = phd.stats.compute_hpd(fc, 0.95)
                
                # Plot the fit. 
                _ = ax[i, j].fill_between(c_range, cred_region[0, :], 
                                   cred_region[1, :], color=color[m], alpha=0.2) 

_  = ax[0, 2].legend(fontsize=8, bbox_to_anchor=(1.04, 0.95))
plt.subplots_adjust(wspace=0.05, hspace=0.05)

if MODEL == 'KaKi_only':
    name = 'ch7_figS17'
else:
    name = 'ch7_figS18'
plt.savefig(f'../figs/{name}.pdf', bbox_inches='tight')
plt.savefig(f'../figs/{name}.png', bbox_inches='tight')

    

# %%
