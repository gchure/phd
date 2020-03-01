#%%
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import phd.viz
import phd.thermo
import phd.stats
constants = phd.thermo.load_constants()
colors, palette= phd.viz.phd_style()


# Load and prune data and deltaF
data = pd.read_csv('../../data/ch3_mutants/Chure2019_summarized_data.csv', comment='#')
data = data[data['class']=='DBL'].copy()
deltaF = pd.read_csv('../../data/ch3_mutants/Chure2019_empirical_F_statistics.csv')
deltaF = deltaF[deltaF['class']=='DBL'].copy()

# Load the sampling information
epRA_samples = pd.read_csv('../../data/ch3_mutants/Chure2019_DNA_binding_energy_samples.csv')
epRA_samps = epRA_samples[(epRA_samples['operator']=='O2') & 
                           (epRA_samples['repressors']==260)].copy()
kaki_epai_samples = pd.read_csv('../../data/ch3_mutants/Chure2019_KaKi_epAI_samples.csv')
kaki_epai_samps = kaki_epai_samples[kaki_epai_samples['operator']=='O2'].copy()

#%%
# ##############################################################################
# FIGURE INSTANTIATION AND FORMATTING 
# ##############################################################################
fig, ax = plt.subplots(3, 7, figsize=(7, 3))
phd.viz.despine(ax.ravel())

for i in [1, 2, 5, 6]:
    for j in [0, 1,  2]:
        ax[j, i].spines['left'].set_visible(False)
        ax[j, i].set_yticks([])
    

# Define the mutant axes
DNA_idx = {'Y20I':0, 'Q21A':1, 'Q21M':2}
IND_idx = {'F164T': 0, 'Q294V':1, 'Q294K':2}

for a in ax.ravel():
    a.xaxis.set_tick_params(labelsize=6)
    a.yaxis.set_tick_params(labelsize=6)
    a.set_xscale('symlog', linthreshx=1E-3)
    a.set_xticks([0, 1E-1, 1E1, 1E3])
    a.set_xlim([-0.0005, 1E4])

# Disable middle column
for i in range(3):
    ax[i, 3].axis('off')

# Set scaling
for i in range(3):
    for j in range(3):
        # ax[i, j].set_xscale('symlog', linthreshx=1E-3, linscalex=0.5)
        ax[i, j].set_ylim([-0.1, 1.2])
        ax[i, j+4].set_ylim([-8, 8])
        ax[i, j].set_yticks([0, 0.5, 1])

# Set the titles and axis labels. 
IND_titles = {'F161T':0, 'Q291V':1, 'Q291K':2}
for m, idx in IND_titles.items():
        # phd.viz.titlebox(ax[0, idx].set_title(m, fontsize=6, y=1.04, 
                        #    backgroundcolor=colors['pale_yellow'])
        # ax[0, idx + 4].set_title(m, fontsize=6, y=1.04, 
                        #    backgroundcolor=colors['pale_yellow'])
        ax[-1, idx].set_xlabel('IPTG [µM]', fontsize=6)
        ax[-1, idx + 4].set_xlabel('IPTG [µM]', fontsize=6)

# Remove unnecessary ticklabels
for i in range(2):
    ax[i,0].set_xticklabels([])
    ax[i,4].set_xticklabels([])
    ax[-1, i+1].set_yticklabels([])
    ax[-1, i+5].set_yticklabels([])
    for j in range(2):
        ax[j, i+1].set_xticklabels([])
        ax[j, i+1].set_yticklabels([])
        ax[j, i+5].set_xticklabels([])
        ax[j, i+5].set_yticklabels([])

# Mutant Identifiers
DNA_titles = {'Y17I':0, 'Q18A':1, 'Q18M':2}
for m, idx in DNA_titles.items():
        # ax[idx, 0].text(-0.8, 0.62, m, fontsize=6,rotation='vertical', 
            # backgroundcolor=colors['pale_yellow'], transform=ax[idx, 0].transAxes)
        # ax[idx, 4].text(-0.8, 0.62, m, fontsize=6,rotation='vertical', 
            # backgroundcolor=colors['pale_yellow'], transform=ax[idx, 4].transAxes)
        ax[idx, 0].set_ylabel('fold-change', fontsize=6, labelpad=0.1)
        ax[idx, 4].set_ylabel('$\Delta F$ [$k_BT$]', fontsize=6, labelpad=0.01)

# ##############################################################################
# INDUCTION DATA 
# ##############################################################################
for dna, dna_idx in DNA_idx.items():
    for ind, ind_idx in IND_idx.items():
        # Get the mutant from the summarized data. 
        m_data = data[data['mutant']==f'{dna}-{ind}']

        ax[dna_idx, ind_idx].errorbar(m_data['IPTGuM'], m_data['mean'], m_data['sem'],
            fmt='o', lw=1, capsize=1, linestyle='none', markeredgecolor='white',
            ms=3.5, color=colors['blue'], markeredgewidth=0.5)

# ##############################################################################
# DELTA F DATA
# ##############################################################################
# compute the reference
for dna, dna_idx in DNA_idx.items():
    for ind, ind_idx in IND_idx.items():
        _data = deltaF[(deltaF['mutant']==f'{dna}-{ind}')] 
        for g, d in _data.groupby(['IPTGuM']):
                d
                _d = d[d['parameter']=='delta_bohr']
                _mu = d[d['parameter']=='fc_mu']['median'].values[0]
                _sig = d[d['parameter']=='fc_sigma']['median'].values[0]
                if (_mu < _sig) | (1 - _mu < _sig):
                        c='gray'
                        a = 0
                else:
                        c=colors['blue']
                        a = 1
                if g == 0:
                   cap_max = 1E-4
                   cap_min = -0.0001
                else:
                   cap_max = 1.5 * g
                   cap_min = 0.6 * g
                ax[dna_idx, ind_idx + 4].plot(g, _d['median'], 'o', 
                                color=c, alpha=a, ms=3.5, markeredgewidth=0.5,
                                markeredgecolor=colors['pale_blue'])
                ax[dna_idx, ind_idx + 4].vlines(g, _d['hpd_min'], _d['hpd_max'],
                                        lw=1, color=c, alpha=a)
                ax[dna_idx, ind_idx + 4].hlines(_d['hpd_min'], cap_min, cap_max, 
                                        lw=0.5, color=c, alpha=a)
                ax[dna_idx, ind_idx + 4].hlines(_d['hpd_max'], cap_min, cap_max, 
                                        lw=0.5, color=c, alpha=a)

# ##############################################################################
# THEORY CURVES
# ##############################################################################
c_range = np.logspace(-3, 4, 200)
c_range[0] = 0
ref_bohr = phd.thermo.SimpleRepression(R=260, ep_r=constants['O2'], 
                        ka=constants['Ka'], ki=constants['Ki'], ep_ai=constants['ep_AI'],
                        effector_conc=c_range).bohr_parameter()
n_draws = int(1E4)
for dna, dna_idx in DNA_idx.items():
    for ind, ind_idx in IND_idx.items():
        _allo_samps = kaki_epai_samps[kaki_epai_samps['mutant']==ind]
        _epRA_samps = epRA_samps[epRA_samps['mutant']==dna]
        epRA_draws = np.random.choice(_epRA_samps['ep_RA'].values, replace=True,
                                      size=n_draws)
        allo_draws = _allo_samps.sample(n=n_draws, replace=True)
        fc_cred_region = np.zeros((2, len(c_range)))
        bohr_cred_region = np.zeros((2, len(c_range)))
        for i, c in enumerate(c_range):
            arch = phd.thermo.SimpleRepression(R=260, ep_r=epRA_draws,
                            ka=allo_draws['Ka'], ki=allo_draws['Ki'], 
                            ep_ai=allo_draws['ep_AI'], n_sites=2, effector_conc=c)
            fc_cred_region[:, i] = phd.stats.compute_hpd(arch.fold_change(), 0.95)
            bohr_cred_region[:, i] = phd.stats.compute_hpd(arch.bohr_parameter() - ref_bohr[i], 0.95)

        ax[dna_idx, ind_idx].fill_between(c_range, fc_cred_region[0, :], 
                                        fc_cred_region[1, :], alpha=0.3, 
                                        color=colors['blue'])
        ax[dna_idx, ind_idx + 4].fill_between(c_range, bohr_cred_region[0, :], 
                                        bohr_cred_region[1, :], alpha=0.3, 
                                        color=colors['blue'])
plt.subplots_adjust(wspace=0.1, hspace=0.1)
plt.savefig('../figs/fig5_plots.svg', bbox_inches='tight')

# %%
