#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import phd.stats
import phd.thermo
import phd.viz
colors, palette = phd.viz.phd_style()
constants = phd.thermo.load_constants()

# %%
# Load, restrict, and clean the various data sets. 
data = pd.read_csv('../../data/ch4_growth/analyzed_foldchange.csv', comment='#')
data = data[(data['fold_change'] > 0) & (data['repressors'] > 0) &
            (data['size']=='large')]
data = data.groupby(['date', 'carbon', 'temp', 
                    'atc_ngml'])[['fold_change', 'repressors']].mean()
data = data.groupby(['carbon', 'temp', 'atc_ngml']).agg(
                    ('mean', 'sem')).reset_index()

ind_data = pd.read_csv('../../data/ch2_induction/RazoMejia_2018.csv', comment='#')
ind_data = ind_data[ind_data['repressors'] > 0]
ind_data.rename(columns={'fold_change_A':'fold_change', 'IPTG_uM':'IPTGuM'}, 
                inplace=True)
ind_data['repressors'] *= 2

# Group and summarize the induction data.
ind_data = ind_data.groupby(['operator', 
                            'repressors', 'IPTGuM'])['fold_change'].agg(
                             ('mean', 'sem')).reset_index()

old_data = pd.read_csv('../../data/other/Garcia2011_Brewster2014.csv', comment='#')
mut_data = pd.read_csv('../../data/ch3_mutants/Chure2019_summarized_data.csv', comment='#')
mut_data  = mut_data[mut_data['repressors'] > 0]

# Load the mutants summary data. 
mutDNA_summary = pd.read_csv('../../data/ch3_mutants/Chure2019_DNA_binding_energy_summary.csv', 
                             comment='#')
mutDNA_summary = mutDNA_summary[(mutDNA_summary['operator']=='O2') & 
                                (mutDNA_summary['repressors']==260)]
mutIND_summary = pd.read_csv('../../data/ch3_mutants/Chure2019_KaKi_epAI_summary.csv', comment='#')
mutIND_summary = mutIND_summary[(mutIND_summary['operator']=='O2') &
                                (mutIND_summary['repressors']==260)]

# Load the entropic parameter inference summary. 
entropy_samps = pd.read_csv('../../data/ch4_growth/pooled_entropic_parameter_samples.csv', comment='#')

#%% Define the parameters for the master curve
bohr_param = np.linspace(-10, 10, 200)
master_curve = (1 + np.exp(-bohr_param))**-1

# Define colors and glyphs
glyphs = {'brewster':'^', 'garcia':'p', 'razomejia':'s', 'chure':'.'}
auth_colors = {'brewster':colors['purple'], 'garcia':colors['orange'], 
                'razomejia':colors['green'], 'chure':colors['blue']}
fills = {'glucose':colors['purple'], 'glycerol':colors['green'],
         'acetate':colors['brown'], 37:colors['purple'],
         32:colors['blue'], 42:colors['red']}

#%% 
# Instantiate the figure
fig, ax = plt.subplots(1, 1, figsize=(3.42, 3.5), dpi=100)
phd.viz.despine(ax)
ax.set_xlabel('free energy [$k_BT$]', fontsize=8)
ax.set_ylabel('fold-change')

ax.plot(bohr_param, master_curve, 'k-', lw=1, label='theoretical prediction', 
        zorder=100)


# ##############################################################################
# BREWSTER AND GARCIA DATA
# ##############################################################################
for g, d in old_data.groupby('author'):
    if g == 'garcia':
        label = 'operator sequence'
        bohr = phd.thermo.SimpleRepression(R=d['repressor'], ep_r=d['energy'],
                                      ka=0.53, ki=139, ep_ai=4.5, 
                                      effector_conc=0).bohr_parameter()
        ax.plot(bohr, d['fold_change'], marker=glyphs[g], alpha=0.4,
                color=auth_colors[g], linestyle='none', ms=5, label=label, zorder=99)


# ##############################################################################
# RAZO MEJIA DATA
# ##############################################################################
ops = np.array([constants[o] for o in ind_data['operator'].values])
bohr = phd.thermo.SimpleRepression(R=ind_data['repressors'], ep_r=ops, 
                                ka=constants['Ka'], ki=constants['Ki'],
                                effector_conc=ind_data['IPTGuM'], 
                                ep_ai=constants['ep_AI']).bohr_parameter()
ax.errorbar(bohr, ind_data['mean'], ind_data['sem'], fmt=glyphs['razomejia'],
            color=auth_colors['razomejia'], ms=3, linestyle='none', alpha=0.4,
            label='inducer concentration')


# ############################################################################## # MUTANTS DATA
# ##############################################################################
mut_data = mut_data[mut_data['mutant'] != 'wt']
labeled = False
for g, d in mut_data.groupby(['mutant', 'class', 'IPTGuM', 'repressors', 'operator']):
    if labeled == False:
        label = 'amino acid sequence'
        labeled = True
    else:
        label = '__nolegend__'
    if g[1] == 'WT':
        ep_R = constants[g[-1]]
        ka = constants['Ka']
        ki = constants['Ki']
        ep_ai = constants['ep_AI']

    elif g[1] == 'DBL':
        DNA, IND = g[0].split('-')
        ep_R = mutDNA_summary[(mutDNA_summary['mutant']==DNA) & 
                              (mutDNA_summary['parameter']=='ep_RA')]['median'].values[0]
        ka = mutIND_summary[(mutIND_summary['mutant']==IND) & 
                              (mutIND_summary['parameter']=='Ka')]['median'].values[0]
        ki = mutIND_summary[(mutIND_summary['mutant']==IND) & 
                              (mutIND_summary['parameter']=='Ki')]['median'].values[0]
        ep_ai = mutIND_summary[(mutIND_summary['mutant']==IND) & 
                              (mutIND_summary['parameter']=='ep_AI')]['median'].values[0]
    elif g[1] == 'DNA':
        ep_R = mutDNA_summary[(mutDNA_summary['mutant']==g[0]) & 
                              (mutDNA_summary['parameter']=='ep_RA')]['median'].values[0]
        ka = constants['Ka']
        ki = constants['Ki']
        ep_ai = constants['ep_AI']

    elif g[1] == 'IND':
        ka = mutIND_summary[(mutIND_summary['mutant']==g[0]) & 
                              (mutIND_summary['parameter']=='Ka')]['median'].values[0]
        ki = mutIND_summary[(mutIND_summary['mutant']==g[0]) & 
                              (mutIND_summary['parameter']=='Ki')]['median'].values[0]
        ep_ai = mutIND_summary[(mutIND_summary['mutant']==g[0]) & 
                              (mutIND_summary['parameter']=='ep_AI')]['median'].values[0]
        ep_R = constants[g[-1]]
        
    bohr = phd.thermo.SimpleRepression(R=g[-2], ep_r=ep_R, 
                                ka=ka, ki=ki, effector_conc=g[2], 
                                ep_ai=ep_ai).bohr_parameter()
    ax.errorbar(bohr, d['mean'], d['sem'], fmt=glyphs['chure'],
            color=auth_colors['chure'], ms=5, linestyle='none', alpha=0.4,
            label=label)


# Add a marker for the conditions. Doing it here for proper ordering
ax.errorbar([], [], [], color='k', markerfacecolor='white', markeredgewidth=0.75,
            linestyle='none', capsize=1, lw=0.75, label='growth conditions', 
            fmt='.', ms=8)

# ##############################################################################
# CARBON DATA
# ##############################################################################
carbon = data[data['temp']==37]
for g, d in carbon.groupby(['carbon']):
    bohr = phd.thermo.SimpleRepression(R=d['repressors']['mean'], ep_r=constants['O2'], 
                                ka=constants['Ka'], ki=constants['Ki'], effector_conc=0, 
                                ep_ai=constants['ep_AI']).bohr_parameter()
    ax.errorbar(bohr, d['fold_change']['mean'], d['fold_change']['sem'], fmt='.',
            color=fills[g], ms=8, linestyle='none', markeredgecolor=colors['grey'],
            label=f'{g}, 37° C', markeredgewidth=0.75, lw=0.75, capsize=1,
            zorder=101)

# ##############################################################################
# TEMP DATA
# ##############################################################################
temp = data[data['temp']!=37]
for g, d in temp.groupby(['temp']):
    epRA_star = entropy_samps[(entropy_samps['temp']==g) & 
                            (entropy_samps['parameter']=='epRA_star')]['value'].median()
    epAI_star = entropy_samps[(entropy_samps['temp']==g) & 
                            (entropy_samps['parameter']=='epAI_star')]['value'].median()
    bohr = phd.thermo.SimpleRepression(R=d['repressors']['mean'], ep_r=epRA_star, 
                                ka=constants['Ka'], ki=constants['Ki'], effector_conc=0, 
                                ep_ai=epAI_star).bohr_parameter()
    ax.errorbar(bohr, d['fold_change']['mean'], d['fold_change']['sem'], fmt='.',
            color=fills[g], ms=8, linestyle='none', markeredgecolor=colors['grey'],
            label=f'glucose, {int(g)}° C', markeredgewidth=0.75, lw=0.75, 
            capsize=1, zorder=101)

# Plot the master curve. 
ax.legend()
plt.savefig('../figs/ch4_fig6.pdf', bbox_inches='tight')
plt.savefig('../figs/ch4_fig6.png', bbox_inches='tight')


# %%
