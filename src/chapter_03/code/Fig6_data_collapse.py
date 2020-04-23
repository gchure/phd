# -*- coding: utf-8 -*-
#%%
import pandas as pd
import numpy as np
import phd.thermo
import matplotlib.pyplot as plt
import phd.viz
import phd.thermo
colors, palette = phd.viz.phd_style()
constants = phd.thermo.load_constants()

# Add DNA binding energy for Oid, not queried in this work
constants['Oid'] = -17
mut_colors = phd.viz.color_selector('mut')

# Load the data from the mutants work.
data = pd.read_csv('../../data/ch3_mutants/Chure2019_summarized_data.csv', comment='#')
data = data[data['mutant']!='wt']
epRA_stats = pd.read_csv('../../data/ch3_mutants/Chure2019_DNA_binding_energy_summary.csv')
epRA_stats = epRA_stats[epRA_stats['repressors']==260]
allo_stats = pd.read_csv('../../data/ch3_mutants/Chure2019_KaKi_epAI_summary.csv')
allo_stats = allo_stats[allo_stats['operator']=='O2']

# Load the data from the old gods
old_gods = pd.read_csv('../../data/other/Garcia2011_Brewster2014.csv', comment='#')
new_gods = pd.read_csv('../../data/ch2_induction/RazoMejia_2018.csv', comment='#')
new_gods['repressors'] *= 2
new_gods.rename(columns={'IPTG_uM':'IPTGuM', 'fold_change_A':'fold_change'},
    inplace=True)
new_gods = new_gods[new_gods['repressors'] > 0]

# Define plotting constants
bohr_range = np.linspace(-10, 10, 200)

# Define colors and glyphs
glyphs = {'garcia':'>', 'brewster':'s'}
_color = {'garcia': colors['red'], 'brewster':colors['blue']}
legend = {'garcia': 'Garcia & Phillips 2011', 'brewster':'Brewster et al. 2014'}

# ##############################################################################
# FIGURE INSTANTIATION
# ##############################################################################
#%%
fig, ax = plt.subplots(1, 1, figsize=(3.42, 3))
phd.viz.despine(ax)
ax.xaxis.set_tick_params(labelsize=6)
ax.yaxis.set_tick_params(labelsize=6)
ax.set_ylim([-0.05, 1.25])
ax.set_xlabel('free energy [$k_BT$]', fontsize=8)
ax.set_ylabel('fold-change', fontsize=8)

# ##############################################################################
# COLLAPSE CURVE
# ##########################3333333#3333########################################
ax.plot(bohr_range, (1 + np.exp(-bohr_range))**-1, 'k-', lw=1)

# ##############################################################################
# GARCIA AND BREWSTER DATA
# ##############################################################################
for g, d in old_gods.groupby('author'):
    ops = [constants[o] for o in d['operator'].values]
    bohr = phd.thermo.SimpleRepression(R=d['repressor'], ep_r=ops,
                                       ka=constants['Ka'], ki=constants['Ki'],
                                       ep_ai=constants['ep_AI'], 
                                       effector_conc=0).bohr_parameter()
    ax.plot(bohr, d['fold_change'], marker=glyphs[g], color=_color[g], 
            linestyle='none', label=legend[g], markeredgewidth=0.75, alpha=0.75, 
            ms=3)

# ##############################################################################
# RAZO-MEJIA ET AL 2018
# ##############################################################################
new_gods = new_gods.groupby(['operator', 'repressors', 'IPTGuM']).agg(('mean', 'sem')).reset_index()
ops = [constants[o] for o in new_gods['operator'].values]
bohr = phd.thermo.SimpleRepression(R=new_gods['repressors'], ep_r=ops,
                                   ka=constants['Ka'], ki=constants['Ki'],
                                   ep_ai=constants['ep_AI'], 
                             effector_conc=new_gods['IPTGuM']).bohr_parameter() 

ax.errorbar(bohr, new_gods['fold_change']['mean'], 
            new_gods['fold_change']['sem'], fmt='o',  color=colors['green'],
                    markeredgewidth=0.75, alpha=0.5,
                    linestyle='none', lw=0.75, capsize=1, label='Razo-Mejia et al. 2018',
                    ms=3)

# ##############################################################################
# DNA MUTS
# ##############################################################################
op_glyphs = {'O1':'^', 'O2':'v', 'O3':'D'}
for g, d in data[data['class']=='DNA'].groupby(['mutant']):
    if g == 'Q21A':
        _g = 'Q18A'
    elif g == 'Q21M':
        _g = 'Q18M'
    elif g == 'Y20I':
        _g = 'Y17I'
    ep_RA = epRA_stats[(epRA_stats['mutant']==g) &
                (epRA_stats['parameter']=='ep_RA')]['median'].values[0]
    bohr = phd.thermo.SimpleRepression(R=d['repressors'], ep_r=ep_RA, 
                                           ka=constants['Ka'], ki=constants['Ki'],
                                           ep_ai=constants['ep_AI'],
                                           effector_conc=d['IPTGuM']).bohr_parameter()
    ax.errorbar(bohr, d['mean'], d['sem'], fmt='^', color=mut_colors[g],
                label=_g, ms=5, markeredgewidth=0.75, markerfacecolor='w',
                lw=0.75)

# ##############################################################################
#  IND MUTS
# ##############################################################################
for g, d in data[data['class']=='IND'].groupby(['mutant']):
    if g == 'Q294K':
        _g = 'Q291K'
    elif g == 'Q294R':
        _g = 'Q291K'
    elif g == 'Q294V':
        _g = 'Q291V'
    elif g == 'F164T':
        _g = 'F161T'
    _stats = allo_stats[(allo_stats['mutant']==g)]
    ka = _stats[_stats['parameter']=='Ka']['median'].values[0]
    ki = _stats[_stats['parameter']=='Ki']['median'].values[0]
    ep_AI = _stats[_stats['parameter']=='ep_AI']['median'].values[0]
    ops = [constants[o] for o in d['operator'].values]
    bohr = phd.thermo.SimpleRepression(R=d['repressors'], ep_r=ops, 
                                           ka=ka, ki=ki,
                                           ep_ai=ep_AI,
                                           effector_conc=d['IPTGuM']).bohr_parameter()
    ax.errorbar(bohr, d['mean'], d['sem'], fmt='p', color=mut_colors[g],
                label=g, ms=5, markeredgewidth=1,  markerfacecolor='w',
                lw=0.75)

# ##############################################################################
# DBL MUTS
# ##############################################################################
for g, d in data[data['class']=='DBL'].groupby(['mutant']):
    _dna, _ind =  g.split('-')
    orig = _dna[0]
    end= _dna[-1]
    pos = int(_dna[1:-1])-3
    new_dna = orig + str(pos) + end
    orig = _ind[0]
    end= _ind[-1]
    pos = int(_ind[1:-1])-3
    new_ind = orig + str(pos) + end
    mut_name = f'{new_dna}-{new_ind}'
    ep_RA = epRA_stats[(epRA_stats['mutant']==g.split('-')[0]) & 
                           (epRA_stats['parameter']=='ep_RA')]['median'].values[0]
    _stats = allo_stats[(allo_stats['mutant']==g.split('-')[1])]
    ka = _stats[_stats['parameter']=='Ka']['median'].values[0]
    ki = _stats[_stats['parameter']=='Ki']['median'].values[0]
    ep_AI = _stats[_stats['parameter']=='ep_AI']['median'].values[0]
    bohr = phd.thermo.SimpleRepression(R=d['repressors'], ep_r=ep_RA, 
                                           ka=ka, ki=ki,
                                           ep_ai=ep_AI,
                                           effector_conc=d['IPTGuM']).bohr_parameter()
    ax.errorbar(bohr, d['mean'], d['sem'], fmt='v', color=mut_colors[g],
                label=mut_name, ms=4, markeredgewidth=0.75,  markerfacecolor='w',
                 lw=0.75)

# ##############################################################################
# LEGEND DETAILS AND FIGURE SAVING
# ##############################################################################
ax.legend(loc='upper left', fontsize=5)
plt.tight_layout()
plt.savefig('../figs/ch3_fig6.pdf', bbox_inches='tight')
plt.savefig('../figs/ch3_fig6.png', bbox_inches='tight')

# %%
