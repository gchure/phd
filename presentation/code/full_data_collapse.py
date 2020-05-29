#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import phd.thermo
import phd.viz
colors, palette = phd.viz.phd_style()
INDUCTION = False
MUTS = False
GROWTH = False

# Load the compiled datasets
data = pd.read_csv('../../src/data/other/Garcia2011_Brewster2014_RazoMejia2018_Chure2019.csv', 
                   comment='#')
growth_data = pd.read_csv('../../src/data/ch4_growth/analyzed_foldchange.csv', comment="#")
entropy_samples = pd.read_csv('../../src/data/ch4_growth/pooled_entropic_parameter_samples.csv', comment='#')

# Summarize the growth data. 
growth_data = growth_data.groupby(['atc_ngml', 'date', 'run_number', 
                                    'carbon', 'temp'])[
                                    ['repressors', 'fold_change']].mean().reset_index()
growth_data = growth_data.groupby(['atc_ngml', 'carbon', 'temp']).agg(('mean', 'sem')).reset_index()
growth_data

# Compute the bohr parameter
dfs = []
for g, d in growth_data.groupby(['carbon', 'temp']):
    if g[1] != 37:
        _samples = entropy_samples[(entropy_samples['temp']==g[1])]
        epRA = _samples[_samples['parameter']=='epRA_star']['value'].mean()
        epAI = _samples[_samples['parameter']=='epAI_star']['value'].mean()
    else:
        epRA = -13.9
        epAI = 4.5
    
    bohr = phd.thermo.SimpleRepression(d['repressors']['mean'] * 1.16, epRA, ka=139,
                           ki=0.53, ep_ai=epAI, effector_conc=0).bohr_parameter() 
    d['bohr_parameter'] = bohr
    dfs.append(d)
growth = pd.concat(dfs, sort=False)

# Set up the figure canvas
fig, ax = plt.subplots(1, 1, figsize=(5, 5))
phd.viz.despine(ax)
bohr_range = np.linspace(-12, 12, 200)
scaling_fn = (1 + np.exp(-bohr_range))**-1
_ = ax.plot(bohr_range, scaling_fn, 'k-', label='scaling function')
ax.set_xlabel('free energy [$k_BT$]')
ax.set_ylabel('fold-change')
ax.set_ylim([-0.05, 1.25])
ax.set_xlim([-12, 12])

induction = data[data['author']=='razo-mejia'] 

rep_colors = {22: colors['red'], 
              60: colors['brown'],  
              124: colors['green'], 
              260: colors['orange'], 
              1220: colors['purple'], 
              1740: colors['blue']}
op_glyphs = {'O1':'^', 'O2':'o', 'O3':'v'}
class_glyphs = {'DNA': 'h', 'IND':'p', 'DBL':'X'}
mut_colors = {'Y20I':colors['red'], 'Q21A':colors['purple'], 'Q21M':colors['light_purple'],
              'F164T':colors['green'], 'Q294R':colors['brown'], 'Q294V':colors['orange'], 
              'Q294K':colors['light_grey'], 'Y20I-F164T':colors['black'], 
              'Y20I-Q294V':colors['dark_brown'], 'Y20I-Q294K':colors['dark_red'],
              'Q21A-F164T':colors['dark_blue'], 'Q21A-Q294V':colors['dark_purple'],
              'Q21A-Q294K':colors['dark_green'], 'Q21M-Q294V':colors['orange'],
              'Q21M-F164T':colors['pale_green'], 'Q21M-Q294K':colors['light_red'],
              'Q21M-Q294V':colors['light_purple']}

temp_colors = {32:colors['blue'], 42:colors['red']}
carbon_colors = {'glycerol':colors['green'], 'acetate':colors['brown']}

name = 'scaling_fn'
if INDUCTION:
    name += '_INDUCTION'
    for o, g in op_glyphs.items():
        _ = ax.plot([], [], f'{g}', color=colors['black'], 
                    markeredgecolor=colors['grey'], markeredgewidth=0.5,
                    label=f'operator {o}', ms=4)
    for r, c in rep_colors.items():
        _ = ax.plot([], [], '-', color=c, label=f'{int(r)} repressors',
        ms=4)
    for g, d in induction.groupby(['repressors', 'operator']):
        marker = op_glyphs[g[1]]
        _color = rep_colors[g[0]]
        _ = ax.errorbar(d['bohr_parameter'], d['mean'], d['sem'], color=_color,
                        markeredgecolor=colors['grey'], fmt=marker, linestyle='none', 
                        ms=4, markeredgewidth=0.5, label='__nolegend__')

DNA_muts = ['Q21M', 'Y20I', 'Q21A']
IND_muts =['F64T', 'Q294K', 'Q294V']
if MUTS:
    name += '_MUTS'
    for c, g in class_glyphs.items():
        if c == 'DNA':
            label = 'DNA binding mutant'
            _color = colors['purple']
        elif c == 'IND':
            label = 'inducer binding mutant'
            _color = colors['green']
        else:
             label = 'double mutant'
             _color = colors['blue']

        _ = ax.plot([], [], f'{g}', color=_color, markeredgecolor=colors['grey'],
                    markeredgewidth=0.5, label=label, ms=4)
    for g, d in data[(data['author']=='chure') & (data['mutant']!='wt')].groupby(['mutant']):
        if (g in DNA_muts) & ('-' not in g):
            marker = class_glyphs['DNA']
            _color = colors['purple']
        if (g in IND_muts) & ('-' not in g):
            marker = class_glyphs['IND']
            _color = colors['green']
        elif '-' in g:
            marker = class_glyphs['DBL']
            _color = colors['blue']

        _ = ax.errorbar(d['bohr_parameter'], d['mean'], d['sem'], fmt=marker, ms=4,
                    color=_color, markeredgewidth=0.5, markeredgecolor=colors['grey'],
                    linestyle='none', label='__nolegend__') 

if GROWTH:
    name += '_GROWTH'
    _ = ax.plot([], [], '>', color=colors['red'], label='grown at 42°C', ms=4, 
                markeredgewidth=0.5, markeredgecolor=colors['grey'])
    _ = ax.plot([], [], '>', color=colors['blue'], label='grown at 32°C', ms=4, 
                markeredgewidth=0.5, markeredgecolor=colors['grey'])
    _ = ax.plot([], [], '<', color=colors['green'], label='grown in glycerol', ms=4, 
                markeredgewidth=0.5, markeredgecolor=colors['grey'])
    _ = ax.plot([], [], '<', color=colors['brown'], label='grown in acetate', ms=4, 
                markeredgewidth=0.5, markeredgecolor=colors['grey'])
    for g, d in growth.groupby(['temp', 'carbon']):
        if (g[0] != 37) | (g[1]!='glucose'):
            if g[0] != 37:
                marker = '>'
                _color = temp_colors[g[0]]
            else:
                marker = '<'
                _color = carbon_colors[g[1]]
            _ = ax.errorbar(d['bohr_parameter'], d['fold_change']['mean'], 
                        d['fold_change']['sem'], fmt=marker, ms=4,
                        color=_color, linestyle='none', markeredgewidth=0.5, 
                        markeredgecolor=colors['grey'], label='__nolegend__',
                        zorder=1000)
ax.legend(loc='upper left')


plt.savefig(f'../figs/collapse_{name}.pdf')

# %%



# %%
