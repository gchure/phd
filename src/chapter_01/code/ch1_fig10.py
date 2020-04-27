#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import phd.viz
colors, palette = phd.viz.phd_style()

# Load the growth curve data.
plates = pd.read_csv('../../data/ch4_growth/compiled_growth_plates.csv', comment='#')
stats = pd.read_csv('../../data/ch4_growth/compiled_growth_statistics.csv', comment="#")
plates['temp_C'] = np.round(plates['temp_C'])
plates['time_min'] = np.round(plates['time_min'])
plates = plates[plates['time_min'] <= 600]

# %%
# reform the stats 
stats = stats[((stats['carbon']=='glucose') | (stats['carbon']=='acetate') |
                (stats['carbon']=='glycerol')) & 
                ((stats['temp']==37) |  (stats['temp']==32) | 
                (stats['temp']==42))] 

tidy_stats = pd.DataFrame([])
for g, d in stats.groupby(['date', 'carbon', 'temp', 'run_number']):
    growth_rate = d[d['parameter']=='max df']['value'].values[0]
    growth_err = d[d['parameter']=='max df std']['value'].values[0]
    dbl_time = d[d['parameter']=='inverse max df']['value'].values[0]
    dbl_err = d[d['parameter']=='inverse max df std']['value'].values[0]
   

    tidy_stats = tidy_stats.append({'date':g[0], 'carbon':g[1], 'temp_C':g[2], 'run_number':g[3],
                                    'growth_rate':growth_rate, 
                                    'dbl_time':dbl_time,
                                    'growth_err':growth_err,
                                    'dbl_err':dbl_err}, 
                                    ignore_index=True)
tidy_stats['growth_rate'] *= 60
tidy_stats['growth_err'] *= 60
#%%
# Restrict only to conditions studied here.
plates = plates[((plates['carbon']=='LB') | (plates['carbon']=='glucose') | (plates['carbon']=='acetate') |
                (plates['carbon']=='glycerol') | (plates['carbon']=='blank')) & 
                ((plates['temp_C']==37) |  (plates['temp_C']==32) | 
                (plates['temp_C']==42))] 

#%% Subtract the blank from each date, run, and time point.
sub_plates = []
for g, d in plates.groupby(['date', 'run_number', 'time_min']):
    d = d.copy()
    blank = d[d['carbon']=='blank']['od_600nm'].mean()
    d['od_sub'] = d['od_600nm'] - blank
    sub_plates.append(d)
sub_plates = pd.concat(sub_plates)
_sub_plates = []
for g, d in sub_plates.groupby(['date', 'carbon', 'temp_C']):
    d = d.copy()
    d = d[d['od_sub']>0]
    a0 = d[d['time_min']==6]['od_sub'].mean()
    d['rel_od'] = d['od_sub'] / a0
    _sub_plates.append(d)
sub_plates = pd.concat(_sub_plates)
#%%
# Choose representable samples. 
#Define the colors
fill_colors = {'acetate': colors['brown'], 'glycerol': colors['green'],
               'glucose':colors['purple'], 37: colors['purple'],
              32:colors['blue'], 42:colors['red']}

fig, ax = plt.subplots(2, 1, figsize=(3, 3.5), dpi=300)
phd.viz.despine(ax)

labels = ['acetate, 37°C', 'glycerol, 37°C','glucose, 37°C', 'glucose, 32°C', 'glucose, 42°C']
carbs = ['acetate', 'glycerol', 'glucose', 'glucose', 'glucose']
vals = ['acetate', 'glycerol', 'glucose', 32,  42]
temps = [37, 37, 37, 32, 42]
for c, t, l, v in zip(carbs, temps, labels, vals):
    for g, d in sub_plates[(sub_plates['carbon']==c) & (sub_plates['temp_C']==t)].groupby('carbon'):
        grp = d.groupby('time_min')['rel_od'].agg(('mean', 'std')).reset_index()
        _stats = tidy_stats[(tidy_stats['carbon']==c) & (tidy_stats['temp_C']==t)]
        if (t == 37) & (c != 'glucose'):
            axes = [ax[0]]
        if (t!=37):
            axes = [ax[1]]
        if (t==37) & (c == 'glucose'):
            axes = ax
        for a in axes:
            a.plot(grp['time_min'], grp['mean'], '-o', ms=2.5,
               lw=.75, markerfacecolor=fill_colors[v],
               color=fill_colors[v], markeredgewidth=0.4, markeredgecolor=colors['grey'],
               label=f'{c}, {t}° C' + '\n' + r'$t_{double}=$' + f"{int(_stats['dbl_time'].mean())} ± {int(_stats['dbl_time'].std())} min")

for a in ax:
    a.legend(handlelength=1, fontsize=5.5)
    a.set_yscale('log')
    a.set_xlim([0, 400])
    a.set_xlabel('time [min]')
    a.set_ylabel('relative optical density')
plt.tight_layout()
plt.savefig('../figs/ch1_fig10_growth_curves.svg', bbox_inches='tight')

# %%
