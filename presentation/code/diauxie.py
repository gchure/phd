#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import phd.viz
import scipy.interpolate
colors, palette = phd.viz.phd_style()

data = pd.read_csv('../../src/data/ch1_introduction/Monod1941_Fig1_Fig2.csv',
                    comment='#')


glucose = False
arabinose = False

_glucose = data[data['secondary_sugar']=='glucose']
_arabinose = data[data['secondary_sugar']=='arabinose']

fig, ax = plt.subplots(1, 1, figsize=(4, 3))
phd.viz.despine(ax)

name = 'canvas'
if glucose:
    name += '_glucose'
    spline = scipy.interpolate.UnivariateSpline(_glucose['hours'], _glucose['optical_density'])    
    spline.set_smoothing_factor(30)
    time = np.linspace(0, _glucose['hours'].max(), 100)
    ax.plot(time, spline(time), '-', color=colors['blue'], lw=0.75, label='__nolegend__')
    ax.plot(_glucose['hours'], _glucose['optical_density'], 'o', color=colors['blue'],
        markeredgecolor=colors['grey'], markeredgewidth=0.5, label=f'sucrose + glucose',
        ms=5, linewidth=0.5)

if arabinose:
    name += '_arabinose'
    spline = scipy.interpolate.UnivariateSpline(_arabinose['hours'], _arabinose['optical_density'])    
    spline.set_smoothing_factor(10)
    time = np.linspace(0, _arabinose['hours'].max(), 100)
    ax.plot(time, spline(time), '-', color=colors['green'], lw=0.75, label='__nolegend__')
    ax.plot(_arabinose['hours'], _arabinose['optical_density'], 'o', color=colors['green'],
        markeredgecolor=colors['grey'], markeredgewidth=0.5, label=f'sucrose + arabinose',
        ms=5, linewidth=0.5)
ax.set_ylim([0, 90])
ax.set_xlim([0, 9.5])
ax.set_ylabel('optical density')
ax.set_xlabel('time [hours]')
ax.legend(loc='upper left')
plt.savefig(f'../figs/monod1941_{name}.pdf', bbox_inches='tight')
# %%
