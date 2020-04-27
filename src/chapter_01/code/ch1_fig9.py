#%% 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import phd.viz
colors, palette = phd.viz.phd_style()

# Load the data
triauxie = pd.read_csv('../../data/ch1_introduction/Monod1947_Fig11.csv')

# Manually label the points by their growth phase. 
triauxie['carbon'] = 'adaptation'
triauxie.loc[(triauxie['time_hrs'] < 2.5), 'carbon'] = 'glucose'
triauxie.loc[(triauxie['time_hrs'] >= 3) & (triauxie['time_hrs'] < 3.8), 'carbon'] = 'sorbitol'
triauxie.loc[(triauxie['time_hrs'] >= 4), 'carbon'] = 'glycerol'

# Compute the splines
spline = scipy.interpolate.UnivariateSpline(triauxie['time_hrs'], triauxie['optical_density'])
spline.set_smoothing_factor(0.5)

# Instantiate the figure canvas.
fig, ax = plt.subplots(1,1, figsize=(3, 2.5))
phd.viz.despine(ax)
ax.set_xlabel('hours')
ax.set_ylabel('optical density')
ax.set_xlim([0, 5.75])
ax.set_ylim([-0.1, 79])
phd.viz.titlebox(ax, 'Monod 1947', color=colors['black'], boxsize='12%', size=6, 
                pad=0.05)

# Define the colors for the carbon points
carbon_colors = {'adaptation': colors['black'], 'glucose':colors['blue'], 
                'sorbitol': colors['orange'], 'glycerol':colors['green']}

for g, d in triauxie.groupby('carbon'):
    if d['time_hrs'].max() == triauxie['time_hrs'].max():
        final = d['time_hrs'].max()
    else:
        final = d['time_hrs'].max() + 0.1
    time_range = np.linspace(d['time_hrs'].min() - 0.3, final, 100)
    ax.plot(time_range, spline(time_range), '-', lw=0.75, color=carbon_colors[g])

    ax.plot(d['time_hrs'], d['optical_density'], 'o', ms=4.5, 
            markeredgecolor=colors['grey'], markeredgewidth=0.5, 
            color=carbon_colors[g], zorder=1000)

# Add textual labels. 
ax.text(1, 20, '  glucose\nutilization', fontsize=6, color=colors['blue'])
ax.text(3.5, 38, '  sorbitol\nutilization', fontsize=6, color=colors['orange'])
ax.text(4.5, 58, '  glycerol\nutilization', fontsize=6, color=colors['green'])

# Add stripes for adaptation periods
ax.vlines(2.75, -0.1, 80, lw=15, color='white', alpha=0.75)
ax.vlines(3.97, -0.1, 80, lw=5.5, color='white', alpha=0.75)
plt.savefig('../figs/ch1_fig9.pdf', bbox_inches='tight')
plt.savefig('../figs/ch1_fig9.png', bbox_inches='tight')

# %%
