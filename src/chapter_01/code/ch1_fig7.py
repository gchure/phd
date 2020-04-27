#%%
import numpy as np
import matplotlib.pyplot as plt
import phd.viz
colors, palette = phd.viz.phd_style()


# %%
# Define the master curve
bohr_range = np.linspace(-6, 6, 200)
collapse_curve = (1 + np.exp(-bohr_range))**-1

# Define the points in free energy space
f_points = np.array([-1, 3, -4, 0])
fc_points = (1 + np.exp(-f_points))**-1
point_colors = [colors['light_grey'], colors['purple'], colors['orange'], colors['blue']]

# Set up the figure canvas
fig, ax = plt.subplots(1, 1, figsize=(3, 2))
phd.viz.despine(ax)


# Format the axes
ax.set_xlabel('free energy [$k_BT$]')
ax.set_ylabel('fold-change')
ax.set_ylim([-0.1, 1.1])


ax.plot(bohr_range, collapse_curve, '-', color=colors['black'], lw=1)
for i, c in enumerate(point_colors):
    ax.plot(f_points[i], fc_points[i], 'o', color=c, ms=7, 
            markeredgecolor=colors['grey'], markeredgewidth=0.75, zorder=100)
# Add the steps
ax.step(f_points[:2], fc_points[:2], '-', color=colors['light_purple'],     
        where='post')
ax.step([f_points[0], f_points[2]], [fc_points[0], fc_points[2]],
        '-', color=colors['orange'], where='post')
ax.step([f_points[0], f_points[3]], [fc_points[0], fc_points[3]], '-',
        color=colors['blue'])
# Add the text labels
ax.text(-0.9, 0.18, '$F_{wt}$', fontsize=6, color=colors['black'])
ax.text(-3.8, -0.05, '$F_1$', fontsize=6, color=colors['orange'])
ax.text(-4.2, 0.30, r'$\Delta F_1 = F_1 - F_{wt}$', color=colors['orange'], fontsize=6)
ax.text(3.2, 0.89, r'$F_2$', color=colors['purple'], fontsize=6)
ax.text(-0.5, 0.3, r'$\Delta F_2 = F_2 - F_{wt}$', color=colors['purple'], fontsize=6)
ax.text(0.52, 0.48,  r'$F_{1,2}$', color=colors['blue'], fontsize=6, zorder=1000)

ax.text(-4, 0.55,  r'$\Delta F_{1,2} = \Delta F_1 + \Delta F_2$', 
        color=colors['blue'], fontsize=6) #, backgroundcolor=colors['grey'])
plt.savefig('../figs/ch1_fig6_plot.svg', bbox_inches='tight')
# %%


# %%
