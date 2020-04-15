
#%%
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import phd.viz
colors, color_list = phd.viz.phd_style()

# Load specific examples of GP processing. 
gp_output = pd.read_csv('../../data/ch8_growth_si/glycerol_example_gp_output.csv', comment='#')
plate = pd.read_csv('../../data/ch8_growth_si/example_growth_plate.csv', comment='#')
glyc = plate[(plate['carbon']=='glycerol') & (plate['time_min']<1200)]

# %%
# Instantiate figure
fig, ax = plt.subplots(1, 3, figsize=(6, 2.5))
phd.viz.despine(ax)
# Add axis labels
for i in range(2):
        ax[i].set_xlabel('time [min]')
        ax[i].set_ylabel('optical density (600 nm)')
ax[0].set_ylim([0, 0.45])
ax[2].set_xlabel('time [min]')
ax[2].set_ylabel('growth rate [hr$^{-1}$]')

# Plot the growth curves and growth rates
ax[0].plot(glyc['time_min'], glyc['od_sub'], '.', color=colors['green'],
           ms=1, alpha=0.15, label='__nolegend__')
ax[0].fill_betweenx(np.linspace(0, 0.6, 100), gp_output['time'].min(), gp_output['time'].max(), 
                    color=colors['light_grey'], alpha=0.25, label='exponential region')
ax[1].plot(gp_output['time'], gp_output['OD_raw_data'], '.', 
        color=colors['green'], 
        markeredgecolor='white', alpha=0.5, markeredgewidth=0.1, 
        label='__nolegend__', ms=3)
ax[1].plot(gp_output['time'], 
        np.exp(gp_output['log(OD)_fit']), '-', color=colors['green'], lw=0.5,
        label='estimated value')
ax[2].fill_between(gp_output['time'],
            (gp_output['growth_rate'] - gp_output['growth_rate_std']) * 60, 
            (gp_output['growth_rate'] + gp_output['growth_rate_std']) * 60,
            color=colors['light_green'], label='standard deviation')
ax[2].plot(gp_output['time'],
            gp_output['growth_rate'] * 60, '-', color=colors['green'], 
            label='mean')
# Add the legends
ax[0].legend()
ax[1].legend()       
ax[2].legend()

# Add the titleboxes
phd.viz.titlebox(ax[0], 'growth curve', size=8, bgcolor='white',
                boxsize="12%", color=colors['black'])
phd.viz.titlebox(ax[1], 'exponential region', size=8, bgcolor='white',
                boxsize="12%", color=colors['black'])
phd.viz.titlebox(ax[2], 'growth rate', size=8, bgcolor='white',
                boxsize="12%", color=colors['black'])
plt.tight_layout()

# Add panel labels
fig.text(0.009, 0.9, '(A)', fontsize=8)
fig.text(0.35, 0.9, '(B)', fontsize=8)
fig.text(0.67, 0.9, '(C)', fontsize=8)

# Add titles
plt.savefig('../figs/ch8_figS1.pdf', bbox_inches='tight')
plt.savefig('../figs/ch8_figS1.png', bbox_inches='tight')
