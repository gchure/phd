#%%
import numpy as np
import matplotlib.pyplot as plt
import phd.viz 
import phd.thermo
colors, palette = phd.viz.phd_style()
constants = phd.thermo.load_constants()

# Define the parameter ranges. 
ep_range = np.linspace(-18, -8, 200)
rep_colors = {22: colors['red'], 
              60: colors['brown'],  
              124: colors['green'], 
              260: colors['orange'], 
              1220: colors['purple'], 
              1740: colors['blue']}

#%% Set up the figure canvas.
fig, ax = plt.subplots(1, 3, figsize=(6, 2))
phd.viz.despine(ax)

# Format the axes. 
titles = ['leakiness', 'saturation', 'dynamic range']
for i, a in enumerate(ax):
    a.set_xlabel('DNA binding energy [$k_BT$]')
    a.set_ylabel(titles[i])

for r, c in rep_colors.items():
    props = phd.thermo.SimpleRepression(R=r, ep_r=ep_range, 
                                       ka=constants['Ka'], ki=constants['Ki'],
                                       ep_ai=constants['ep_AI'], effector_conc=0
                                       ).compute_properties()

    ax[0].plot(ep_range, props['leakiness'], color=c, label=r, lw=0.75) 
    ax[1].plot(ep_range, props['saturation'], color=c, lw=0.75)
    ax[2].plot(ep_range, props['dynamic_range'], color=c, lw=0.75)


# Add legend and panel labels. 
leg = ax[0].legend(title='repressors\n   per cell', fontsize=6, handlelength=1)
leg.get_title().set_fontsize(6)
plt.tight_layout()
fig.text(0.02, 0.95, '(A)', fontsize=8)
fig.text(0.35, 0.95, '(B)', fontsize=8)
fig.text(0.68, 0.95, '(C)', fontsize=8)
plt.savefig('../figs/ch6_figS18.pdf', bbox_inches='tight')
plt.savefig('../figs/ch6_figS18.png', bbox_inches='tight')
# %

# Make an analogous figure for the effecitve hill and EC50.

#%% Set up the figure canvas.
fig, ax = plt.subplots(1, 2, figsize=(5, 2))
phd.viz.despine(ax)

# Format the axes. 
titles = ['EC$_{50}$ [µM]', 'effective Hill coefficient']
for i, a in enumerate(ax):
    a.set_xlabel('DNA binding energy [$k_BT$]')
    a.set_ylabel(titles[i])
ax[0].set_yscale('log')

for r, c in rep_colors.items():
    props = phd.thermo.SimpleRepression(R=r, ep_r=ep_range, 
                                       ka=constants['Ka'], ki=constants['Ki'],
                                       ep_ai=constants['ep_AI'], effector_conc=0
                                       ).compute_properties()

    ax[0].plot(ep_range, props['EC50'], color=c, label=r, lw=0.75) 
    ax[1].plot(ep_range, props['effective_hill'], color=c, lw=0.75)



# Add legend and panel labels. 
leg = ax[0].legend(title='repressors\n   per cell', fontsize=6, handlelength=1)
leg.get_title().set_fontsize(6)
plt.tight_layout()
fig.text(0.02, 0.95, '(A)', fontsize=8)
fig.text(0.5, 0.95, '(B)', fontsize=8)
plt.savefig('../figs/ch6_figS19.pdf', bbox_inches='tight')
plt.savefig('../figs/ch6_figS19.png', bbox_inches='tight')
#
# %%
