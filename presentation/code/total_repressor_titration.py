#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import phd.viz
import phd.thermo
colors, palette = phd.viz.phd_style()

# %%
# Define the number of repressors 
rep_range = np.logspace(0, 3, 200)

# Define the DNA binding energies. 
ep_range = np.array([-16, -14, -12, -10])

# Mesh and compute the fold-change. 
c_range = np.logspace(-2, 4, 200)
r, ep, c = np.meshgrid(rep_range, ep_range, c_range)
arch = phd.thermo.SimpleRepression(R=r, ep_r=ep, ep_ai=5, effector_conc=c,
                                ka=200, ki=1)
fc = arch.fold_change()

# %%
# Set up the figure canvas.
fig, ax = plt.subplots(2, 1, figsize=(4, 5))
phd.viz.despine(ax)
for a in ax:
    a.set_xscale('log')
    a.set_ylabel('fold-change')
ax[0].set_xlim([1, 1E3])
ax[1].set_xlim([1E-2/200, 1E4/200])
ax[0].set_yscale('log')
ax[0].set_xlabel('total repressors per cell')
ax[1].set_xlabel('inducer concentration relative to $K_A$ ($c / K_A$)')

# Plot the foldchange curves. 
ep_colors = [colors['dark_red'], colors['dark_orange'], colors['orange'], 
            colors['light_orange']]

for i in range(4):
    ax[0].plot(rep_range, fc[i, :, 0], lw=1, color=ep_colors[i], label=int(ep_range[i]))
    ax[1].plot(c_range/200, fc[i, 135, :], lw=1, color=ep_colors[i], label=int(ep_range[i]))

leg = ax[0].legend(title=r'$\Delta\varepsilon_{RA}$ [$k_BT$]')
leg.get_title().set_fontsize(6)

plt.subplots_adjust(hspace=0.5)
plt.savefig('../figs/total_repressor_titration.svg', bbox_inches='tight') 

#

# %%
