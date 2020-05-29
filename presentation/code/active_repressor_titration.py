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
r, ep = np.meshgrid(rep_range, ep_range)
arch = phd.thermo.SimpleRepression(R=r, ep_r=ep, ep_ai=100, effector_conc=0,
                                ka=200, ki=0.5)
fc = arch.fold_change()

# %%
# Set up the figure canvas.
fig, ax = plt.subplots(1, 1, figsize=(4, 3))
phd.viz.despine(ax)

# Format the axes
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('active repressor copy number')
ax.set_ylabel('fold-change in gene expression')
ax.set_xlim([1, 1E3])

# Plot the foldchange curves. 
ep_colors = [colors['red'], colors['dark_orange'], colors['orange'], 
            colors['light_orange']]

for i in range(4):
    ax.plot(rep_range, fc[i, :], lw=1, color=ep_colors[i], label=int(ep_range[i]))
l = ax.legend(title=r'$\Delta\varepsilon_{RA}$ [$k_BT$]')
l.get_title().set_fontsize(6)

plt.savefig('../figs/active_repressor_titration.pdf', bbox_inches='tight') 



# %%
