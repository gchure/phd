#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import phd.viz
import phd.thermo
colors, palette = phd.viz.phd_style()

# Define a concentration range
c_range = np.logspace(-2, 4, 200)
pact = phd.thermo.MWC(ka=200, ki=1, ep_ai=5, effector_conc=c_range).pact()

# %%
fig, ax = plt.subplots(1, 1, figsize=(3, 2))
phd.viz.despine(ax)
ax.set_xscale('log')
ax.plot(c_range, pact, '-', color=colors['blue'], lw=1)
ax.set_xlabel('inducer concentration')
ax.set_ylabel('probability of being active')
ax.set_ylim([0, 1.1])
ax.set_xticklabels([])
plt.savefig('../figs/fig3_plot.svg', bbox_inches='tight')
# %%
