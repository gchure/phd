#%%
import numpy as np
import pandas as pd
import scipy.stats
import matplotlib.pyplot as plt
import seaborn as sns
import phd.viz
import phd.flow
import fcsparser
colors, palette = phd.viz.phd_style()

# Load example flow cytometry data
_, data = fcsparser.parse('../../data/ch2_induction/example_flow/20160813_r1_wt_O2_RBS1027_0uMIPTG.fcs')

blues = sns.color_palette('Blues_r', n_colors=7)

# Set the range  of alpha.
alpha_range = [0.8, 0.6, 0.4, 0.25, 0.05]

# Generate an understandable legend.
fig, ax = plt.subplots(1,1, figsize=(4, 3))
phd.viz.despine(ax)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('forward scatter [a.u.]', fontsize=8)
ax.set_ylabel('side scatter [a.u.]', fontsize=8)
ax.plot(data['FSC-A'], data['SSC-A'], 'k.', rasterized=True, alpha=0.5,
ms=0.1, label='__nolegend__')
ax.plot([], [], 'k.', label=1)
for i, a in enumerate(alpha_range):
    gated = phd.flow.gaussian_gate(data, alpha=a)
    ax.plot([], [], '.', label=a, color=blues[i])
    ax.plot(gated['FSC-A'], gated['SSC-A'], '.', color=blues[i], rasterized=True, alpha=0.5,
    ms=0.1, label='__nolegend__')

leg = ax.legend(title=r'gating fraction $\alpha$', fontsize=8)
leg.get_title().set_fontsize(8)
plt.savefig('../figs/ch6_figS8.pdf', bbox_inches='tight')
plt.savefig('../figs/ch6_figS8.png', bbox_inches='tight')


# %%
