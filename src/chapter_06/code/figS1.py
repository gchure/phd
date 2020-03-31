#%%
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import phd.viz
colors, palette = phd.viz.phd_style()

# Load in the datasets
data_a = pd.read_csv('../../data/ch6_induction_si/figS1_partA.csv')
data_b = pd.read_csv('../../data/ch6_induction_si/figs1_partB.csv')
data_O2 = pd.read_csv('../../data/ch2_induction/RazoMejia_2018.csv', comment='#')
data_O2['repressors'] *= 2
data_O2 = data_O2[(data_O2['operator'] == 'O2') & (data_O2['repressors'] == 260)] 

# Define necessary parameters.
ka = data_a[data_a.parameter == 'logKA']
ki = data_a[data_a.parameter == 'logKI']
c_range = np.logspace(-8, -2, 500)
R = 260
ep_r = -13.9
_colors = sns.color_palette('magma', n_colors=len(data_b)+1)

def foldchange(c, R, ep_ai, ep_r, ep_a, ep_i):
    mwc_term = (1 + c / ep_a)**2 * (1 + np.exp(-ep_ai)) / ((1 + c / ep_a)**2 +
                                                           np.exp(-ep_ai) * (1 + c / ep_i)**2)
    fc = (1 + mwc_term * (R / 4.6E6) * np.exp(-ep_r))**-1
    return fc
# Define the figure axis and labels.
fig, ax = plt.subplots(1, 2, figsize=(6, 2.25))
phd.viz.despine(ax)
_ = ax[0].set_xlabel(r'allosteric parameter $\Delta\varepsilon_{AI}\,(k_BT)$')
_ = ax[0].set_ylabel(r'best-fit parameter value')
_ = ax[1].set_xscale('log')
_ = ax[1].set_xlabel('IPTG [µM]')
_ = ax[1].set_ylabel('fold-change')

# Plots for panel (A)
_ = ax[0].plot(ka.ep, ka.bestfit, '-', color=colors['orange'],
               label=r'$\mathrm{log}\, \frac{K_A}{1\mathrm{M}}$')
_ = ax[0].plot(ki.ep, ki.bestfit, '-', color=colors['purple'],
               label=r'$\mathrm{log}\, \frac{K_A}{1\mathrm{M}}$')

# Plot the curves
for i in range(len(data_b)):
    ep_ai = data_b.iloc[i]['ep_ai']
    ka = np.exp(data_b.iloc[i]['log_ka'])
    ki = np.exp(data_b.iloc[i]['log_ki'])
    fc = foldchange(c_range, R, ep_ai, ep_r, ka, ki)
    _ = ax[1].plot(c_range * 1E6, fc, label=ep_ai, color=_colors[i])

# plot the data.
grouped = data_O2.groupby(['IPTG_uM']).fold_change_A
for group, data in grouped:
    mean_fc = np.mean(data)
    mean_sem = np.std(data) / np.sqrt(len(data))
    _ = ax[1].errorbar(group, mean_fc, mean_sem,
                       linestyle='none', color=colors['orange'], 
                       fmt='o', markeredgecolor='white', markeredgewidth=0.5,
                       label='__nolegend__', ms=4.5)


# Add the legends and labels.
_ = ax[0].legend(loc='lower left')
leg = ax[1].legend(loc='upper left', title=r"""allosteric parameter
       $\Delta\varepsilon_{AI}$ $(k_BT)$""", bbox_to_anchor=(1, 1),
       fontsize=6)
leg.get_title().set_fontsize(6)
fig.text(0, 0.95, '(A)', fontsize=8)
fig.text(0.45, 0.95, '(B)', fontsize=8)

# Format and save the figure.
plt.tight_layout()
plt.savefig('../figs/ch6_figS1.pdf', bbox_inches='tight')
plt.savefig('../figs/ch6_figS1.png', bbox_inches='tight')


# %%
