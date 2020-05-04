
#%%
import numpy as np
import pandas as pd
import scipy.stats
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import phd.thermo
import phd.viz
colors, palette = phd.viz.phd_style()


# Load the prior predictive check data. 
data = pd.read_csv('../../data/ch7_mutants_si/epRA_prior_predictive_checks.csv')

# Define the percentiles to plot. 
percs = [99, 95, 80, 50, 20, 10, 5]

# Compute the percentiles of the simulations. 
grouped = data.groupby(['IPTGuM'])
df = pd.DataFrame([], columns=['percentile', 'IPTGuM', 'fc_low', 'fc_high'])
for g, d in grouped:
    for p in percs:     
        remainder = 100 - p
        low = remainder / 2
        upper = p + remainder / 2 
        _percs = np.percentile(d['fc_draw'], [low, upper])
        df = df.append({'percentile': p,
                         'IPTGuM': g,
                         'fc_low':_percs[0],
                         'fc_high': _percs[1]},
                      ignore_index=True)

# Compute the PDFs
epRA_range = np.linspace(-38, 10, 500)
sigma_range = np.linspace(0, 1.2, 500)
epRA_pdf = scipy.stats.norm.pdf(epRA_range, loc=-10, scale=5)
sigma_pdf = scipy.stats.norm.pdf(sigma_range, loc=0, scale=0.1)

# Determine the number of simulations performed and generate range values for drawn samples
n_sims = len(data['ep_RA'].unique())
y = np.random.normal(0, 0.002, size=n_sims)
       
# Set up the figure canvas. 
fig = plt.figure(figsize=(6, 3))
gs = gridspec.GridSpec(7, 7)
ax0 = fig.add_subplot(gs[:3, :3])
ax1 = fig.add_subplot(gs[4:, :3])
ax2 = fig.add_subplot(gs[:, 4:])

# Apply special formatting
ax = [ax0, ax1, ax2]
phd.viz.despine(ax)

ax0.set_xlim([-40, 10])
ax0.set_ylim([0, 0.1])
ax1.set_xlim([0, 0.5])
ax1.set_ylim([0, 5])
ax2.set_xscale('log')
ax2.set_xlim([1E-2,5E3])

# Add appropriate axis labels. 
for i, a in enumerate(ax[:-1]):
    a.set_yticks([])
    a.set_ylabel('probability')
ax0.set_xlabel(r'DNA binding energy ($\Delta\varepsilon_{RA}$) [$k_BT$]')
ax1.set_xlabel('standard deviation ($\sigma$)')
ax2.set_xlabel('IPTG [µM]')
ax2.set_ylabel('fold-change')


# Add panel labels. 
ax0.text(-0.2, 1.05, '(A)', fontsize=6.5, transform=ax0.transAxes)
ax1.text(-0.2, 1.05, '(B)', fontsize=6.5, transform=ax1.transAxes)
ax2.text(-0.3, 1.02, '(C)', fontsize=6.5, transform=ax2.transAxes)

# Plot the PDFs.
_ = ax0.plot(epRA_range, epRA_pdf, color=colors['purple'], lw=0.75)
_ = ax0.fill_between(epRA_range, epRA_pdf, color=colors['light_purple'],
     alpha=0.25)
_ = ax1.plot(sigma_range, sigma_pdf, color=colors['purple'], lw=0.75)
_ = ax1.fill_between(sigma_range, sigma_pdf, color=colors['light_purple'],
    alpha=0.25)

# Plot the PDFs
_ = ax0.hlines(0.095,-40, 10, color='w', lw=16)
_y = np.ones(n_sims) * 0.095 + np.random.normal(0, 0.002, size=n_sims)
_ = ax0.plot(data['ep_RA'].unique(), _y, ',', color='k')
_ = ax1.hlines(4.55, 0, 0.5, color='w', lw=14)
_y = np.ones(n_sims) * 4.5 + np.random.normal(0, 0.1, size=n_sims)
_ = ax1.plot(data['sigma'].unique(), _y, ',', color='k')


# Plot the percentiles of the fold-change. 
perc_colors = sns.color_palette('Purples', n_colors=len(percs))
for c, p in zip(perc_colors, percs):
    _d = df[df['percentile']==p]
    _ = ax2.fill_between(_d['IPTGuM'].unique(), _d['fc_low'], _d['fc_high'],
                        color=c, label=p, alpha=0.5)    
_leg = ax2.legend(title='percentile', handlelength=1, fontsize=6, loc='lower right')
_leg.get_title().set_fontsize(6)
plt.savefig('../figs/ch7_figS2.pdf', bbox_inches='tight')
plt.savefig('../figs/ch7_figS2.png', bbox_inches='tight')

# %%
