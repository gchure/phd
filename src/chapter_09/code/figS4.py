# %%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import pandas as pd
import scipy.stats
import phd.viz
import phd.stats
colors, palette = phd.viz.phd_style()
# Load in the MCMC data.
data = pd.read_csv('../../data/ch9_mscl_si/complete_mcmc_traces.csv')

fig = plt.figure(figsize=(6, 3))
gs = gridspec.GridSpec(10, 10)
ax0 = fig.add_subplot(gs[:, 0:5])
ax1 = fig.add_subplot(gs[:4, 6:])
ax2 = fig.add_subplot(gs[6:, 6:])
ax0.set_ylim([2500, 4100])
ax0.set_xlim([4.75, 6])
fig.text(0.01, 0.9, '(A)', fontsize=8)
fig.text(0.55, 0.9, '(B)', fontsize=8)
fig.text(0.55, 0.4, '(C)', fontsize=8)

# Format the axes
phd.viz.despine([ax0, ax1, ax2])

# Plot the KDE if the samples.
_ = sns.kdeplot(data['hyper_A_mu'], data['hyper_alpha_mu'], ax=ax0, shade=True,
cmap=plt.cm.Purples)
ax0.set_ylabel('calibration factor [a.u. / MscL channel]', fontsize=8)
ax0.set_xlabel('average cell area [µm$^2$]', fontsize=8)

# Add a fake legend
ax1.plot([], [], color=colors['purple'], lw=1, label='replicate\n parameter')
ax1.plot([], [], color=colors['orange'], lw=1, label='hyper-\nparameter')
ax1.legend(fontsize=6, handlelength=1, loc='upper right')

# Evaluate the KDE for the low-level parameters
alpha_range = np.linspace(2000, 8000, 500)
area_range = np.linspace(4, 6.5, 500)
for i in range(6):
    # Evaluate the KDE and normalize.
    alpha_kernel = scipy.stats.gaussian_kde(data['alpha__{}'.format(i)])
    area_kernel = scipy.stats.gaussian_kde(data['avg_A__{}'.format(i)])
    alpha_fit = alpha_kernel(alpha_range)
    area_fit = area_kernel(area_range)
    alpha_fit *= np.sum(alpha_fit)**-1
    area_fit *= np.sum(area_fit)**-1

    # Plot the distributions.
    _ = ax1.plot(alpha_range, alpha_fit, color=colors['purple'])
    _ = ax2.plot(area_range, area_fit, color=colors['purple'])
    _ = ax1.fill_between(alpha_range, alpha_fit, alpha=0.2, color=colors['light_purple'])
    _ = ax2.fill_between(area_range, area_fit, alpha=0.2, color=colors['light_purple'])


# Plot the hyper parameters
hyper_alpha_kernel = scipy.stats.gaussian_kde(data['hyper_alpha_mu'])
hyper_area_kernel = scipy.stats.gaussian_kde(data['hyper_A_mu'])
hyper_alpha_fit = hyper_alpha_kernel(alpha_range)
hyper_area_fit = hyper_area_kernel(area_range)
hyper_alpha_fit *= np.sum(hyper_alpha_fit)**-1
hyper_area_fit *= np.sum(hyper_area_fit)**-1
_ = ax1.plot(alpha_range, hyper_alpha_fit, color=colors['orange'], lw=0.75) 
_ = ax1.fill_between(alpha_range, hyper_alpha_fit, color=colors['orange'], alpha=0.4, zorder=100) 
_ = ax2.plot(area_range, hyper_area_fit, color=colors['orange'], lw=0.75) 
_ = ax2.fill_between(area_range, hyper_area_fit, color=colors['orange'], alpha=0.4, zorder=100) 
ax1.set_xlabel('calibration factor [a.u./channel]')
ax2.set_xlabel('average cell area [µm$^2$]')
ax1.set_ylabel('$\propto$ probability')
ax2.set_ylabel('$\propto$ probability')
for a in (ax1, ax2):
    a.set_yticks([])
plt.tight_layout()
plt.savefig('../figs/ch9_figS4.png', bbox_inches='tight')
plt.savefig('../figs/ch9_fig4.pdf', bbox_inches='tight')


# %%
