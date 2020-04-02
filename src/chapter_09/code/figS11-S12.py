# -*- coding: utf-8 -*-
# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import phd.stats
import phd.viz 
import phd.mscl
colors, palette = phd.viz.phd_style()

def compute_mean_sem(df):
    n_tot = len(df)
    n_surv = np.sum(df['survival'].astype(int))
    mean_chan = df['effective_channels'].mean()
    prob = n_surv / n_tot
    sem_chan = df['effective_channels'].std() / np.sqrt(len(df))
    prob_err = n_surv * (n_tot - n_surv) / n_tot**3
    out_dict = {'mean_chan': mean_chan, 'p_survival': prob, 'chan_err': sem_chan,
                'prob_err': np.sqrt(prob_err)}
    return pd.Series(out_dict)

# Load the data with the three shock logistic regression.
data = pd.read_csv('../../data/ch9_mscl_si/mscl_survival_data_three_shock.csv')

# Load the samples and the stats.
samples = pd.read_csv('../../data/ch9_mscl_si/three_shock_complete_mcmc_traces.csv')
pooled_samples = pd.read_csv('../../data/ch9_mscl_si/pooled_complete_mcmc_traces.csv')
indiv_samples = pd.read_csv('../../data/ch9_mscl_si/complete_mcmc_shock_rate_idx.csv')
pooled_data = pd.read_csv('../../data/ch9_mscl_si/pooled_mscl_survival_data.csv')

idx = {0: 'slow', 1: 'medium', 2: 'fast'}

# %%
# Instantiate the complicated figure
fig = plt.figure(figsize=(6,5))
gs = gridspec.GridSpec(22, 30)

# Add panel labels.
fig.text(0.05, 0.9, '(A)', fontsize=8)
fig.text(0.5, 0.9, '(B)', fontsize=8)
fig.text(0.05, 0.48, '(C)', fontsize=8)
fig.text(0.5, 0.48, '(D)', fontsize=8)
# Top panels
surv_ax1 = fig.add_subplot(gs[0, 0:14])
surv_ax2 = fig.add_subplot(gs[0, 16:])
death_ax1 = fig.add_subplot(gs[8, 0:14])
death_ax2 = fig.add_subplot(gs[8, 16:])
surv_ax3 = fig.add_subplot(gs[12, 0:14])
surv_ax4 = fig.add_subplot(gs[12, 16:])
death_ax3 = fig.add_subplot(gs[21, 0:14])
death_ax4 = fig.add_subplot(gs[21, 16:])
points = [surv_ax1, surv_ax2, surv_ax3, surv_ax4,
          death_ax1, death_ax2, death_ax3, death_ax4]

# Survival prob axes.
prob_ax1 = fig.add_subplot(gs[1:8, 0:14])
prob_ax2 = fig.add_subplot(gs[1:8, 16:])
prob_ax3 = fig.add_subplot(gs[13:21, 0:14])
prob_ax4 = fig.add_subplot(gs[13:21, 16:])
curves = [prob_ax1, prob_ax2, prob_ax3, prob_ax4]

# Format the various axes
for i, ax in enumerate(points):
    # ax.set_xticks([])
    if i <= 3:
        ax.set_xticks([])
        ax.spines['bottom'].set_visible(False)
    ax.set_yticks([])
    ax.set_facecolor('none')
    ax.set_xlim([0, 1E3])
    ax.tick_params(labelsize=8)
    # ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)

for i, ax in enumerate(curves):
    ax.tick_params(labelsize=8)
    ax.set_xlim([0, 1E3])
    ax.set_ylim([-0.02, 1.02])
    ax.set_xticklabels([])
    ax.spines['bottom'].set_visible(False)
    ax.set_xticks([])
    if i % 2 == 1:
        ax.set_ylabel('')
        ax.set_yticklabels([])
    else:
        ax.set_ylabel('survival probablity', fontsize=8)
    if i < 2:
        ax.set_xticklabels([])

bottoms = [death_ax1, death_ax2, death_ax3, death_ax4]
for i, ax in enumerate(bottoms):
    ax.set_xlabel('effective channel number', fontsize=8)

# Set the titles
tops = [surv_ax1, surv_ax2, surv_ax3, surv_ax4]
titles = ['slow shock ($< 0.5$ Hz)', 'intermediate shock (0.5 - 1.0 Hz)',
          'fast shock ($>$ 1.0 Hz)', 'all shocks']

for i, ax in enumerate(tops):
    ax.set_title(
        titles[i], fontsize=8, y=1.1,
        bbox={'boxstyle':'square', 'facecolor':'white'})
    ax.set_xticklabels([])

#  Plot the survival curves and credible regions.
chan_range = np.logspace(0, 3, 500)
for i in range(3):

    # Set the credible regions.
    cred_region = np.zeros((2, len(chan_range)))
    for j, c in enumerate(chan_range):
        logit = samples['beta_0__{}'.format(
            i)] + samples['beta_1__{}'.format(i)] * np.log(c)
        prob = (1 + np.exp(-logit))**-1
        cred_region[:, j] = phd.stats.compute_hpd(prob, 0.95)

    # Fill the shaded cred regions.
    _ = curves[i].fill_between(
        chan_range, cred_region[0, :], cred_region[1, :], color=colors['light_red'], alpha=0.5)

    # Draw the outlines.
    _ = curves[i].plot(chan_range, cred_region[0, :],
                       color=colors['red'], lw=0.5)
    _ = curves[i].plot(chan_range, cred_region[1, :],
                       color=colors['red'], lw=0.5)

# Plot the data.
rates = ['slow', 'medium', 'fast']
for i, r in enumerate(rates):
    _data = data[data['shock_class'] == r]
    surv = _data[_data['survival'] == True]
    death = _data[_data['survival'] == False]
    ys = np.random.normal(0, 0.1, len(surv))
    yd = np.random.normal(0, 0.1, len(death))
    _ = tops[i].plot(surv['effective_channels'], ys, 'k.', ms=1, alpha=0.75)
    _ = bottoms[i].plot(death['effective_channels'],
                        yd, 'k.', ms=1, alpha=0.75)

    # Plot the binned data.
    bin_width = 50
    binned = phd.mscl.density_binning(
        _data, bin_width=bin_width, groupby='shock_class', input_key='effective_channels', min_cells=20)
    grouped = binned.groupby('bin_number').apply(compute_mean_sem)
    _ = curves[i].errorbar(grouped['mean_chan'], grouped['p_survival'], xerr=grouped['chan_err'],
                           yerr=grouped['prob_err'], color=colors['black'], lw=0.75, linestyle='none',
                           marker='o', ms=4,  label='{} channels/bin'.format(bin_width), zorder=100,
                           markeredgewidth=0.5, markeredgecolor='white')

    # Plot the strain data.
    grouped = _data.groupby(['rbs']).apply(compute_mean_sem)
    _ = curves[i].errorbar(grouped['mean_chan'], grouped['p_survival'], 
                           xerr=grouped['chan_err'], yerr=grouped['prob_err'],
                           color=colors['blue'], lw=0.75, ms=4, marker='o', 
                           label='1 SD mutant / bin', linestyle='none', zorder=100,
                           markeredgecolor='white', markeredgewidth=0.5)

# Plot the pooled data.
surv = data[data['survival'] == True]
death = data[data['survival'] == False]
ys = np.random.normal(0, 0.1, len(surv))
yd = np.random.normal(0, 0.1, len(death))
surv_ax4.plot(surv['effective_channels'], ys, 'k.', ms=1, alpha=0.75)
death_ax4.plot(death['effective_channels'], yd, 'k.', ms=1, alpha=0.75)


# Plot the pooled data credible regions.
cred_region = np.zeros((2, len(chan_range)))
for i, c in enumerate(chan_range):
    logit = pooled_samples['beta_0'] + pooled_samples['beta_1'] * np.log(c)
    prob = (1 + np.exp(-logit))**-1
    cred_region[:, i] = phd.stats.compute_hpd(prob, 0.95)
_ = prob_ax4.fill_between(
    chan_range, cred_region[0, :], cred_region[1, :], color=colors['light_red'], alpha=0.5, label='__nolegend__')
_ = prob_ax4.plot(
    chan_range, cred_region[0, :], color=colors['red'], lw=0.75, label='__nolegend__')
_ = prob_ax4.plot(
    chan_range, cred_region[1, :], color=colors['red'], lw=0.75, label='__nolegend__')


# Plot the binned pooled data
pooled_data['idx'] = 1
binned = phd.mscl.density_binning(
    pooled_data, bin_width=bin_width, groupby='idx', input_key='effective_channels', min_cells=20)
grouped = binned.groupby('bin_number').apply(compute_mean_sem)
_ = prob_ax4.errorbar(grouped['mean_chan'], grouped['p_survival'], xerr=grouped['chan_err'],
                      yerr=grouped['prob_err'], linestyle='none', lw=0.75, ms=4, marker='o', color=colors['black'], 
                      label='{} channels / bin'.format(bin_width), zorder=100,
                      markeredgecolor='white', markeredgewidth=0.5)

# Plot the strain bins.
grouped = pooled_data.groupby(['rbs']).apply(compute_mean_sem)
_ = prob_ax4.errorbar(grouped['mean_chan'], grouped['p_survival'], xerr=grouped['chan_err'], yerr=grouped['prob_err'],
                           color=colors['blue'], lw=0.75, ms=4, marker='o', 
                           label='1 SD mutant / bin', linestyle='none', zorder=100,
                           markeredgecolor='white', markeredgewidth=0.5)

_leg = prob_ax1.legend(fontsize=8, handlelength=0.75)
plt.subplots_adjust(hspace=0.8)
plt.savefig('../figs/ch9_figS11.pdf', bbox_inches='tight')
plt.savefig('../figs/ch9_figS11.png', bbox_inches='tight')



# %%
flow_rates = np.sort(data['flow_rate'].unique())


# Set up the complicated figure.
fig, ax = plt.subplots(3, 4, figsize=(8, 5), sharex=False, sharey=True)
phd.viz.despine(ax.ravel())

# Properly format and label.
for i in range(4):
    ax[-1, i].set_xlabel('effective channel number', fontsize=8)
    if i < 3:
        ax[i, 0].set_ylabel('survival probability', fontsize=8)
ax = ax.ravel()
for i in range(len(ax) - 5):
    ax[i].set_xticklabels([])

ax[7].set_xlabel('effective channel number', fontsize=8)
for a in ax:
    a.tick_params(labelsize=8)
    a.set_xlim([0, 1000])
    a.set_ylim([-0.15, 1.15])
ax[-1].axis('off')

# Iterate through each shock group and plot the prediction, credible region, and bins.
for i, r in enumerate(flow_rates):
    # Plot the credible regions.
    cred_region = np.zeros((2, len(chan_range)))
    for j, c in enumerate(chan_range):
        logit = indiv_samples['beta_0__{}'.format(i)] + indiv_samples['beta_1__{}'.format(i)] * np.log(c)
        prob = (1 + np.exp(-logit))**-1
        cred_region[:, j] = phd.stats.compute_hpd(prob, 0.95)
    _ = ax[i].fill_between(chan_range, cred_region[0, :], cred_region[1, :], color=colors['light_red'], alpha=0.5)
    _ = ax[i].plot(chan_range, cred_region[0, :], color=colors['light_red'])
    _ = ax[i].plot(chan_range, cred_region[1, :], color=colors['light_red'])
    _ = ax[i].hlines(1.1, 0, 1E3, color='w', lw=11)
    _ = ax[i].hlines(-0.1, 0, 1E3, color='w', lw=11)
    # Add the title for the shock rate.
    # Isolate the data to be plotted.
    _data = data[data['flow_rate']==r]
    surv = _data[_data['survival']==True]
    death = _data[_data['survival']==False]
    # Plot the points.
    ys = np.random.normal(loc=1.1, scale=0.01, size=len(surv))
    yd = np.random.normal(loc=-0.1, scale=0.01, size=len(death))
    _ = ax[i].plot(surv['effective_channels'], ys, 'k.', ms=1, alpha=0.5, label='__nolegend__')
    _ = ax[i].plot(death['effective_channels'], yd, 'k.', ms=1, alpha=0.5, label='__nolegend__')
    rbs_binned = _data.groupby('rbs').apply(compute_mean_sem)
    rbs_df = pd.DataFrame(rbs_binned).reset_index()
    _ = ax[i].errorbar(rbs_df['mean_chan'], rbs_df['p_survival'], xerr=rbs_df['chan_err'], yerr=rbs_df['prob_err'],
        ms=4, color=colors['blue'], lw=1, linestyle='none', fmt='o', zorder=1000, label='1 SD mutant / bin',
        markeredgecolor='white', markeredgewidth=0.5)
    if i == 0:
        phd.viz.titlebox(ax[i], '{:.3f} Hz; N = {}'.format(r, len(_data)), size=8, color=colors['black'], bgcolor='white',
        pad=0.05, boxsize='12%')
    else:
        phd.viz.titlebox(ax[i], '{:.2f} Hz; N = {}'.format(r, len(_data)), size=8, color=colors['black'], bgcolor='white',
        pad=0.05, boxsize='12%')


_leg = ax[0].legend(fontsize=6)
plt.savefig('../figs/ch9_figS12.pdf', bbox_inches='tight')
plt.savefig('../figs/ch9_figS12.png', bbox_inches='tight')

# %%
