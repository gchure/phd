
# %%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pystan
import pandas as pd
import phd.viz
import phd.stats
colors, palette = phd.viz.phd_style()

# Load the necessary data sets.
data = pd.read_csv('../../data/ch5_mscl/mscl_survival_data.csv')
model = pystan.StanModel('./generic_logistic_regression.stan')

# %%
# Use area as a predictor variable.
data.loc[data['shock_class'] == 'slow', 'idx'] = 1
data.loc[data['shock_class'] == 'fast', 'idx'] = 2
data_dict = {'J': int(2), 'N': len(data), 'trial': data['idx'].values.astype(
    int), 'output': data['survival'].astype(int), 'predictor': data['area']}
samples = model.sampling(data=data_dict, iter=5000, chains=4)
area_samples = samples.to_dataframe()

# %%
# Use log shock rate as a predictor
data_dict = {'J': 1, 'N': len(data), 'trial': np.ones(len(data)).astype(
    int), 'output': data['survival'].astype(int), 'predictor': data['flow_rate'].values}
samples = model.sampling(data=data_dict, iter=5000, chains=4)
rate_samples = samples.to_dataframe()
# %%
# Use both channel number and log shock rate as a predictor variable
bivariate_model = pystan.StanModel(
    './bivariate_logistic_regression.stan')
data_dict = {'N': len(data), 'output': data['survival'].astype(int), 'predictor_1': np.log(
    data['effective_channels'].values), 'predictor_2': data['flow_rate'].values}
samples = bivariate_model.sampling(data=data_dict, iter=5000, chains=4)
bivariate_samples = samples.to_dataframe()

# %%
# Generate the plots.
fig, ax = plt.subplots(2, 2, figsize=(6, 5))
phd.viz.despine(ax.ravel())
# Add figure panel labels.
fig.text(0, 0.95, '(A)', fontsize=8)
fig.text(0.5, 0.95, '(B)', fontsize=8)
fig.text(0, 0.45, '(C)', fontsize=8)
fig.text(0.5, 0.45, '(D)', fontsize=8)

ax = ax.ravel()
# Add labels
for i in range(4):
    ax[i].tick_params(labelsize=8)
    if i < 3:
        ax[i].set_ylabel('survival probability')
        ax[i].set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax[0].set_xlabel('cell area [µm$^2$]')
ax[1].set_xlabel('cell area [µm$^2$]')
ax[2].set_xlabel('shock rate [Hz]')
ax[3].set_xlabel('effective channel number')
ax[3].set_ylabel('shock rate [Hz]')
area_range = np.linspace(0, 40, 500)
titles = ['slow shock (< 1.0 Hz)', 'fast shock ($\geq$ 1.0 Hz)']
color_key = [colors['purple'], colors['blue']]
fill_colors = [colors['light_purple'], colors['light_blue']]
for i in range(2):
    ax[i].set_ylim([-0.2, 1.2])
    ax[i].set_xlim([0, 40])
    phd.viz.titlebox(ax[i], titles[i], size=8, pad=0.05, boxsize='12%', color=colors['black'],
                    bgcolor='white')
    _ = ax[i].set_xlabel('cell area [µm$^2$]', fontsize=8)
    # Set the title

    # Plot the credible regions.
    cred_region = np.zeros((2, len(area_range)))
    for j, a in enumerate(area_range):
        beta_0 = area_samples['beta_0[{}]'.format(i+1)]
        beta_1 = area_samples['beta_1[{}]'.format(i+1)]
        logit = beta_0 + beta_1 * a
        prob = (1 + np.exp(-logit))**-1
        cred_region[:, j] = phd.stats.compute_hpd(prob, 0.95)
    _ = ax[i].fill_between(area_range, cred_region[0, :], cred_region[1, :], color=fill_colors[i], alpha=0.5)
    _ = ax[i].plot(area_range, cred_region[0, :], color=color_key[i], lw=0.75, label='__nolegend__')
    _ = ax[i].plot(area_range, cred_region[1, :], color=color_key[i], lw=0.75, label='__nolegend__')

    # Set the top and bottom stripes for plotting the measurements 
    _ = ax[i].hlines(1.15, 0, 40, color='w', lw=10, label='__nolegend__')
    _ = ax[i].hlines(-0.15, 0, 40, color='w', lw=10, label='__nolegend__')

    #  Plot the actual data.
    _data = data[data['idx']==i + 1]
    _surv = _data[_data['survival'] == True]['area']
    _death  = _data[_data['survival'] == False]['area']
    ys = np.random.normal(loc=1.15, scale=0.01, size=len(_surv))
    yd = np.random.normal(loc=-0.15, scale=0.01, size=len(_death))
    _ = ax[i].plot(_surv, ys, 'k.', ms=1, alpha=0.5)
    _ = ax[i].plot(_death, yd, 'k.', ms=1, alpha=0.5)


# Plot the shock rate prediction.
ax[2].set_xlim([0, 2.5])
ax[2].set_ylim([-0.2, 1.2])
rate_range = np.linspace(0, 2.5, 500)

# Compute the credible region.
cred_region = np.zeros((2, len(rate_range)))
for i, r in enumerate(rate_range):
    logit = rate_samples['beta_0[1]'] + rate_samples['beta_1[1]'] * r
    prob = (1 + np.exp(-logit))**-1
    cred_region[:, i] = phd.stats.compute_hpd(prob, 0.95)
_ = ax[2].fill_between(rate_range, cred_region[0, :], cred_region[1, :], color='slategray', alpha=0.5)
_ = ax[2].plot(rate_range, cred_region[0, :], color='k', lw=0.75, label='__nolegend__')
_ = ax[2].plot(rate_range, cred_region[1, :], color='k', lw=0.75, label='__nolegend__')

# Plot the cell measurements.
_ = ax[2].hlines(1.15, 0, 2.5, lw=10, color='w')
_ = ax[2].hlines(-0.15, 0, 2.5, lw=10, color='w')
surv = data[data['survival']==True]['flow_rate']
death = data[data['survival']==False]['flow_rate']
ys = np.random.normal(loc=1.15, scale=0.01, size=len(surv))
yd = np.random.normal(loc=-0.15, scale=0.01, size=len(death))
_ = ax[2].plot(surv, ys, 'k.', ms=2, alpha=0.5, label='__nolegend__')
_ = ax[2].plot(death, yd, 'k.', ms=2, alpha=0.5, label='__nolegend__')

# Plot the bivariate regression.
chan_range = np.logspace(0, 3, 500)
X, Y = np.meshgrid(chan_range, rate_range)
beta_0 = bivariate_samples['beta_0'].median()
beta_1 = bivariate_samples['beta_1'].median()
beta_2 = bivariate_samples['beta_2'].median()
logit = beta_0 + beta_1 * np.log(X) + beta_2 * Y
prob = (1 + np.exp(-logit))**-1

# Plot the contours of probability. 
_ = ax[3].contourf(X, Y, prob, cmap='viridis', levels=[0, 0.2, 0.4, 0.6, 0.8, 0.9, 1])
_cont = ax[3].contour(X, Y, prob, colors='w', levels = [0, 0.2, 0.4, 0.6, 0.8, 0.9, 1])

# Plot the points of cells.
_ = ax[3].vlines(-10, -0.2, 2.7, color='w', lw=15)
_ = ax[3].vlines(1005, -0.2, 2.7, color='w', lw=15)
_ = ax[3].hlines(-0.15, 0, 1010, color='w', lw=15)
_ = ax[3].hlines(2.65, 0, 1010, color='w', lw=15)
_ = ax[3].set_ylim([-0.2, 2.7])
_ = ax[3].set_xlim([-15, 1020])

_surv = data[data['survival']==True]
_death = data[data['survival']==False]
ys = np.random.normal(2.65, 0.01, size=len(surv))
yd = np.random.normal(-0.15, 0.01, size=len(death))
xs = np.random.normal(1005, 5, size=len(_surv))
xd = np.random.normal(-5, 2, size=len(_death))
_ = ax[3].plot(_surv['effective_channels'], ys, 'k.', ms=2, alpha=0.75, label='__nolegend__')
_ = ax[3].plot(_death['effective_channels'], yd, 'k.', ms=2, alpha=0.75, label='__nolegend__')
_ = ax[3].plot(xs, _surv['flow_rate'], 'k.', ms=1, alpha=0.5, label='__nolegend__')
_ = ax[3].plot(xd, _death['flow_rate'], 'k.', ms=1, alpha=0.5, label='__nolegend__')
label_pos = [(100, 0.5), (200, 1.0), (400, 1.5), (600, 1.5), (850, 0.5)]
plt.clabel(_cont, inline=1, fontsize=8,  manual=label_pos, fmt='%0.2f')

plt.tight_layout()
plt.savefig('../figs/ch9_figS10.pdf', bbox_inches='tight')
plt.savefig('../figs/ch9_figS10.png', bbox_inches='tight')
# %%
