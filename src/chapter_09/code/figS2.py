# -*- coding: utf-8 -*-
# %%
import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as plt
import statsmodels.tools.numdiff  as smnd
import fcsparser
import scipy.optimize
import phd.flow
import phd.stats
import phd.viz
from datetime import datetime
colors, palette = phd.viz.phd_style()

#%%
# Start difference
T_START = 138# in sec 

# Define the data directory
files = np.sort(glob.glob('../../data/ch9_mscl_si/flow/RP*_.*.fcs'))
means = []
time = []
for f in files:
    m, data = fcsparser.parse(f)
    gated = phd.flow.gaussian_gate(data, 0.4)
    mean_FITC = gated['FITC-A'].mean()
    time.append(m['$ETIM']) 
    means.append(mean_FITC)


#%% Compute the proper time vector.
pattern = "%H:%M:%S"
dt = [T_START/60] 
time_init = datetime.strptime(time[0], pattern)
for i in range(1, len(time)):
    t = datetime.strptime(time[i], pattern)
    diff = t - time_init
    dt.append((T_START + diff.seconds )/ 60)

# %% Compute the growth rate
growth_data = pd.read_csv('../../data/ch9_mscl_si/20180720_sfGFP_lb500_growth.csv')

# Trim the data to the linear region
lin_data = growth_data[(growth_data['time_min'] > 20) & (growth_data['time_min'] < 120)]
lin_data['log_A_A0'] = np.log(lin_data['od600'].values / lin_data.iloc[0]['od600'])
lin_data['adj_time'] = lin_data['time_min'] - lin_data.iloc[0]['time_min']
# Define a function to compute the grwoth rate. 
def log_post(lam, log_A_A0, time, neg=True):
    if neg == True:
        prefactor = -1
    else:
        prefactor = 1
    k = len(time)
    theo = time * lam
    return prefactor * -0.5 * np.log(k) * np.log(np.sum((log_A_A0 - theo)**2))

popt = scipy.optimize.minimize(log_post, 0.001, args=(lin_data['log_A_A0'], lin_data['adj_time']), 
                                method='powell')
hess = smnd.approx_hess(np.array([float(popt.x)]), log_post, args=(lin_data['log_A_A0'], lin_data['adj_time'], False))
cov = -np.linalg.inv(hess)
std = np.sqrt(cov)
mean_lam = float(popt.x)
err = np.sqrt(cov)[0][0]


# %%

# Instantiate the figure
fig, ax = plt.subplots(1, 2, figsize=(6, 3))
phd.viz.despine(ax)

# Format axes 
ax[0].set_ylim([0.3, 1.15])
ax[1].set_xlim([-5, 80])

# Add labels
ax[0].set_xlabel('time [min]')
ax[0].set_ylabel('normalized fluorescence')
ax[1].set_xlabel('time [min]')
ax[1].set_ylabel('$\log (A / A_0)$')
ax[0].text(-0.2, 1.0, '(A)', transform=ax[0].transAxes)
ax[1].text(-0.24, 1.0, '(B)', transform=ax[1].transAxes)
ax[0].set_xlabel('time [min]')
ax[0].set_ylabel('normalized fluorescence')

# Plot the maturation data
_ = ax[0].plot(dt, np.array(means)/means[-1], '-o', lw=1, color=colors['purple'], 
                       markeredgecolor='white', ms=4.5, markeredgewidth=0.5, 
                       label='measured\nintensity')
_ = ax[0].fill_betweenx(np.linspace(0, 1.15, 500), 30, 40, color='white', zorder=-1, 
                label='observed\n growth rate')

time_range = np.linspace(-5, 80, 500)
_ = ax[1].plot(time_range, mean_lam * time_range, '-', color=colors['purple'], lw=0.75, 
                label='best fit')
_ = ax[1].fill_between(time_range, (mean_lam - std)[0][0] * time_range,
                    (mean_lam + std)[0][0] * time_range, color=colors['purple'],
                    alpha=0.5, label='__nolegend__')
_ = ax[1].plot(lin_data['adj_time'], lin_data['log_A_A0'], 'o', color=colors['purple'], ms=4,
               label='data', markeredgecolor='white', markeredgewidth=0.5)
leg = ax[1].legend(loc='upper left', fontsize=6, title='$t_{double} = 35 \pm 1$ min') 
leg.get_title().set_fontsize(6)
plt.tight_layout()
plt.savefig('../figs/figS2.pdf', bbox_inches='tight')
plt.savefig('../figs/figS2.png', bbox_inches='tight')

# %%
