#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import phd.viz
import phd.stats
import phd.bayes
import pystan
colors, plaette = phd.viz.phd_style()

#%%
# Load the fluctuation data and restrict only to the glucose 37 measurements. 
data = pd.read_csv('../../data/ch8_growth_si/analyzed_fluctuations.csv', comment='#')
data = data[(data['carbon']=='glucose') & (data['temp']==37)]

#%% Load the stan model .
model = pystan.StanModel('./calibration_factor.stan')

# %%
# Instantiate storage dataframes for the sampling results. 
samples = []

# Iterate through each replicate, perform the inference, add identifiers and
# store 
for g, d in data.groupby(['date', 'run_no']):
    # Assemble the data dictionary.
    data_dict = {'N':len(d), 'I1':d['I_1'], 'I2':d['I_2']}
    sample = model.sampling(data_dict, iter=5000)
    sample = sample.to_dataframe()

    # Add identifiers
    sample['date'] = g[0]
    sample['run_no'] = g[1]
    samples.append(sample)

# Concanate the sample data frames. 
samples = pd.concat(samples)


# %%
# Instantiate the figure
color_cycle = ['dark_purple', 'dark_orange', 'dark_blue', 'dark_green']
fig, ax = plt.subplots(1, 2, figsize=(6, 2.5), dpi=100, sharex=True)
phd.viz.despine(ax)
ax[0].set_xlabel('calibration factor [a.u. / LacI]')
ax[1].set_xlabel('calibration factor [a.u. / LacI]')
# ax[0].set_xscale('log')
ax[0].set_ylabel('empirical cumulative distribution')
ax[1].set_ylabel('posterior probability')
ax[0].set_xlim([100, 1200])
phd.viz.titlebox(ax[0], 'glucose, 37°C', color=colors['black'], size=6, 
                bgcolor='white', boxsize="10%")
phd.viz.titlebox(ax[1], 'glucose, 37°C', color=colors['black'], 
                bgcolor='white', size=6, boxsize="10%")

_samples = samples[(samples['date'] == 20181002) | (samples['date'] == 20190102) |
                   (samples['date'] == 20181011)]

# Iterate through each identifier in the sampling and compute the ecdf. 
iter = 0

for g, d in _samples.groupby(['date', 'run_no']):
    x, y = np.sort(d['alpha']), np.arange(0, len(d), 1) / len(d)
    x[0] = 100
    x[-1] = 3000
    ax[0].step(x, y, color=colors[color_cycle[iter]], lw=1)
    iter += 1

# Define the number of bins for a histogram. 
bins = np.linspace(200, 1200, 100)

# Iterate through each identifier and plot the histogram. 
iter = 0
for g, d in _samples.groupby(['date', 'run_no']):
    hist, _ = np.histogram(d['alpha'], bins=bins, density=True)
    ax[1].step(bins[:-1], hist , lw=1, color=colors[color_cycle[iter]], alpha=0.5)
    ax[1].fill_between(bins[:-1], np.zeros(len(hist)), hist, 
                            lw=1, color=colors[color_cycle[iter]], step='pre',
                            alpha=0.25)
    iter += 1
plt.tight_layout()
fig.text(0.01, 0.9, '(A)', fontsize=8) 
fig.text(0.5, 0.9, '(B)', fontsize=8) 
plt.savefig('../figs/ch8_figS5.pdf', bbox_inches='tight')
plt.savefig('../figs/ch8_figS5.png', bbox_inches='tight')

# %%
