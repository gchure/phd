#%%
import numpy as np
import pandas as pd
import scipy.stats as st
import matplotlib.pyplot as plt
import glob
import phd.viz
import phd.flow
colors, palette = phd.viz.phd_style()

# Load the fold-change measurements
flow_data = pd.read_csv('../../data/ch2_induction/RazoMejia_2018.csv', comment="#")
mic_data = pd.read_csv('../../data/ch6_induction_si/mic_comparison_data.csv')

# Load the intensity distribution measurements. 
flow_intensities = []
replicate = 1
for f in glob.glob('../../data/ch6_induction_si/5000*/*uM*.csv'):
    _df = pd.read_csv(f, comment='#')
    gate = phd.flow.gaussian_gate(_df, 0.4)
    gate['replicate'] = replicate
    replicate += 1
    flow_intensities.append(gate)
flow_intensities = pd.concat(flow_intensities, sort=False)

mic_intensities = pd.concat([pd.read_csv(f, comment='#') for f in glob.glob('../../data/ch6_induction_si/5000uM_IPTG_distributions/*microscopy*.csv')])
mic_intensities = mic_intensities[(mic_intensities['IPTG_uM']==5000) & 
                                  (mic_intensities['rbs']=='RBS1027') & 
                                  (mic_intensities['operator']=='O2')]

# Set up a new data frame for the fold-change comparison
_flow_data = flow_data[['IPTG_uM', 'operator', 'fold_change_A', 'repressors']]
_flow_data['method'] = 'flow'
_flow_data.rename(columns={'fold_change_A':'fold_change'}, inplace=True)
_mic_data = mic_data[['IPTG_uM', 'operator', 'repressors', 'fold_change']]
_mic_data['method'] = 'microscopy'
comparison_df = pd.concat([_flow_data, _mic_data])
comparison_df = comparison_df[(comparison_df['repressors'] == 130) &
                              (comparison_df['operator']=='O2')]
comparison_df = comparison_df.groupby(['method', 'IPTG_uM', 'operator', 
                       'repressors'])['fold_change'].agg(('mean', 'sem')).reset_index()

# %%
# Set up the figure canvas, 
fig, ax = plt.subplots(1, 3, figsize=(8, 3))
phd.viz.despine(ax)

ax[0].set_xlabel('fold-change (microscopy)')
ax[0].set_ylabel('fold-change\n(flow cytometry)')
ax[0].set_xlim([-0.35, 1.1])
ax[0].set_ylim([-0.35, 1.1])
ax[1].set_xlim([-0.5, 0.75])
ax[2].set_xticks([0, 1, 2])
ax[2].set_xlim([-0.5, 2.5])
ax[2].set_xticklabels(['variance', 'skewness', 'kurtosis'])
ax[2].set_yscale('log')
ax[1].set_xlabel('normalized intensity\nabout mean')
ax[1].set_ylabel('cumulative distribution')
ax[2].set_ylabel('central moment value')


# Plot the fold-change agreement. 
ax[0].plot([-1.5, 1.5], [-1.5, 1.5], 'k-', label='__nolegend__', lw=0.75)
mic = comparison_df[comparison_df['method']=='microscopy']
flow = comparison_df[comparison_df['method']=='flow']
ax[0].errorbar(mic['mean'], flow['mean'], xerr=mic['sem'], yerr=flow['sem'],
    fmt='o', markersize=5, color=colors['blue'], markeredgecolor='white',
    markeredgewidth=0.75, linestyle='none', lw=0.75, label='__nolegend__')

# Plot the normalized moments of the intensities (microscopy) 
label = False
for g, d in mic_intensities.groupby(['date']):
    x, y  = np.sort(d['mean_intensity'].values), np.arange(len(d)) / len(d)
    x_norm = (x - x[0]) / (x[-1] - x[0])
    x_norm -= x_norm.mean()
    if label == False:
        label = 'microscopy'
    else:
        label = '__nolegend__'
    # ax[1].step(x_norm, y, color=colors['purple'], lw=0.75, alpha=0.75, label=label)
    
    # compute the moments. 
    variance = st.moment(x_norm, 2)
    skewness = st.moment(x_norm, 3)
    kurtosis = st.moment(x_norm, 4)

    ax[2].plot(np.random.normal(0, 0.1), variance, 'o',
                color=colors['purple'], markeredgecolor='white', ms=5,
                markeredgewidth=0.75, label=label, alpha=0.75)
    ax[2].plot(np.random.normal(1, 0.1), skewness, 'o', ms=5,
                color=colors['purple'], markeredgecolor='white', 
                markeredgewidth=0.75, label='__nolegend__', alpha=0.75)
    ax[2].plot(np.random.normal(2, 0.1), kurtosis, 'o', ms=5,
                color=colors['purple'], markeredgecolor='white', 
                markeredgewidth=0.75, label='__nolegend__', alpha=0.75)


_x = np.sort(flow_intensities['FITC-A'].values)
_x_norm = (_x - _x[0]) / (_x[-1] - _x[0])
_x_norm -= _x_norm.mean()
bins = np.linspace(-0.5, 0.5, 45)
ax[1].hist(_x_norm, bins, density=True, alpha=0.5, label='flow cytometry',
           color=colors['red'])
_x = np.sort(mic_intensities['mean_intensity'].values)
_x_norm = (_x - _x[0]) / (_x[-1] - _x[0])
_x_norm -= _x_norm.mean()
ax[1].hist(_x_norm, bins, density=True, alpha=0.5, label='microscopy',
        color=colors['blue'])

# Plot the normalized moments of the intensities (flow) 
label = False
for g, d in flow_intensities.groupby(['replicate']):
    x, y  = np.sort(d['FITC-A'].values), np.arange(len(d)) / len(d)
    x_norm = (x - x[0]) / (x[-1] - x[0])
    x_norm -= x_norm.mean()
    if label == False:
        label = 'flow cytometry'
    else:
        label = '__nolegend__'
    # ax[1].step(x_norm, y, color=colors['orange'], lw=0.75, alpha=0.75, label=label)
    
    # compute the moments. 
    variance = st.moment(x_norm, 2)
    skewness = st.moment(x_norm, 3)
    kurtosis = st.moment(x_norm, 4)

    ax[2].plot(np.random.normal(0, 0.1), variance, 'o',
                color=colors['orange'], markeredgecolor='white', ms=5,
                markeredgewidth=0.75, label=label, alpha=0.75)
    ax[2].plot(np.random.normal(1, 0.1), skewness, 'o', ms=5,
                color=colors['orange'], markeredgecolor='white', 
                markeredgewidth=0.75, label='__nolegend__', alpha=0.75)
    ax[2].plot(np.random.normal(2, 0.1), kurtosis, 'o', ms=5,
                color=colors['orange'], markeredgecolor='white', 
                markeredgewidth=0.75, label='__nolegend__', alpha=0.75)

ax[1].legend(fontsize=6, handlelength=1)
ax[2].legend(fontsize=6)
plt.tight_layout()
plt.savefig('../figs/ch6_figS10_histogram.pdf', bbox_inches='tight')
# plt.savefig('../figs/ch6_figS10.png', bbox_inches='tight')

# %%j