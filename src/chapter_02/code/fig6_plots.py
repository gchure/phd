#%%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec 
import pandas as pd
import seaborn as sns
import phd.viz
import phd.stats
import pickle 
colors, palette = phd.viz.phd_style()
constants = phd.thermo.load_constants()

# Load the data set
data = pd.read_csv('../../data/ch2_induction/RazoMejia_2018.csv', comment='#')

# Load the flatchains for the prediction measurements. 
with open('../../data/ch2_induction/mcmc/SI_I_O2_R260.pkl', 'rb') as file:
    unpickler = pickle.Unpickler(file)
    gauss_flatchain = unpickler.load()
    gauss_flatlnprobability = unpickler.load()
ka_fc = np.exp(-gauss_flatchain[:, 0])[::100]
ki_fc = np.exp(-gauss_flatchain[:, 1])[::100]


#%%
# Compute the theoretical property curves. 
rep_range = np.logspace(0, 4, 200)
prop_df = pd.DataFrame([])
ops = {'O1':constants['O1'], 'O2':constants['O2'], 'O3':constants['O3']}
for op, ep in ops.items():
    for i, r in enumerate(rep_range):
        arch = phd.thermo.SimpleRepression(r, ep, ka=ka_fc, ki=ki_fc, 
                                          ep_ai=constants['ep_AI'], effector_conc=0)
        props = arch.compute_properties()
        for prop, val in props.items():
            if prop == 'leakiness':
                hpd_min, hpd_max = val, val
            else:
                hpd_min, hpd_max = phd.stats.compute_hpd(val, 0.95)
            prop_df = prop_df.append({'operator':op, 'binding_energy':ep,
                                      'property':prop, 'repressors':r,
                                      'hpd_min':hpd_min, 'hpd_max':hpd_max},
                                      ignore_index=True)

# %%
#Instantiate the figure canvas. 
fig, ax = plt.subplots(2, 3, figsize=(6, 3.5), dpi=100)
phd.viz.despine(ax.ravel())
# Format axes as needed. 
for i, a in enumerate(ax.ravel()):
    a.set_xscale('log')
    a.set_xlabel('repressors per cell')

ax[0, 0].set_yscale('log')
ax[0, 0].set_ylabel('leakiness')
ax[0, 1].set_ylabel('saturation')
ax[0, 2].set_ylabel('dynamic range')
ax[1, 0].set_ylabel(r'$EC_{50}$ [µM]')
ax[1, 0].set_yscale('log')
ax[1, 1].set_ylabel('effective Hill coefficient')
ax[1, 2].axis('off')

# Define the color palette
op_palette = sns.color_palette('viridis', n_colors=4)
op_color = {'O1':op_palette[0], 'O2': op_palette[1], 'O3':op_palette[2]}

# Assign axes to properties
axes = {'leakiness':ax[0, 0], 'saturation':ax[0, 1], 'dynamic_range':ax[0, 2],
        'EC50':ax[1, 0], 'effective_hill':ax[1, 1]}

# iterate through each property and operator
for g, d in prop_df.groupby(['operator', 'property']):
    if g[1] == 'leakiness':
        axes[g[1]].plot(rep_range, d['hpd_min'], color=op_color[g[0]], lw=1,
                        label=op)       
    else:
        axes[g[1]].fill_between(rep_range, d['hpd_min'].values, d['hpd_max'].values,
                            color=op_color[g[0]], lw=0.5, ec=colors['black'],
                            alpha=0.5)
leg = ax[0,0].legend(loc='lower left', fontsize=6, title='operator')
leg.get_title().set_fontsize(6)

# plot the data on top.of the measured properties 
for g, d in data[(data['repressors']>0) & 
    (data['operator']!='Oid')].groupby(['operator', 'repressors']):

    # Leakiness
    leak_d = d[d['IPTG_uM'] == 0]
    ax[0,0].errorbar(2 * g[1], leak_d['fold_change_A'].mean(), leak_d['fold_change_A'].sem(),
        fmt='o', linestyle='none', color=op_color[g[0]], ms=4,
        alpha=1, markeredgewidth=0.5, capsize=1, markeredgecolor='white')

    # Saturation 
    sat_d = d[d['IPTG_uM'] == d['IPTG_uM'].max()]
    ax[0,1].errorbar(2 * g[1], sat_d['fold_change_A'].mean(), sat_d['fold_change_A'].sem(),
        fmt='o', linestyle='none', color=op_color[g[0]], ms=4,
        alpha=1, markeredgewidth=0.5, capsize=1, markeredgecolor='white')

    # Dynamic range. 
    min_len = np.min(np.array([len(leak_d), len(sat_d)]))
    
    dyn_mean = np.mean(sat_d['fold_change_A'].values[:min_len] -\
                    leak_d['fold_change_A'].values[:min_len])
    dyn_sem = np.std(sat_d['fold_change_A'].values[:min_len] -\
                    leak_d['fold_change_A'].values[:min_len]) / min_len
    ax[0,2].errorbar(2 * g[1], dyn_mean, dyn_sem, fmt='o', linestyle='none', 
                    color=op_color[g[0]], alpha=1, markeredgewidth=0.5, capsize=1,
                    ms=4, markeredgecolor='white', linewidth=0.5)


    # Compute the inferred properties
    # Load the flat chains for the particular fit. 
    with open(f'../../data/ch2_induction/mcmc/SI_I_{g[0]}_R{2 * g[1]}.pkl', 'rb') as f:
        unpickler = pickle.Unpickler(f)
        flatchain = unpickler.load()
        lnprob = unpickler.load()
        ka = np.exp(-flatchain[:, 0])
        ki = np.exp(-flatchain[:, 1])
        max_idx = np.argmax(lnprob, axis=0)
        ea, ei, _= flatchain[max_idx]
        ka_mode = np.exp(-ea)
        ki_mode = np.exp(-ei)
    
    # Compute the appropriate properties. 
    arch_cred = phd.thermo.SimpleRepression(g[1], ep_r=constants[g[0]], ka=ka, ki=ki,
                        ep_ai=4.5, effector_conc=0).compute_properties()
    arch_mode = phd.thermo.SimpleRepression(g[1], ep_r=constants[g[0]], ka=ka_mode, 
                        ki=ki_mode, ep_ai=4.5, effector_conc=0).compute_properties()

    for i, p in enumerate(['EC50', 'effective_hill']):
        hpd_min, hpd_max = phd.stats.compute_hpd(arch_cred[p], 0.95)
        mode = arch_mode[p]
        ax[1, i].vlines(g[1], hpd_min, hpd_max, color=op_color[g[0]], lw=1)
        ax[1, i].plot(g[1], mode, marker='s', markerfacecolor=op_color[g[0]], 
                      markeredgecolor='w', lw=0.5, ms=4, markeredgewidth=0.5)

plt.tight_layout()
plt.savefig('../figs/fig6_properties.svg', bbox_inches='tight')

# %%
