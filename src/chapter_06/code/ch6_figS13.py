#%%
import glob
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import phd.viz
import phd.thermo
import seaborn as sns
colors, palette = phd.viz.phd_style()
constants = phd.thermo.load_constants()

# Load the master data file
datadir = '../../data/'
df = pd.read_csv('../../data/ch2_induction/RazoMejia_2018.csv', comment='#')

# Now we remove the autofluorescence and delta values
df = df[(df.rbs != 'auto') & (df.rbs != 'delta')]
df['repressors'] *= 2


# Load the flat-chain  used for parameter estimation
with open('../../data/ch6_induction_si/main_text_KaKi.pkl', 'rb') as file:
    unpickler = pickle.Unpickler(file)
    gauss_flatchain = unpickler.load()
    gauss_flatlnprobability = unpickler.load()

# map value of the parameters
max_idx = np.argmax(gauss_flatlnprobability, axis=0)
ea, ei, sigma = gauss_flatchain[max_idx]
ka, ki = np.exp(-ea), np.exp(-ei)
# Convert the flatchains to units of concentration.
ka_fc = np.exp(-gauss_flatchain[:, 0])
ki_fc = np.exp(-gauss_flatchain[:, 1])

# Separate the data for calculation of other properties.
operators = df.operator.unique()
energies = {'O1': -15.3, 'O2': -13.9, 'O3': -9.7}

# Compute the dynamic range
drs = []
for g, d in df.groupby(['operator']):
    unique_IPTG = d.IPTG_uM.unique()
    min_IPTG = np.min(unique_IPTG)
    max_IPTG = np.max(unique_IPTG)
    # Group the new data by repressors.
    grouped_rep = d.groupby(['rbs', 'date', 'username'])
    rbs_ind = {'HG104': 0, 'RBS1147': 1, 'RBS446': 2, 'RBS1027': 3,
               'RBS1': 4, 'RBS1L': 5}
    rep_dr = [[], [], [], [], [], []]
    rep_std = []

    # Deal with one notational mistake (?)
    for g_rep, d_rep in grouped_rep:
        if g_rep[2] != 'sloosbarnes':
            dr = d_rep[d_rep.IPTG_uM == max_IPTG].fold_change_A.values - \
                d_rep[d_rep.IPTG_uM == min_IPTG].fold_change_A.values
            rep_dr[rbs_ind[g_rep[0]]].append(dr[0])

    # Compute the means.
    for i, dr in enumerate(rep_dr):
        rep_dr[i] = np.mean(dr)
        rep_std.append(np.std(dr) / np.sqrt(len(dr)))

    reps = np.sort(df.repressors.unique())
    dr_df = pd.DataFrame([reps, rep_dr, rep_std]).T
    dr_df.columns = ['repressors', 'dynamic_range', 'err']
    dr_df.insert(0, 'operator', g)
    drs.append(dr_df)
drng = pd.concat(drs, axis=0)


# Load the global fit chains,
with open('../../data/ch6_induction_si/SI_E_global.pkl', 'rb') as file:
    unpickler = pickle.Unpickler(file)
    flatchain = unpickler.load()
    flat_lnprob = unpickler.load()


# Set the indexing for the MCMC dataframe.
index = ['log_ka', 'log_ki', 'sigma', 'HG104', 'RBS1147', 'RBS446',
         'RBS1027', 'RBS1', 'RBS1L', 'Oid', 'O1', 'O2', 'O3']

global_df = pd.DataFrame(flatchain, columns=index)
global_df['Ka'] = np.exp(-global_df['log_ka'])
global_df['Ki'] = np.exp(-global_df['log_ki'])
index = global_df.columns


#%%

# Define parameter ranges. 
rep_range = np.logspace(0, 4, 300)
c = np.logspace(-2, 4, 500)

# Make a new figure and plot the properties.
fig, ax = plt.subplots(2, 3, figsize=(6, 4))
ax = ax.ravel()
phd.viz.despine(ax)

# Format the axes
for a in ax:
    a.set_xlabel('repressors per cell')
    a.set_xscale('log')
for a in [ax[0], ax[3]]:
    a.set_yscale('log')
ax[0].set_ylabel('leakiness')
ax[1].set_ylabel('saturation')
ax[2].set_ylabel('dynamic range')
ax[3].set_ylabel('EC$_{50}$ [µM]')
ax[4].set_ylabel('effective Hill coefficient')
ax[-1].set_axis_off()


ops = ['O1', 'O2', 'O3']
op_colors = {op:c for op, c in zip(ops, sns.color_palette('viridis', n_colors=len(ops) + 1))}
props = ['leakiness', 'saturation', 'dynamic_range', 'EC50', 'effective_hill']
prop_ax = {p:a for p, a in zip(props, ax)}
for i, op in enumerate(ops):   

    # Compute the credible regions
    creds = {p: np.zeros((2, len(rep_range))) for p in props}

    # Iterate through each sampling and compute the properties. 
    for j, r in enumerate(rep_range):
        properties = phd.thermo.SimpleRepression(R=r, ep_r=global_df[op].values[::10],
                                                 ka=global_df['Ka'].values[::10], ki=global_df['Ki'].values[::10],
                                                 ep_ai=4.5, effector_conc=0).compute_properties()

        for p, v in properties.items():
           low, high = phd.stats.compute_hpd(v, 0.95)
           creds[p][0, j] = low
           creds[p][1, j] = high


    for p, v in creds.items():
        if p == 'leakiness':
            prop_ax[p].plot(rep_range, v[0, :], lw=0.75, label=op, 
                            color=op_colors[op]) 
        else:
            prop_ax[p].fill_between(rep_range, v[0, :], v[1, :], alpha=0.5,
                                    color=op_colors[op]) 

# Plot the data measurements. 
for g, d in df.groupby(['operator', 'repressors']):
    d.sort_values(['IPTG_uM'], inplace=True)
    leak = d[d['IPTG_uM']==0].groupby('repressors')['fold_change_A'].agg(('mean', 'sem')).reset_index()
    sat = d[d['IPTG_uM']==d['IPTG_uM'].max()].groupby('repressors')['fold_change_A'].agg(('mean', 'sem')).reset_index()
    _leak = d[d['IPTG_uM']==0]
    _sat = d[d['IPTG_uM']==d['IPTG_uM'].max()]
    min_val = np.min([len(_leak), len(_sat)])
    dyn_rng = _sat['fold_change_A'].values[:min_val] -\
              _leak['fold_change_A'].values[:min_val]
    dyn_rng_mean, dyn_rng_sem = np.mean(dyn_rng), np.std(dyn_rng) / np.sqrt(len(dyn_rng))

    ax[0].errorbar(g[1], leak['mean'], leak['sem'], markersize=4, color=op_colors[g[0]],
                   markeredgewidth=0.75, markeredgecolor='white', label='__nolegend__',
                   fmt='o')
    ax[1].errorbar(g[1], sat['mean'], sat['sem'], markersize=4, color=op_colors[g[0]],
                   markeredgewidth=0.75, markeredgecolor='white', label='__nolegend__',
                   fmt='o')
    ax[2].errorbar(g[1], dyn_rng_mean, dyn_rng_sem, markersize=4, color=op_colors[g[0]],
                   markeredgewidth=0.75, markeredgecolor='white', label='__nolegend__',
                   fmt='o')

    # To compute the EC50 and effective hill, use the inferred values for each
    # specific strain. 
    with open(f'../../data/ch2_induction/mcmc/SI_I_{g[0]}_R{g[1]}.pkl', 
              'rb') as f:
        unpickler = pickle.Unpickler(f)
        flatchain = unpickler.load()
        flat_lnprob = unpickler.load()

    # map value of the parameters
    max_idx = np.argmax(flat_lnprob, axis=0)
    ea, ei, _ = flatchain[max_idx]
    strain_ka_mode, strain_ki_mode = np.exp(-ea), np.exp(-ei)

    # Convert the flatchains to units of concentration.
    strain_ka_fc = np.exp(-flatchain[:, 0])
    strain_ki_fc = np.exp(-flatchain[:, 1])

    # Given the data, compute the effective hill and ec50
    _arch_mode = phd.thermo.SimpleRepression(R=g[1], ep_r=constants[g[0]],
                                             ka=strain_ka_mode, ki=strain_ki_mode,
                                             ep_ai=4.5, effector_conc=0).compute_properties()
    _arch_fc = phd.thermo.SimpleRepression(R=g[1], ep_r=constants[g[0]],
                                          ka=strain_ka_fc, ki=strain_ki_fc,
                                          ep_ai=4.5, effector_conc=0).compute_properties()
    arch_creds = {'EC50':np.zeros((2)), 'effective_hill':np.zeros(2)}
    for p in ['EC50', 'effective_hill']:
        arch_creds[p][:] = phd.stats.compute_hpd(_arch_fc[p], 0.95)

    # Plot the modes and cred regions. 
    ax[3].vlines(g[1], arch_creds['EC50'][0], arch_creds['EC50'][1], lw=0.75,
                color=op_colors[g[0]])
    ax[3].plot(g[1], _arch_mode['EC50'], 's', markersize=4, markeredgewidth=0.5,
              markeredgecolor='white', color=op_colors[g[0]])
    ax[4].vlines(g[1], arch_creds['effective_hill'][0], 
                arch_creds['effective_hill'][1], lw=0.75,
                color=op_colors[g[0]])
    ax[4].plot(g[1], _arch_mode['effective_hill'], 's', markersize=4, 
                markeredgewidth=0.5, markeredgecolor='white', color=op_colors[g[0]])

ax[0].legend(fontsize=6)
plt.tight_layout()
plt.savefig('../figs/figS13_plots.svg')

# %%
