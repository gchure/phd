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

# Load the master data file
datadir = '../../data/'
df = pd.read_csv('../../data/ch2_induction/RazoMejia_2018.csv', comment='#')

# Now we remove the autofluorescence and delta values
df = df[(df.rbs != 'auto') & (df.rbs != 'delta')]


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


# Compute the mode and HPD for every parameter.
# max_idx = np.argmax(flat_lnprob, axis=0)
# param_fit = global_df.ix[max_idx, :]
# param_fit = param_fit.to_frame(name='mode')
# param_hpd = pd.DataFrame(columns=['hpd_min', 'hpd_max'])

# # Loop through and generate the dataframe.
# for column in global_df:
#     param_hpd = param_hpd.append(pd.Series(mwc.hpd(global_df[column], 0.95),
#                                            index=['hpd_min', 'hpd_max'], name=column))

# param_fit = pd.concat([param_fit, param_hpd], axis=1)

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
    leak = d.iloc[0].groupby('repressors')['fold_change_A'].agg(('mean', 'sem')).reset_index()
    sat = d.iloc[-1].groupby('repressors')['fold_change_A'].agg(('mean', 'sem')).reset_index()
    dyn_rng = d.iloc[-1]['fold_change_A'].values() - d.iloc[0]['fold_change_A'].values
    dyn_rng_mean, dyn_rng_sem = np.mean(dyn_rng), np.std(dyn_rng) / np.sqrt(len(dyn_rng))

    ax[0].errorbar(g[1], leak['mean'], leak['sem'], markersize=4, color=op_colors[g[0]],
                   markeredgewidth=0.75, markeredgecolor='white', label='__nolegend__')
    ax[1].errorbar(g[1], sat['mean'], sat['sem'], markersize=4, color=op_colors[g[0]],
                   markeredgewidth=0.75, markeredgecolor='white', label='__nolegend__')
    ax[2].errorbar(g[1], dyn_rng_mean, dyn_rng_sem, markersize=4, color=op_colors[g[0]],
                   markeredgewidth=0.75, markeredgecolor='white', label='__nolegend__')


#%%
    leak = leakiness(num_rep, ep_r, 4.5)

    sat = saturation(num_rep, ep_r, 4.5, ka_gf / ki_gf)
    dyn_rng = dyn_range(num_rep, ep_r, ka_gf / ki_gf, ep_ai=4.5)
    ec50 = EC50(ka_gf, ki_gf, 4.5, num_rep, ep_r)
    hill = effective_Hill(ka_gf, ki_gf, 4.5, num_rep, ep_r)
    ax[0].plot(num_rep, leak, '-', color=en_colors[op], lw=1)
    ax[1].plot(num_rep, sat, '-', color=en_colors[op], lw=1)
    ax[2].plot(num_rep, dyn_rng, '-', color=en_colors[op], lw=1)
    ax[3].plot(num_rep, ec50 / 1E6, '-', color=en_colors[op], lw=1)
    ax[4].plot(num_rep, hill, '-', color=en_colors[op], lw=1)

    # Compute and plot the credible regions.
    leak_cred = leakiness_cred_region(num_rep, global_df[op], 4.5)
    sat_cred = saturation_cred_region(
        num_rep, global_df[op], 4.5, global_df['Ka'], global_df['Ki'])
    dyn_cred = dyn_cred_region(
        num_rep, global_df['Ka'], global_df['Ki'], global_df[op])
    ec50_cred = ec50_cred_region(num_rep, global_df[op], 4.5, global_df['Ka'],
                                 global_df['Ki'])
    hill_cred = effective_hill_cred(num_rep, global_df[op], 4.5,
                                    global_df['Ka'], global_df['Ki'])

    ax[0].fill_between(num_rep, leak_cred[0, :], leak_cred[1, :],
                       color=en_colors[op], alpha=0.4)

    ax[1].fill_between(num_rep, sat_cred[0, :], sat_cred[1, :],
                       color=en_colors[op], alpha=0.4)
    ax[2].fill_between(num_rep, dyn_cred[0, :], dyn_cred[1, :],
                       color=en_colors[op], alpha=0.4)
    ax[3].fill_between(num_rep, ec50_cred[0, :] / 1E6, ec50_cred[1, :] / 1E6,
                       color=en_colors[op], alpha=0.4)
    ax[4].fill_between(num_rep, hill_cred[0, :], hill_cred[1, :],
                       color=en_colors[op], alpha=0.4)

    # Compute the EC50 and effective hill and plot.
    for j, R in enumerate(reps):
        # Load the single fit data.
        chain = glob.glob(
            '../../data/mcmc/SI_G_{0}_{1}.pkl'.format(op, rbs[j]))
        with open(chain[0], 'rb') as file:
            unpickler = pickle.Unpickler(file)
            flatchain2 = unpickler.load()
            flat_lnprob2 = unpickler.load()
        max_idx = np.argmax(flat_lnprob2)
        ea, ei, _, rep, ep = flatchain2[max_idx]
        rep *= 2
        rep_cred = 2 * mwc.hpd(flatchain2[:, 3], 0.95)
        print(ep)
        ka_fc, ki_fc = np.exp(-flatchain2[:, 0]), np.exp(-flatchain2[:, 1])
        ep_fc = flatchain2[:, 4]
        ka, ki = np.exp(-ea), np.exp(-ei)

        ec502 = EC50(ka, ki, 4.5, rep, ep)
        rep = np.array([rep, ])
        ec50_cred = ec50_cred_region(rep, ep_fc, 4.5, ka_fc,
                                     ki_fc)
        hill = effective_Hill(ka, ki, 4.5, rep, ep)
        hill_cred = effective_hill_cred(rep, ep_fc, 4.5, ka_fc,
                                        ki_fc)

        ax[3].plot(rep, ec502 / 1E6, 's', markerfacecolor='w',
                   markeredgecolor=en_colors[op], markersize=4, markeredgewidth=1)

        ax[4].plot(rep, hill, 's', markerfacecolor='w',
                   markeredgecolor=en_colors[op], markersize=4, markeredgewidth=1)
        ax[3].vlines(rep, ec50_cred[0] / 1E6,
                     ec50_cred[1] / 1E6, color=en_colors[op])
        ax[4].vlines(rep, hill_cred[0], hill_cred[1], color=en_colors[op])

        # Load data for leakiness and saturation plots
        dat = df[(df['repressors'] == (R / 2)) & (df['operator'] == op)]
        iptg = dat.IPTG_uM.unique()
        grouped = dat.groupby('IPTG_uM').fold_change_A.mean()

        unique_IPTG = dat['IPTG_uM'].unique()

        #   Slice the min and max IPTG values.
        leak_vals = dat[dat['IPTG_uM'] == np.min(unique_IPTG)].fold_change_A
        sat_vals = dat[dat['IPTG_uM'] == np.max(unique_IPTG)].fold_change_A

        # Compute the mean and standard errors of reach.
        mean_leak = np.mean(leak_vals)
        sem_leak = np.std(leak_vals) / np.sqrt(len(leak_vals))
        mean_sat = np.mean(sat_vals)
        sem_sat = np.std(sat_vals) / np.sqrt(len(sat_vals))

        # Plot the data for every point.
        rep = rep[0]
        ax[0].plot(rep, mean_leak, 'o', color=en_colors[op], ms=4)
        ax[0].errorbar(rep, mean_leak, yerr=sem_leak, color=en_colors[op])
        ax[1].plot(rep, mean_sat, 'o', color=en_colors[op], ms=4)
        ax[1].errorbar(rep, mean_sat, yerr=sem_sat, color=en_colors[op])
        d = drng[(drng.repressors == (0.5 * R)) & (drng.operator == op)]
        ax[2].plot(rep, d.dynamic_range, 'o', color=en_colors[op], ms=4)
        ax[2].errorbar(rep, d.dynamic_range, yerr=d.err, color=en_colors[op])


ylabs = ['leakiness', 'saturation', 'dynamic range', '$[EC_{50}]$ (M)',
         'effective Hill coefficient', '']
for i, a in enumerate(ax):
    a.set_xscale('log')
    a.set_xlabel('repressors per cell', fontsize=12)
    a.set_ylabel(ylabs[i], fontsize=12)
    a.set_xlim([1, 1E4])

plt.figtext(0., .96, '(A)', fontsize=12)
plt.figtext(0.33, .96, '(B)', fontsize=12)
plt.figtext(0.66, .96, '(C)', fontsize=12)
plt.figtext(0.0, .5, '(D)', fontsize=12)
plt.figtext(0.33, .5, '(E)', fontsize=12)
ax[0].set_yscale('log')
ax[3].set_yscale('log')
mwc.scale_plot(fig, 'two_row')
leg = ax[0].legend(title='    binding\n energy ($k_BT$)', loc='lower left',
                   fontsize=5)
leg.get_title().set_fontsize(6)
plt.tight_layout()
plt.savefig('../../figures/SI_figs/figS20.svg', bbox_inches='tight')
