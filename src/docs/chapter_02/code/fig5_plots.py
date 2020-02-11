#%%
import numpy as np
import matplotlib.pyplot as plt
import  matplotlib.gridspec as gridspec 
import pandas as pd
import phd.viz
import phd.stats
import pickle 
colors, palette = phd.viz.phd_style()

data = pd.read_csv('../../data/ch2_induction/RazoMejia_2018.csv', comment='#')
params = pd.read_csv('../../data/ch2_induction/RazoMejia_KaKi_estimates.csv')

# Now we remove the autofluorescence and delta values
data = data[(data['rbs'] != 'auto') & 
            (data['rbs'] != 'delta') & 
            (data['operator']!='Oid')]

# Load the flat-chain
with open('../../data/ch2_induction/mcmc/SI_I_O2_R260.pkl', 'rb') as file:
    unpickler = pickle.Unpickler(file)
    gauss_flatchain = unpickler.load()
    gauss_flatlnprobability = unpickler.load()

# Map value of the parameters
max_idx = np.argmax(gauss_flatlnprobability, axis=0)
ea, ei, sigma = gauss_flatchain[max_idx]
lnprob =  gauss_flatchain[:, 2][::100]
ka_fc = np.exp(-gauss_flatchain[:, 0])[::100]
ki_fc = np.exp(-gauss_flatchain[:, 1])[::100]

#%%
# Define the IPTG concentrations to evaluate
c_range = np.logspace(-2, 4, 200)
c_range[0] = 0 
rep_range = np.logspace(0, 4, 200)

# Set the colors for the strains
rep_colors = {22: colors['light_red'], 
              60: colors['light_brown'],  
              124: colors['light_green'], 
              260: colors['light_orange'], 
              1220: colors['light_purple'], 
              1740: colors['light_blue']} 
edge_colors = {22: colors['dark_red'], 
              60: colors['dark_brown'],  
              124: colors['dark_green'], 
              260: colors['dark_orange'], 
              1220: colors['dark_purple'], 
              1740: colors['dark_blue']} 

# Define the operators and their respective energies
operators = ['O1', 'O2', 'O3']
energies = {'O1': -15.3, 'O2': -13.9, 'O3': -9.7}


#%%
sampling_df = pd.DataFrame(np.array([ka_fc, ki_fc, lnprob]).T, 
                          columns=['ka', 'ki', 'logprob'])
sampling_df.sort_values('logprob', inplace=True)

# ##############################################################################
# PREDICTED INDUCTION PROFILES
# ##############################################################################
# Set up the dataframe for the fold-change. 
fc_df = pd.DataFrame()

# Iterate through each operator, repressor, and IPTG concentration to calculate
# the fold-change
for op, op_en in energies.items():
    for r, _ in rep_colors.items():
        for i, c in enumerate(c_range):
            arch = phd.thermo.SimpleRepression(r, op_en, ka=ka_fc, ki=ki_fc, 
                                               ep_ai=4.5, effector_conc=c)
            
            fc_min, fc_max = phd.stats.compute_hpd(arch.fold_change(), 0.95)
            fc_df = fc_df.append({'fc_min': fc_min,
                                  'fc_max': fc_max,
                                  'repressors': r,
                                  'IPTGuM': c,
                                  'operator': op,
                                  'binding_energy': op_en}, ignore_index=True)

#%%

# Set up the canvas for the operator plots. 
fig = plt.figure(figsize=(5, 4))
gs = gridspec.GridSpec(8, 8)
ax0 = fig.add_subplot(gs[:4, :4])
ax1 = fig.add_subplot(gs[:4, 4:])
ax2 = fig.add_subplot(gs[4:, :4])
ax3 = fig.add_subplot(gs[4:6, 4:])
ax4 = fig.add_subplot(gs[6:, 4:])

ax3.set_xticklabels([])
ax3.set_ylim([0.5E-3, 10])
ax4.set_ylim([0.5E-2, 1E4])
ax4.set_xticklabels(['', '22', '60', '124', '260', '1220', '1740'])
ax4.set_xlabel('repressors per cell')
ax3.set_ylabel('$K_I$ [µM]')
ax4.set_ylabel('$K_A$ [µM]')
for a in [ax3, ax4]:
    a.set_yscale('log')
    a.set_xlim([0, 7])

ax = [ax0, ax1, ax2, ax3, ax4]
op_ax = {'O1':ax0, 'O2':ax1, 'O3':ax2}

for o, a in op_ax.items():
    a.set_xscale('symlog', linthreshx=1E-2)
    a.set_xlabel('IPTG [µM]')
    a.set_ylabel('fold-change')
    a.set_ylim([-0.05, 1.2])
    phd.viz.titlebox(a, f'operator {o}', color=colors['black'], 
                    bgcolor='white', boxsize="12%")

# Plot the data. 
for g, d in data.groupby(['operator', 'repressors', 'IPTG_uM']):
    if (g[0] == 'O2') & (g[1] == 130):
        face = 'white'
        edge = colors['dark_orange']
        ew = 0.5
    else:
        face = rep_colors[2 * g[1]]
        edge = edge_colors[2 * g[1]]
        ew = 0.25
    mean_fc = d['fold_change_A'].mean()
    sem_fc = d['fold_change_A'].std() / np.sqrt(len(d))
    _ = op_ax[g[0]].errorbar(g[-1], mean_fc, sem_fc, fmt='.', lw=0.5, 
        markerfacecolor=face, color=edge, markeredgecolor=edge, 
        markeredgewidth=ew, ms=6, label='__nolegend__')
plt.subplots_adjust(hspace=0.8, wspace=0.65)

axes={'O1':ax[0], 'O2':ax[1], 'O3':ax[2]}
for _, a in axes.items():
    a.set_xscale('symlog', linthreshx=1E-2)
    a.set_ylim([-0.01, 1.15])
    a.set_ylabel('fold-change')
    a.set_xlabel('IPTG [µM]')

for g, d in fc_df.groupby('operator'):
    for r, _d in d.groupby('repressors'):
        _ = axes[g].fill_between(_d['IPTGuM'], _d['fc_min'], _d['fc_max'],
                                 alpha=0.75, facecolor=rep_colors[r], label=int(r),
                                 edgecolor=edge_colors[r], lw=0.25)
    phd.viz.titlebox(axes[g], f'operator {g}', color=colors['black'], 
                    boxsize='12%')

# Isolate the data for the fit strain
fit_strain = data[(data['operator']=='O2') & 
                  (data['repressors']==130)].groupby(['IPTG_uM']
                  ).agg(('mean', 'sem')).reset_index()


# Plot the estimated parameters
op_glyphs = {'O1':'s', 'O2':'o', 'O3':'^'}
fudge = {'O1': -0.25, 'O2': 0, 'O3': 0.25}
rep_pos = {22:1, 60:2, 124:3, 260:4, 1220:5, 1740:6}
params.sort_values(by='repressors', inplace=True)
for g, d in params.groupby(['operator', 'repressors']):
    ka = d[d['parameter']=='ka']
    ki = d[d['parameter']=='ki']
    ax3.vlines(rep_pos[g[1]] + fudge[g[0]], ki['hpd_min'], ki['hpd_max'],
        color=edge_colors[g[1]], 
        label='__nolegend__')
    ax3.plot(rep_pos[g[1]] + fudge[g[0]], ki['mode'], op_glyphs[g[0]],
            markerfacecolor='w', color=edge_colors[g[1]], ms=4,
            label='__nolegend__')
    ax4.vlines(rep_pos[g[1]] + fudge[g[0]], ka['hpd_min'], ka['hpd_max'],
        color=edge_colors[g[1]], label='__nolegend__')
    ax4.plot(rep_pos[g[1]] + fudge[g[0]], ka['mode'], op_glyphs[g[0]],
            markerfacecolor='w', color=edge_colors[g[1]], ms=4,
             label='__nolegend__')

for op, m in op_glyphs.items():
    ax3.plot([], m, markerfacecolor='w', markeredgecolor='k', label=op,
    ms=3)
ax3.hlines(0.53, 0, 7, zorder=1, linestyle='--', color='k')
ax4.hlines(139, 0, 7, zorder=1, linestyle='--', color='k')
ax3.legend(loc='lower right', ncol=3, fontsize=6)
# Add the legend. 
leg = ax[0].legend(loc='upper left', title='rep. per cell', fontsize=6, handlelength=1.5)
leg.get_title().set_fontsize(6)
# plt.tight_layout()
plt.savefig('../figs/fig5_induction_profiles.svg', bbox_inches='tight')

# %%
