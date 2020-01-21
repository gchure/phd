#%%
import os
import glob
import pickle
import re
import numpy as np
import pandas as pd
import sys
import phd.viz
import phd.stats
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import matplotlib.colors as plc
import seaborn as sns
import scipy.stats
colors, palette = phd.viz.phd_style()

df = pd.read_csv('../../data/ch2_induction/RazoMejia_2018.csv', comment='#')

# Now we remove the autofluorescence and delta values
df = df[(df['rbs'] != 'auto') & (df['rbs'] != 'delta')]

#===============================================================================
# O2 RBS1027
#===============================================================================
# Load the flat-chain
with open('../../data/ch2_induction/mcmc/SI_I_O2_R260.pkl', 'rb') as file:
    unpickler = pickle.Unpickler(file)
    gauss_flatchain = unpickler.load()
    gauss_flatlnprobability = unpickler.load()

# map value of the parameters
max_idx = np.argmax(gauss_flatlnprobability, axis=0)
ea, ei, sigma = gauss_flatchain[max_idx]

ka_fc = np.exp(-gauss_flatchain[:, 0])
ki_fc = np.exp(-gauss_flatchain[:, 1])

#%%
#===============================================================================
# Plot the theory vs data for all 4 operators with the credible region
#===============================================================================

# Define the IPTG concentrations to evaluate
c_range = np.logspace(-2, 4, 200)
c_range[0] = 0


# Set the colors for the strains
rep_colors = {22: 'red', 
              60: 'brown',  
              124: 'green', 
              260: 'orange', 
              1220: 'purple', 
              1740:'blue'} 

# Define the operators and their respective energies
operators = ['O1', 'O2', 'O3']
energies = {'O1': -15.3, 'O2': -13.9, 'O3': -9.7, 'Oid': -17}

# Initialize the figure.
fig, ax = plt.subplots(3, 3, figsize=(6, 6))
ax = ax.ravel()


# Plot the predictions.
for i, op in enumerate(operators):
    data = df[df.operator == op]
    # loop through RBS mutants
    for j, rbs in enumerate(df.rbs.unique()):
        # plot the theory using the parameters from the fit.
        if (op == 'O2') & (rbs == 'RBS1027'):
            label = None
        else:
            label = df[df.rbs == rbs].repressors.unique()[0] * 2

        R = df[df['rbs']==rbs]['repressors'].unique()[0] * 2

        # plot 95% HPD region using the variability in the MWC parameters
        fc_cred_region = np.zeros((2, len(c_range)))
        for k, c in enumerate(c_range):
            arch = phd.thermo.SimpleRepression(R=R, ep_r=energies[op],
                                               ka=ka_fc[::100], ki=ki_fc[::100], 
                                               ep_ai=4.5, effector_conc=c
                                               ).fold_change()
            cred_region[:, k] = phd.stats.compute_hpd(arch, 0.95)
j
        ax[i].fill_between(c_range, cred_region[0, :], cred_region[1, :],
                           alpha=0.6, color=colors[f'{rep_colors[R]}'])  

        fc_mean = data[data.rbs == rbs].groupby('IPTG_uM')['fold_change_A'].agg(('mean', 'sem')).reset_index()

    
        # Distinguish between the fit data and the predictions
        if (op == 'O2') & (rbs == 'RBS1027'):
            face = 'w'
            ecolor = colors[f'dark_{rep_colors[R]}']
            ax[i].errorbar(fc_mean['IPTG_uM'], fc_mean['mean'], yerr=fc_mean['sem'], 
                       linestyle='none', fmt='o', markerfacecolor=face, 
                       markeredgecolor=ecolor, ms=4, markeredgewidth=0.75,
                       lw=1, color=ecolor, capsize=1)


    # Add operator and binding energy labels.
    phd.viz.titlebox(ax[i], 
        f'{op} | ' +  r'$\Delta\varepsilon_{RA} =$'  + f'{energies[op]} ' + '$k_BT$', 
         bgcolor=colors['grey'], color=colors['black'], fontsize=8, pad=0.05, boxsize="16%")
    ax[i].set_xscale('symlog', linthreshx=1E-2, linscalex=0.5)
    ax[i].set_xlabel('IPTG [µM]', fontsize=8)
    ax[i].set_ylabel('fold-change', fontsize=8)
    ax[i].set_ylim([-0.01, 1.1])
    ax[i].set_xticks([0, 1E-1, 1E1, 1E3])
    ax[i].tick_params(labelsize=8)


# Compute and plot the properties. 
rep_range = np.logspace(0, 3, 200)
leak_cred_region = np.zeros((2, len(rep_range)))
sat_cred_region = np.zeros((2, len(rep_range)))
dyn_rng_cred_region = np.zeros((2, len(rep_range)))
ec50_cred_region = np.zeros((2, len(rep_range)))
hill_cred_region = np.zeros((2, len(rep_range)))

# for i, r in enumerate(rep_range):
#     # Define the architecture. 
#     arch = phd.thermo.SimpleRepression(R=r, ep_r) 


plt.tight_layout()

#%%

# def leakiness(num_rep, ep_r, ep_ai, n_ns=4.6E6):
    # pact = 1 / (1 + np.exp(-ep_ai))
    # return (1 + pact * (num_rep / n_ns) * np.exp(-ep_r))**-1
# 
# 
# def saturation(num_rep, ep_r, ep_ai, ka_ki, n_sites=2, n_ns=4.6E6):
    # pact = 1 / (1 + np.exp(-ep_ai) * ka_ki**n_sites)
    # return (1 + pact * (num_rep / n_ns) * np.exp(-ep_r))**-1
# 
# 
# def saturation_cred_region(num_rep, ep_r, ep_ai, ka_flatchain, ki_flatchain,
                        #    n_sites=2, n_ns=4.6E6, mass_frac=0.95):
    # pact = 1 / (1 + np.exp(-ep_ai) * (ka_flatchain / ki_flatchain)**n_sites)
    # cred_region = np.zeros([2, len(num_rep)])
    # for i, R in enumerate(num_rep):
        # fc = (1 + pact * (R / n_ns) * np.exp(-ep_r))**-1
        # cred_region[:, i] = mwc.hpd(fc, mass_frac)
    # return cred_region
# 
# 
# def dyn_range(num_rep, ep_r, ka_ki, ep_ai=4.5, n_sites=2, n_ns=4.6E6):
    # pact_leak = 1 / (1 + np.exp(-ep_ai))
    # pact_sat = 1 / (1 + np.exp(-ep_ai) * (ka_ki)**n_sites)
    # leak = (1 + pact_leak * (num_rep / n_ns) * np.exp(-ep_r))**-1
    # sat = (1 + pact_sat * (num_rep / n_ns) * np.exp(-ep_r))**-1
    # return sat - leak
# 
# 

# 
# def pact(IPTG, K_A, K_I, e_AI):
    # '''
    # Computes the probability that a repressor is active
    # Parameters
    # ----------
    # IPTG : array-like
        # Array of IPTG concentrations in uM
    # K_A : float
        # Dissociation constant for active repressor
    # K_I : float
        # Dissociation constant for inactive repressor
    # e_AI : float
        # Energetic difference between the active and inactive state
    # Returns
    # -------
    # probability that repressor is active
    # '''
    # pact = (1 + IPTG * 1 / K_A)**2 / \
        # (((1 + IPTG * 1 / K_A))**2 + np.exp(-e_AI) * (1 + IPTG * 1 / K_I)**2)
    # return pact
# 
# 
# def fold_change(IPTG, K_A, K_I, e_AI, R, Op):
    # '''
    # Computes fold-change for simple repression
    # Parameters
    # ----------
    # IPTG : array-like
        # Array of IPTG concentrations in uM
    # K_A : float
        # Dissociation constant for active repressor
    # K_I : float
        # Dissociation constant for inactive repressor
    # e_AI : float
        # Energetic difference between the active and inactive state
    # R : float
        # Number of repressors per cell
    # Op : float
        # Operator binding energy
    # Returns
    # -------
    # probability that repressor is active
    # '''
    # return 1 / (1 + R / 5E6 * pact(IPTG, K_A, K_I, e_AI) * np.exp(-Op))
# 
# 
# def dyn_cred_region(num_rep, ka_flatchain, ki_flatchain, ep_r, mass_frac=0.95, epsilon=4.5):
    # cred_region = np.zeros([2, len(num_rep)])
    # ka_ki = ka_flatchain / ki_flatchain
    # for i, R in enumerate(num_rep):
        # drng = dyn_range(R, ep_r, ka_ki, ep_ai=epsilon)
        # cred_region[:, i] = mwc.hpd(drng, mass_frac)
    # return cred_region
# 
# 
# def EC50(K_A, K_I, e_AI, R, Op):
    # '''
    # Computes the concentration at which half of the repressors are in the active state
    # Parameters
    # ----------
    # K_A : float
        # Dissociation constant for active repressor
    # K_I : float
        # Dissociation constant for inactive repressor
    # e_AI : float
        # Energetic difference between the active and inactive state
    # Returns
    # -------
    # Concentration at which half of repressors are active (EC50)
    # '''
    # t = 1 + (R / 4.6E6) * np.exp(-Op) + (K_A / K_I)**2 * \
        # (2 * np.exp(-e_AI) + 1 + (R / 4.6E6) * np.exp(-Op))
    # b = 2 * (1 + (R / 4.6E6) * np.exp(-Op)) + \
        # np.exp(-e_AI) + (K_A / K_I)**2 * np.exp(-e_AI)
    # return K_A * ((K_A / K_I - 1) / (K_A / K_I - (t / b)**(1 / 2)) - 1)
# 
# 
# def ec50_cred_region(num_rep, Op, e_AI, K_A, K_I,
                    #  mass_frac=0.95):
    # cred_region = np.zeros([2, len(num_rep)])
    # for i, R in enumerate(num_rep):
        # t = 1 + (R / 4.6E6) * np.exp(-Op) + (K_A / K_I)**2 * \
            # (2 * np.exp(-e_AI) + 1 + (R / 4.6E6) * np.exp(-Op))
        # b = 2 * (1 + (R / 4.6E6) * np.exp(-Op)) + \
            # np.exp(-e_AI) + (K_A / K_I)**2 * np.exp(-e_AI)
        # ec50_rng = K_A * ((K_A / K_I - 1) / (K_A / K_I - (t / b)**(1 / 2)) - 1)
        # cred_region[:, i] = mwc.hpd(ec50_rng, mass_frac)
    # return cred_region
# 
# 
# def effective_Hill(K_A, K_I, e_AI, R, Op):
    # '''
    # Computes the effective Hill coefficient
    # Parameters
    # ----------
    # K_A : float
        # Dissociation constant for active repressor
    # K_I : float
        # Dissociation constant for inactive repressor
    # e_AI : float
        # Energetic difference between the active and inactive state
    # Returns
    # -------
    # effective Hill coefficient
    # '''
    # c = EC50(K_A, K_I, e_AI, R, Op)
    # return 2 / (fold_change(c, K_A, K_I, e_AI, R, Op) - fold_change(0, K_A, K_I, e_AI, R, Op)) *\
        # (-(fold_change(c, K_A, K_I, e_AI, R, Op))**2 * R / 5E6 * np.exp(-Op) *
        #  2 * c * np.exp(-e_AI) * (1 / K_A * (1 + c / K_A) * (1 + c / K_I)**2 - 1 / K_I *
                                #   (1 + c / K_A)**2 * (1 + c / K_I)) / ((1 + c / K_A)**2 + np.exp(-e_AI) * (1 + c / K_I)**2)**2)
# 
# 
# def effective_hill_cred(num_rep, Op, e_AI, K_A, K_I,
                        # mass_frac=0.95):
    # cred_region = np.zeros([2, len(num_rep)])
    # for i, R in enumerate(num_rep):
        Compute the EC50
        # c = EC50(K_A, K_I, e_AI, R, Op)
        Compute the hill
        # e_hill = 2 / (fold_change(c, K_A, K_I, e_AI, R, Op) - fold_change(0, K_A, K_I, e_AI, R, Op)) *\
            # (-(fold_change(c, K_A, K_I, e_AI, R, Op))**2 * R / 5E6 * np.exp(-Op) *
            #  2 * c * np.exp(-e_AI) * (1 / K_A * (1 + c / K_A) * (1 + c / K_I)**2 - 1 / K_I *
                                    #   (1 + c / K_A)**2 * (1 + c / K_I)) / ((1 + c / K_A)**2 + np.exp(-e_AI) * (1 + c / K_I)**2)**2)
        # cred_region[:, i] = mwc.hpd(e_hill, mass_frac)
# 
    # return cred_region
# 
# 
# rep_range = np.logspace(0, 4, 200)
# ka_ki = np.exp(-ea) / np.exp(-ei)
# en_colors = sns.color_palette('viridis', n_colors=len(operators))
# titles = ['leakiness', 'saturation', 'dynamic range',
        #   'EC50 ($\mu$M)', 'effective Hill coefficient']
# for i, op in enumerate(operators):
    Compute the properties
    # leak = leakiness(rep_range, energies[op], ep_ai=4.5)
    # sat = saturation(rep_range, energies[op], 4.5, np.exp(-ea) / np.exp(-ei))
    # dyn_rng = dyn_range(rep_range, energies[op], ka_ki)
    # ec50 = EC50(np.exp(-ea), np.exp(-ei), 4.5, rep_range, energies[op])
    # e_hill = effective_Hill(np.exp(-ea), np.exp(-ei),
                            # 4.5, rep_range, energies[op])
# 
    # ax[3].plot(rep_range, leak, color=en_colors[i], label=energies[op])
    # ax[4].plot(rep_range, sat, color=en_colors[i], label=energies[op])
    # ax[5].plot(rep_range, dyn_rng, color=en_colors[i], label=energies[op])
    # ax[6].plot(rep_range, ec50 / 1E6, color=en_colors[i])
    # ax[7].plot(rep_range, e_hill, color=en_colors[i])
    # ax[i + 3].set_xlabel('repressors per cell', fontsize=12)
    # ax[i + 3].set_ylabel(titles[i], fontsize=12)
# 
    Plot the credible regions
    # sat_cred = saturation_cred_region(
        # rep_range, energies[op], 4.5, ka_fc, ki_fc)
    # dyn_cred = dyn_cred_region(rep_range,
                            #    ka_fc, ki_fc, epsilon=4.5,
                            #    ep_r=energies[op])
    # ec50_cred = ec50_cred_region(rep_range, energies[op], 4.5, ka_fc,
                                #  ki_fc, mass_frac=0.95)
    # hill_cred = effective_hill_cred(
        # rep_range, energies[op], 4.5, ka_fc, ki_fc, mass_frac=0.95)
    # ax[5].fill_between(rep_range, dyn_cred[0, :], dyn_cred[1, :],
                    #    alpha=0.3, color=en_colors[i])
    # ax[4].fill_between(rep_range, sat_cred[0, :], sat_cred[1, :],
                    #    alpha=0.3, color=en_colors[i])
    # ax[6].fill_between(rep_range, ec50_cred[0, :] / 1E6, ec50_cred[1, :] / 1E6,
                    #    alpha=0.3, color=en_colors[i])
    # ax[7].fill_between(rep_range, hill_cred[0, :], hill_cred[1, :],
                    #    alpha=0.3, color=en_colors[i])
# 
    # ax[i + 3].set_xlim([1, 1E4])
    
# 
# ax[6].set_xlim([1, 1E4])
# ax[7].set_xlim([1, 1E4])
# ax[6].set_ylabel('$[EC_{50}]\,\,$(M)', fontsize=12)
# ax[7].set_ylabel('effective Hill coefficient', fontsize=12)
# leg_1 = ax[0].legend(loc='upper left', title='rep. / cell',
                    #  fontsize=8, handlelength=1)
# leg_2 = ax[3].legend(title='   binding\n energy ($k_BT$)',
                    #  loc='lower left', fontsize=8, handlelength=1)
# leg_1.get_title().set_fontsize(8)
# leg_2.get_title().set_fontsize(8)
# ax[3].set_yscale('log')
# ax[6].set_yscale('log')
# ax[6].set_yticks([1E-6, 1E-5, 1E-4])
# ax[7].set_yticks([1.2, 1.4, 1.6, 1.8])
# ax[8].set_axis_off()
# 
# for i in range(3, len(ax)):
    # ax[i].set_xscale('log')
    # ax[i].set_xticks([1, 10, 100, 1000, 1E4])
    # ax[i].set_xlabel('repressors per cell', fontsize=12)
# 

plt.figtext(0.01, 0.96, '(C)', fontsize=8)
plt.figtext(0.33, 0.96, '(D)', fontsize=8)
plt.figtext(0.64, 0.96, '(E)', fontsize=8)
plt.figtext(0.01, 0.65, '(F)', fontsize=8)
plt.figtext(0.34, 0.65, '(G)', fontsize=8)
plt.figtext(0.64, 0.65, '(H)', fontsize=8)
plt.figtext(0.01, 0.32, '(I)', fontsize=8)
plt.figtext(0.34, 0.32, '(J)', fontsize=8)
# 
# 
# plt.tight_layout()
# plt.savefig('../../figures/main_figs/fig5_curves.svg', bbox_inches='tight')
# 
%%
# 

# lab = ['$K_A\,\,(\mu\mathrm{M})$', '$K_I\,\,(\mu\mathrm{M})$']
# ka_ki_df = pd.DataFrame(np.array([ka_fc, ki_fc]).T, columns=lab)
# inds = np.arange(0, len(ka_fc), 1)
# np.random.seed(666)
# 

# 

# plt.close('all')
# g = sns.JointGrid(lab[0], lab[1], ka_ki_df, xlim=(100, 200), ylim=(0.45, 0.625), space=0.05,
                #   size=2)
# 
# g.ax_joint.plot(ka_fc, ki_fc, '.', color='#937D69', ms=2,
                # alpha=0.05, rasterized=True, zorder=1)
# g.plot_joint(sns.kdeplot, cmap=sns.cubehelix_palette(n_colors=10, as_cmap=True, reverse=True), zorder=10, linewidth=1, n_levels=5, shade=True, alpha=0.5, shade_lowest=False,
            #  )
# 
# 

# ind = np.where(gauss_flatlnprobability == gauss_flatlnprobability.max())[0]
# ka_mode = ka_fc[ind][0]
# ki_mode = ki_fc[ind][0]
# ka_cred = mwc.hpd(ka_fc, mass_frac=0.95)
# ki_cred = mwc.hpd(ki_fc, mass_frac=0.95)
# 
# g.ax_marg_y.plot(6, ki_mode, 'o', color=colors[4])
# g.ax_marg_y.vlines(6, ki_cred[0], ki_cred[1], color=colors[4])
# g.ax_marg_x.plot(ka_mode, 0.015, 'o', color=colors[4])
# g.ax_marg_x.hlines(0.015, ka_cred[0], ka_cred[1], color=colors[4])
# 
# g.ax_joint.set_xlim([100, 230])
# g.ax_joint.set_ylim([0.45, 0.65])
# g.plot_marginals(sns.kdeplot, shade=True,
                #  color=colors[4], zorder=1, linewidth=1)
# 
# Plot the mode and HPD for each marginal distribution
# g.fig.set_figwidth(5.75)
# g.fig.set_figheight(3.25)

# 
# 

# 

# plt.savefig('../../figures/main_figs/fig5_ka_ki_posterior.pdf',
            # bbox_inches='tight')


# %%
