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
import altair as alt
import seaborn as sns
import scipy.stats
colors, palette = phd.viz.altair_theme()
DPI = 227
data = pd.read_csv('../../data/ch2_induction/RazoMejia_2018.csv', comment='#')

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

# Define the operators and their respective energies
operators = ['O1', 'O2', 'O3']
energies = {'O1': -15.3, 'O2': -13.9, 'O3': -9.7}


#%%
# ##############################################################################
# SAMPLING JOINTPLOT
# ##############################################################################
sampling_df = pd.DataFrame(np.array([ka_fc, ki_fc, lnprob]).T, 
                          columns=['ka', 'ki', 'logprob'])
sampling_df.sort_values('logprob', inplace=True)
sample_base = alt.Chart(sampling_df, width=3 * DPI, height=2 * DPI)
sampling_points = sample_base.mark_point(size=0.5).encode(
                x=alt.X('ki:Q', axis=alt.Axis(title='KI [µM]'),
                        scale=alt.Scale(domain=[0.4, 0.65])),
                y=alt.Y('ka:Q', axis=alt.Axis(title='KA [µM]'),
                        scale=alt.Scale(domain=[90, 230])),
                fill=alt.value(colors['black']),
                strokeWidth=alt.value(0),
                opacity=alt.value(0.5))

sampling_points.save('../figs/fig5_sampling_points.svg')
sampling_points
#%%
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

#  Iterate through each operator and compute the properties
prop_df = pd.DataFrame([])
for op, op_en in energies.items():
    for i, r in enumerate(rep_range):
        arch = phd.thermo.SimpleRepression(r, op_en, ka=ka_fc, ki=ki_fc, ep_ai=4.5,
                                           effector_conc=0).compute_properties()
        for prop, val in arch.items():
            if prop == 'leakiness':
                val_min = val 
                val_max = val 
            else:
                val_min, val_max = phd.stats.compute_hpd(val, 0.95)
            prop_df = prop_df.append({'val_min':val_min,
                              'val_max':val_max,
                              'property': prop,
                              'repressors': r,
                              'operator': op,
                              'binding_energy':op_en}, ignore_index=True)
        

#%%

# Generate the plot of the fit strain.
fit_strain = data[(data['operator']=='O2') & (data['repressors']==130) & 
                  (data['IPTG_uM']>0)]

fit_strain_base = alt.Chart(fit_strain).transform_aggregate(
                            mean_fc='mean(fold_change_A)',
                            sem_fc='stderr(fold_change_A)',
                            groupby=['IPTG_uM']
                            ).transform_calculate(
                            ymin='datum.mean_fc - datum.sem_fc',
                            ymax='datum.mean_fc + datum.sem_fc' 
                            )
fit_strain_plot = fit_strain_base.mark_point(size=30).encode(
                            x=alt.X('IPTG_uM:Q', scale=alt.Scale(type='log')),
                            y=alt.Y('mean_fc:Q', axis=alt.Axis(title='fold-change')),
                            stroke=alt.value(colors['dark_orange']),
                            strokeWidth=alt.value(.75),
                            fill=alt.value('white'))
fit_strain_errors = fit_strain_base.mark_rule().encode(
                            x='IPTG_uM:Q',
                            y='ymin:Q',
                            y2='ymax:Q',
                            strokeWidth=alt.value(2),
                            color=alt.value(rep_colors[260])) 
                    

# Generate the plot for each operator
fc_fill = alt.hconcat()
for g, d in fc_df.groupby('operator'):
    _fc_fill = alt.Chart(data=d[d['IPTGuM']>0], width=200, 
                        height=200).mark_area().encode(
                        x = alt.X('IPTGuM:Q', axis=alt.Axis(title='IPTG [µM]',
                                tickCount=4), scale= alt.Scale(type='log')),   
                        y = alt.Y('fc_min:Q', axis=alt.Axis(title='fold-change'),
                                scale= alt.Scale(domain=[-0.05, 1.1])),
                        y2='fc_max:Q', 
                        fill=alt.Color('repressors:O', scale=alt.Scale(
                                range=list(rep_colors.values())), 
                                legend=alt.Legend(title='repressors per cell')),
                        opacity= alt.value(0.75),
                        strokeWidth = alt.value(0.5),
                        stroke=alt.Color('repressors:O',scale=alt.Scale(
                                range=list(rep_colors.values())))

                    ).properties(
                        title=f'operator {g}'
                    )
            
    if g == 'O2':
        _fc_fill += fit_strain_plot + fit_strain_errors
    fc_fill |= _fc_fill
fc_fill.save('../figs/fig5_induction_plots.svg')

#%%


pred_leak = alt.Chart(prop_df[prop_df['property']=='leakiness'], 
                     width=200, height=200).mark_line(size=3).encode(
                     x = alt.X('repressors:Q', 
                                axis=alt.Axis(title='repressors per cell',
                                              tickCount=4),
                                             scale=alt.Scale(type='log')),
                    y=alt.Y('val_min:Q', axis=alt.Axis(title='leakiness',
                                                       tickCount=4),
                                        scale=alt.Scale(type='log')),
                    opacity=alt.value(0.5),
                    color=alt.Color(field='operator',
                                    type='ordinal',
                                    scale=alt.Scale(scheme='viridis'))
                    )
    
pred_sat = alt.Chart(prop_df[prop_df['property']=='saturation'],
                     width=200, height=200).mark_area().encode(
                     x = alt.X('repressors:Q', 
                              axis=alt.Axis(title='repressors per cell',
                                            tickCount=4),
                                            scale=alt.Scale(type='log')),
                     y=alt.Y('val_min:Q', axis=alt.Axis(title='saturation',
                                                       tickCount=4)), 
                    y2='val_max:Q',
                    opacity=alt.value(0.5),
                    color=alt.Color(field='operator',
                                    type='ordinal',
                                    scale=alt.Scale(scheme='viridis'))


                     )

pred_dyn_rng = alt.Chart(prop_df[prop_df['property']=='dynamic_range'],
                     width=200, height=200).mark_area().encode(
                     x = alt.X('repressors:Q', 
                              axis=alt.Axis(title='repressors per cell',
                                            tickCount=4),
                                            scale=alt.Scale(type='log')),
                     y=alt.Y('val_min:Q', axis=alt.Axis(title='dynamic range',
                                                       tickCount=4)), 
                    y2='val_max:Q',
                    opacity=alt.value(0.5),
                    color=alt.Color(field='operator',
                                    type='ordinal',
                                    scale=alt.Scale(scheme='viridis'))


                     )

pred_ec50 = alt.Chart(prop_df[prop_df['property']=='EC50'],
                     width=200, height=200).mark_area().encode(
                     x = alt.X('repressors:Q', 
                              axis=alt.Axis(title='repressors per cell',
                                            tickCount=4),
                                            scale=alt.Scale(type='log')),
                     y=alt.Y('val_min:Q', axis=alt.Axis(title='EC\u2085\u2080 [µM]',
                                                       tickCount=4),
                                          scale=alt.Scale(type='log')), 
                    y2='val_max:Q',
                    opacity=alt.value(0.5),
                    color=alt.Color(field='operator',
                                    type='ordinal',
                                    scale=alt.Scale(scheme='viridis'))


                     )

pred_e_hill = alt.Chart(prop_df[prop_df['property']=='effective_hill'],
                     width=200, height=200).mark_area().encode(
                     x = alt.X('repressors:Q', 
                              axis=alt.Axis(title='repressors per cell',
                                            tickCount=4),
                                            scale=alt.Scale(type='log')),
                     y=alt.Y('val_min:Q', axis=alt.Axis(title='effective Hill coefficient',
                                                       tickCount=4),
                                          scale=alt.Scale(domain=[1.1, 1.9])), 
                    y2='val_max:Q',
                    opacity=alt.value(0.5),
                    color=alt.Color(field='operator',
                                    type='ordinal',
                                    scale=alt.Scale(scheme='viridis'))


                     )



prop_row1 = alt.hconcat(pred_leak | pred_sat | pred_dyn_rng)
prop_row2 = alt.hconcat(pred_ec50 | pred_e_hill)
fig_col = alt.vconcat(fc_fill, prop_row1, prop_row2)
fig_col
#%%
row = alt.hconcat()
for g, d in fc_df.groupby('operator'):
    base = alt.Chart(d, width=100, height=150)
    fill = o2_base.mark_area().encode(
                    x=alt.X(field='IPTGuM', 
                            axis=alt.Axis(title='IPTG [µM]',
                                          tickCount=6),
                            scale=alt.Scale(type='log'),),
                    y=alt.Y(field='fc_min', 
                            axis=alt.Axis(title='fold-change'),
                            scale=alt.Scale(type='log')),
                    y2='fc_max',
                    color=alt.Color(field='repressors', 
                                    type='ordinal',
                                    scale=alt.Scale(range=list(rep_colors.values()))),
                    opacity=alt.value(0.5)

            ).properties(
                title=f'operator {op}'
            )
    bases.append(base)
    row |= fill

row

#%%
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
