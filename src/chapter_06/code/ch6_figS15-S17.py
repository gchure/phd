#%%
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import phd.viz 
import phd.thermo
import phd.stats
colors, palette = phd.viz.phd_style()
constants = phd.thermo.load_constants()

# Define the colors. 
rep_colors = {22: colors['red'], 
              60: colors['brown'],  
              124: colors['green'], 
              260: colors['orange'], 
              1220: colors['purple'], 
              1740: colors['blue']}

# Load the data
data = pd.read_csv('../../data/ch2_induction/RazoMejia_2018.csv', comment='#')
data['repressors'] *= 2
data = data[data['repressors'] > 0]
data = data.groupby(['operator', 'repressors', 
                   'IPTG_uM'])['fold_change_A'].agg(('mean', 'sem')).reset_index()



# %%
FIG_NO = 15
for g, d in data.groupby(['operator']):
    fig, ax = plt.subplots(6, 6, figsize=(6, 6))
    for a in ax.ravel():
        a.set_xscale('log')
        a.set_xlim([1E-2, 1E4])
        a.set_ylim([-0.05, 1.1])
        phd.viz.despine(a)
    for i in range(6):
        ax[i, 0].set_ylabel('fold-change')
        ax[-1, i].set_xlabel('IPTG [µM]')
    for i in range(5):
        for j in range(1, 6):
            ax[i, j].set_xticks([])
            ax[i, j].set_yticks([])
            ax[i, j].spines['left'].set_visible(False)
            ax[i, j].spines['bottom'].set_visible(False)

    for i in range(5):
        ax[i, 0].spines['bottom'].set_visible(False)
        ax[i, 0].set_xticks([])

    for i in range(1, 6):
        ax[-1, i].spines['left'].set_visible(False)
        ax[-1, i].set_yticks([])

    # Label the prediction and comparison strains
    props = dict(boxstyle='square', facecolor='white')
    fig.text(0, 0.45, '$K_A, K_I$ fit strain', fontsize=8,
             bbox=props, rotation='vertical')
    fig.text(0.45, 0.92, 'comparison strain', fontsize=8,
             bbox=props)


    iter = 0
    c_range = np.logspace(-2, 4, 200)
    for _g, _d in d.groupby(['repressors']):
        # Update the axis labels
        ax[0, iter].set_title(f'R = {_g}', fontsize=8)
        ax[iter, 0].text(-0.7, 0.35, f'R = {_g}', fontsize=8, rotation='vertical', 
                         ha='center', transform=ax[iter, 0].transAxes)

        # Load the flatchains. 
        with open(f'../../data/ch2_induction/mcmc/SI_I_{g}_R{_g}.pkl', 'rb') as f:
            unpickler = pickle.Unpickler(f)
            gauss_flatchain = unpickler.load()
            ka = np.exp(-gauss_flatchain[:, 0])[::10]
            ki = np.exp(-gauss_flatchain[:, 1])[::10]


        iter2 = 0
        for __g, __d in d.groupby(['repressors']):
            # Plot the predicted curves. 
            cred_region = np.zeros((2, len(c_range)))
            for i, c in enumerate(c_range):
                arch = phd.thermo.SimpleRepression(R=__g, ep_r=constants[g],
                                               ka=ka, ki=ki, ep_ai=4.5,
                                               effector_conc=c).fold_change()
                cred_region[:, i] = phd.stats.compute_hpd(arch, 0.95)
            ax[iter, iter2].fill_between(c_range, cred_region[0, :], 
                                        cred_region[1, :], alpha=0.5, color=rep_colors[_g])
            # Plot the measured data
            if iter == iter2:
                face = 'w'
                edge = rep_colors[__g]
            else:
                face = rep_colors[__g]
                edge = 'w'
            ax[iter, iter2].errorbar(__d['IPTG_uM'], __d['mean'], __d['sem'], fmt='o',
                                    color=rep_colors[__g], markerfacecolor=face, 
                                    markeredgecolor=edge,
                                    markeredgewidth=0.5, ms=3, linestyle='none')
            iter2 += 1
        iter += 1
    plt.subplots_adjust(wspace=0.05, hspace=0.05)
    plt.savefig(f'../figs/ch6_figS{FIG_NO}.pdf', bbox_inches='tight')
    plt.savefig(f'../figs/ch6_figS{FIG_NO}.png', bbox_inches='tight')
    FIG_NO += 1


# %%
