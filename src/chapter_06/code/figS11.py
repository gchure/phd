
#%%
# For scientific computing
import numpy as np
import pandas as pd
import scipy.special
import phd.viz
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn as sns
colors, palette = phd.viz.phd_style()


# Variability  in fold-change as parameter change.
def fold_change_oo(Ka, Ki, R, era, eai=4.5, Nns=4.6E6):
    '''
    computes the gene expression fold change for a simple repression architecture
    in the limit where the inducer concentration goes to infinity
    Parameters
    ----------
    Ka, Ki : float.
        Dissociation constants of the ligand to the active and inactive state
        of the repressor respectively.
    R : float.
        Mean repressor copy number per cell
    era : float.
        Repressor-DNA binding energy
    eai : float.
        Energy difference between active and inactive state of repressor
    Nns : float.
        Number of non-specific binding sites.
    Returns
    -------
    fold-change
    '''
    return (1 + 1 / (1 + np.exp(-eai) * (Ka / Ki)**2) * R / Nns *
            np.exp(-era))**-1


# Let us now define the numerical values for all the needed parameters
era_num = np.array([-15.3, -13.9, -9.7])  # kBT
Ka_num = 139.96  # µM
Ki_num = 0.54  # µM


# Let's now plot the change in fold-change as $K_A$ and $K_I$ vary for
# different energies and repressor copy numbers.

# Factor by which the Ka and Ki are varied
factor = 2

Ka_array = np.logspace(np.log10(Ka_num / factor),
                       np.log10(Ka_num * factor), 100)
Ki_array = np.logspace(np.log10(Ki_num / factor),
                       np.log10(Ki_num * factor), 100)


# Initialize plot
fig, ax = plt.subplots(2, 2, figsize=(6, 4))
ax = ax.ravel()
phd.viz.despine(ax)

# Loop through binding energies
_colors = sns.color_palette('viridis', n_colors=4)
rep = 260  # repressors per cell
for i, eRA in enumerate(era_num):
    # compute the ∆fold-change_Ka
    delta_fc = fold_change_oo(Ka_num, Ki_num, rep, eRA) - \
        fold_change_oo(Ka_array, Ki_num, rep, eRA)

    ax[0].plot(np.log10(Ka_array / Ka_num), delta_fc,
               label=r'{:.1f}'.format(eRA), color=_colors[i])

    # compute the ∆fold-change_KI
    delta_fc = fold_change_oo(Ka_num, Ki_num, rep, eRA) - \
        fold_change_oo(Ka_num, Ki_array, rep, eRA)

    ax[1].plot(np.log10(Ki_array / Ki_num), delta_fc,
               label=r'{:.1f}'.format(eRA), color=_colors[i])

# Format Ka plot
ax[0].set_xlabel(r'$\log_{10} \frac{K_A}{K_A^{fit}}$')
ax[0].set_ylabel(r'$\Delta$fold-change$_{K_A}$')

# Format Ki plot
ax[1].set_xlabel(r'$\log_{10} \frac{K_I}{K_I^{fit}}$')
ax[1].set_ylabel(r'$\Delta$fold-change$_{K_I}$')
leg = ax[1].legend(loc='center left', title=r'$\Delta\varepsilon_{RA}$ ($k_BT$)',
             ncol=1, fontsize=6, bbox_to_anchor=(1, 0.5))
leg.get_title().set_fontsize(6)
# Set the colors for the strains


rep_colors = {22: colors['red'], 
              60: colors['brown'],  
              124: colors['green'], 
              260: colors['orange'], 
              1220: colors['purple'], 
              1740: colors['blue']}
eRA = -15.3
for i, rep in enumerate(rep_colors.keys()):
    # compute the ∆fold-change_Ka
    delta_fc = fold_change_oo(Ka_num, Ki_num, rep, eRA) - \
        fold_change_oo(Ka_array, Ki_num, rep, eRA)

    ax[2].plot(np.log10(Ka_array / Ka_num), delta_fc,
               label=str(rep), color=rep_colors[rep])

    # compute the ∆fold-change_KI
    delta_fc = fold_change_oo(Ka_num, Ki_num, rep, eRA) - \
        fold_change_oo(Ka_num, Ki_array, rep, eRA)

    ax[3].plot(np.log10(Ki_array / Ki_num), delta_fc,
               label=str(rep), color=rep_colors[rep])

# Format Ka plot
ax[2].set_xlabel(r'$\log_{10} \frac{K_A}{K_A^{fit}}$')
ax[2].set_ylabel(r'$\Delta$fold-change$_{K_A}$')


# # Format Ki plot
ax[3].set_xlabel(r'$\log_{10} \frac{K_I}{K_I^{fit}}$')
ax[3].set_ylabel(r'$\Delta$fold-change$_{K_I}$')
ax[3].margins(0.01)
leg = ax[3].legend(loc='center left', title=r'repressors / cell', ncol=1,
             fontsize=6, bbox_to_anchor=(1, 0.5))
leg.get_title().set_fontsize(6)

# Label plot
plt.figtext(0.01, .95, '(A)', fontsize=8)
plt.figtext(0.45, .95, '(B)', fontsize=8)
plt.figtext(0.01, .48, '(C)', fontsize=8)
plt.figtext(0.45, .48, '(D)', fontsize=8)
plt.tight_layout()
plt.savefig('../figs/ch6_figS11.pdf', bbox_inches='tight')
plt.savefig('../figs/ch6_figS11.png', bbox_inches='tight')

# %%
