#%%
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import phd.viz
colors, palette = phd.viz.phd_style()

# Define the fold-change functions and the p_act.
def p_act(c_range, k_a, k_i, ep_ai=5, n_sites=int(2)):
    """
    Computes the probability of finding an active transcription factor
    at a given range of ligand concentrations.

    Parameters
    ----------
    c_range : 1d-array
        Range of ligand concentrations over which to compute the probability
        of finding an active transcription factor. This should be in the same
        units as k_a and k_i.
    k_a, k_i : floats
        Ligand dissociation constants for the active and inactive transcription
        factors respectively. These should be in the same units as c_range.
    ep_ai : float
        Difference in energy between the active and inactive states of the
        transcription factor. The default value is -4.5 kBT.
    n_sites : int
        Number of cooperative ligand binding sites on the transcription
        factor. The default value is 2.

    Returns
    -------
    p_act : 1d-array
        The probability of finding an active transcription factor at each
        ligand concentration specified in c_range.

    Raises
    ------
    TypeError:
        This exception is thrown if the provided number of cooperative ligand
        binding sites, n_sites, is not an integer.
    """
    if type(n_sites) is not int:
        raise TypeError('n_sites must be an integer.')

    numer = (1 + c_range / k_a)**n_sites
    denom = numer + np.exp(-ep_ai) * (1 + c_range / k_i)**n_sites
    return numer / denom


def fc_repression(c_range, num_rep,  k_a, k_i, ep_r, ep_ai=5,
                  n_sites=int(2), n_ns=4.6E6):
    """
    Computes the fold-change for a simple repression motif over a range of
    ligand concentrations.

    Parameters
    ----------
    c_range : 1d-array
        Range of ligand concentrations over which to compute the fold-change.
    num_rep : int
        Number of repressors in the cell.
    k_a, k_i : floats
        Ligand dissociation constants for the active and inactive transcription
        factor respectively. These should be in the same units as c_range.
    ep_r : float
        Binding energy of the  transcription factor to the DNA. This should be
        in units of kBT.
    ep_ai : float
        Difference in energy between the active an inactive state of the
        transcription factor. Default value is 4.5 kBT.
    n_sites : int
        Number of cooperative ligand binding sites on the transcription factors
    n_ns : int
        Number of nonspecific binding sites available to teh transcription factor. Default value is 4.6E6, the length of the E. coli genome in
        basepairs.

    Returns
    -------
    fold_change : 1d-array
        The fold-change in gene expression for the simple repression motif
        at each concentration of ligand given in c_range.

    Raises
    ------
    TypeError:
        This is thrown when the provided n_sites variable is not an integer.
        Thrown by the function p_act.
    """
    # Compute the MWC probability.
    mwc_term = p_act(c_range, k_a, k_i, n_sites=n_sites, ep_ai=ep_ai)

    # Compute and return the foldchange.
    repression = 1 + mwc_term * (num_rep / n_ns) * np.exp(-ep_r)
    return 1 / repression


def fc_activation(c_range, num_act, k_a, k_i, ep_a, ep_ap, ep_ai=4.5,
                  n_ns=4.6E6, n_sites=int(2)):
    """
    Computes the fold-change for a simple repression motif over a range of
    ligand concentrations.

    Parameters
    ----------
    c_range : 1d-array
        Range of ligand concentrations over which to compute the fold-change.
    num_act : int
        Number of activators in the cell.
    k_a, k_i : floats
        Ligand dissociation constants for the active and inactive transcription
        factor respectively. These should be in the same units as c_range.
    ep_a : float
        Binding energy of the  transcription factor to the DNA. This should be
        in units of kBT.
    ep_ap : float
        Interaction energy between the active transcription factor and the
        polymerase when both are bound to the promoter..
    ep_ai : float
        Difference in energy between the active an inactive state of the
        transcription factor. Default value is 4.5 kBT.
    n_sites : int
        Number of cooperative ligand binding sites on the transcription factors
    n_ns : int
        Number of nonspecific binding sites available to teh transcription factor. Default value is 4.6E6, the length of the E. coli genome in
        basepairs.

    Returns
    -------
    fold_change : 1d-array
        The fold-change in gene expression for the simple repression motif
        at each concentration of ligand given in c_range.

    Raises
    ------
    TypeError:
        This is thrown when the provided n_sites variable is not an integer.
        Thrown by the function p_act.
    """
    # Compute the MWC probability.
    mwc_term = p_act(c_range, k_a, k_i, n_sites=n_sites, ep_ai=ep_ai)

    # Compute and return the fold-change
    numer = 1 + mwc_term * (num_act / n_ns) * np.exp(-(ep_a + ep_ap))
    denom = 1 + mwc_term * (num_act / n_ns) * np.exp(-ep_a)
    return numer / denom

# Set the parameters of interest.
num_tf = [1, 10, 100, 200, 500, 1000]
ep_tf = [-8, -10, -12, -14, -16, -18]
ep_ap = [0, -1, -2, -2.5, -3, -4]
c_range = np.logspace(-9, -2, 500)
ep_ai = 5
k1 = 200E-6
k2 = 0.5E-6

# Define the colors.
tf_colors = sns.color_palette('viridis', n_colors=len(num_tf))
ep_colors = sns.color_palette('plasma', n_colors=len(ep_tf))
act_colors = sns.color_palette('Blues', n_colors=8)

# Set up the figure axis. and add labels.
plt.close('all')
fig, ax = plt.subplots(2, 3, figsize=(6, 4))
# Format the axes.
axes = ax.ravel()
phd.viz.despine(axes)
axes[2].set_axis_off()
for a in axes:
    a.set_xscale('log')
    a.set_xlabel('$c / K_A$')
    a.set_ylabel('fold-change')
    # Set titles.
phd.viz.titlebox(axes[0], 'varying TF copy number', bgcolor='white', 
                color=colors['black'], pad=0.05, boxsize='12%')
phd.viz.titlebox(axes[1], 'varying TF binding energy', bgcolor='white', 
                color=colors['black'], pad=0.05, boxsize='12%')
phd.viz.titlebox(axes[5], 'varying interaction energy', bgcolor='white', 
                color=colors['black'], pad=0.05, boxsize='12%')

# Add the legend.
for i, _ in enumerate(num_tf):
    axes[0].plot([], [], '-', color=tf_colors[i], label=num_tf[i])
    axes[1].plot([], [], '-', color=ep_colors[i], label=ep_tf[i])
    axes[2].plot([], [], '-', color=act_colors[i], label=ep_ap[i])

leg1 = axes[0].legend(title='transcription factors per cell',
               ncol=3, bbox_to_anchor=(3.8, 1.2), fontsize=6)
leg2 = axes[1].legend(title='binding energy [$k_BT$]',
               ncol=3, bbox_to_anchor=(1.15, 0.8), fontsize=6)
leg3 = axes[2].legend(title='interaction energy [$k_BT$]',
               ncol=3, bbox_to_anchor=(0.99, 0.32), fontsize=6)
for l in [leg1, leg2, leg3]:
    l.get_title().set_fontsize(6)
# Plot Varying TF copy number first.
for i, T in enumerate(num_tf):
    # Compute the fold-change curve.
    corepression_fc = fc_repression(c_range, T, k2, k1, ep_tf[2], -ep_ai)
    activation_fc = fc_activation(c_range, T, k2, k1, ep_tf[2], ep_ap[2],
                                  -ep_ai)
    # Plot the result.
    _ = ax[0, 0].plot(c_range / k2, corepression_fc, '-', color=tf_colors[i])
    _ = ax[1, 0].plot(c_range / k2, activation_fc, '-', color=tf_colors[i])

# Plot varying TF binding energy
for i, ep in enumerate(ep_tf):
    # Compute the fold-change curve.
    corepression_fc = fc_repression(c_range, num_tf[3], k2, k1, ep, -ep_ai)
    activation_fc = fc_activation(c_range, num_tf[3], k2, k1, ep, ep_ap[2],
                                  -ep_ai)
    # Plot the result.
    _ = ax[0, 1].plot(c_range / k2, corepression_fc, '-', color=ep_colors[i])
    _ = ax[1, 1].plot(c_range / k2, activation_fc, '-', color=ep_colors[i])

# Plot varying activation energy.
for i, ep in enumerate(ep_ap):
    # Compute the fold-change curve.
    activation_fc = fc_activation(c_range, num_tf[3], k2, k1, ep_tf[2], ep,
                                  -ep_ai)
    # Plot the result.
    _ = ax[1, 2].plot(c_range / k2, activation_fc, '-', color=act_colors[i])

# Add various textual labels.
props = dict(boxstyle='square', facecolor='white')
fig.text(0.02, 0.625, 'corepression', 
         rotation='vertical', fontsize=8, bbox=props)

fig.text(0.02, 0.19, 'inducible activation', 
         rotation='vertical', fontsize=8, bbox=props)

fig.text(0, 0.93, '(A)', fontsize=8)
fig.text(0, 0.5, '(B)', fontsize=8)
plt.subplots_adjust(hspace=0.4, wspace=0.4)
plt.savefig('../figs/ch6_figS20.pdf', bbox_inches='tight')
plt.savefig('../figs/ch6_figS20.png', bbox_inches='tight')

# %%
