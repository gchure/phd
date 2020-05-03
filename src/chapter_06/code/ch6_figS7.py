#%%
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import fsolve
import phd.viz
colors, palette = phd.viz.phd_style()

# Define functions to be used for figures
def fugacity_leakiness(R, Ns, e_s, e_AI=4.5, Nc=0, e_c=0):
    '''
    Solves for the leakiness of a simple repression construct with
    multiple promoter copies (Ns, with energy e_s) or competitor sites
    (Nc, with energy e_c).
    Parameters
    ----------
    R : float
        Number of repressors per cell
    e_AI : float
        Energetic difference between the active and inactive state
    Ns : float
        Number of specific operators available for repressor binding
    Nc : float
        Number of competitor operators available for repressor binding
    e_s : float
        Binding energy between specific operator and repressor as inferred in
        Garcia 2011
    e_c : float
        Binding energy between competitor operator and repressor
    e_AI : float
        Energetic difference between the active and inactive state
    Returns
    -------
    leakiness
    '''
    NNS = 4.6E6
    p_A = 1 / (1 + np.exp(-e_AI))
    Reff = R * p_A
    leakiness = []
    for R in R:
        Reff = R * p_A

        def func(x): return -Reff + Ns * (x * np.exp(-e_s)) / (1 + x * np.exp(-e_s)) +\
            NNS * (x) / (1 + x) + \
            Nc * (x * np.exp(-e_c)) / (1 + x * np.exp(-e_c))
        lam = fsolve(func, 0)
        leakiness.append(1 / (1 + lam * np.exp(-(e_s))))
    return np.array(leakiness)


def fugacity_saturation(R, Ns, e_s, K_A=139E-6, K_I=0.53E-6, e_AI=4.5, Nc=0,
                        e_c=0):
    '''
    Solves for the saturation of a simple repression construct with
    multiple promoter copies (Ns, with energy e_s) or competitor sites
    (Nc, with energy e_c).
    Parameters
    ----------
    R : float
        Number of repressors per cell
    e_AI : float
        Energetic difference between the active and inactive state
    Ns : float
        Number of specific operators available for repressor binding
    Nc : float
        Number of competitor operators available for repressor binding
    e_s : float
        Binding energy between specific operator and repressor as inferred in
        Garcia 2011
    e_c : float
        Binding energy between competitor operator and repressor
    K_A : float
        Dissociation constant for active repressor
    K_I : float
        Dissociation constant for inactive repressor
    e_AI : float
        Energetic difference between the active and inactive state
    Returns
    -------
    saturation
    '''
    NNS = 4.6E6
    p_A = 1 / (1 + np.exp(-e_AI) * (K_A / K_I)**2)
    saturation = []
    for R in R:
        Reff = R * p_A

        def func(x): return -Reff + Ns * (x * np.exp(-e_s)) / (1 + x * np.exp(-e_s)) +\
            NNS * (x) / (1 + x) + \
            Nc * (x * np.exp(-e_c)) / (1 + x * np.exp(-e_c))
        lam = fsolve(func, 0)
        saturation.append(1 / (1 + lam * np.exp(-(e_s))))
    return np.array(saturation)


def fugacity_dynamic_range(R, Ns, e_s, K_A=139E-6, K_I=0.53E-6, e_AI=4.5, Nc=0,
                           e_c=0):
    '''
    Solves for the dynamic range of a simple repression construct with
    multiple promoter copies (Ns, with energy e_s) or competitor sites
    (Nc, with energy e_c).
    Parameters
    ----------
    R : float
        Number of repressors per cell
    e_AI : float
        Energetic difference between the active and inactive state
    Ns : float
        Number of specific operators available for repressor binding
    Nc : float
        Number of competitor operators available for repressor binding
    e_s : float
        Binding energy between specific operator and repressor as inferred in
        Garcia 2011
    e_c : float
        Binding energy between competitor operator and repressor
    K_A : float
        Dissociation constant for active repressor
    K_I : float
        Dissociation constant for inactive repressor
    e_AI : float
        Energetic difference between the active and inactive state
    Returns
    -------
    dynamic range
    '''
    return fugacity_saturation(R, Ns, e_s, K_A, K_I, e_AI, Nc, e_c) -\
        fugacity_leakiness(R, Ns, e_s, e_AI, Nc, e_c)


# Define parameters
ops = [-15.3, -13.9, -9.7]
op_names = ['O1', 'O2', 'O3']
fig_labels = [['(A)', '(B)', '(C)'], ['(D)', '(E)', '(F)']]
y_labels = ['leakiness', 'saturation', 'dynamic range']
reps = np.logspace(0, 3, 100)
Nc = [10, 100]
e_c = -17.0

# Plot figure
op_colors = sns.color_palette('viridis', n_colors=len(ops))

fig, ax = plt.subplots(ncols=3, nrows=2, sharey=False, figsize=(6, 4))
phd.viz.despine(ax.ravel())
for i, a in enumerate(ax):
    for j, op in enumerate(ops):
        # Produce plots
        a[j].axvline(Nc[i], ls='-', color='white', alpha=0.75, lw=3)
        a[0].plot(reps, fugacity_leakiness(reps, 1, op, Nc=Nc[i], e_c=e_c),
                  color=op_colors[j], label=ops[j])

        a[1].plot(reps, fugacity_saturation(reps, 1, op, Nc=Nc[i], e_c=e_c),
                  color=op_colors[j])

        a[2].plot(reps, fugacity_dynamic_range(reps, 1, op, Nc=Nc[i], e_c=e_c),
                  color=op_colors[j])

        # Format axes
        a[0].set_yscale('log')
        a[j].set_xscale('log')
        a[j].set_xlabel('repressors per cell')
        a[j].set_ylabel(y_labels[j])
        a[j].set_xlim(1E0, 1E3)

        a[0].set_ylim(9E-4, 1.5E0)
        a[1].set_ylim(-0.015, 1.1)
        a[2].set_ylim(-0.015, 1.1)

# Add figure text
        for k in range(3):
            phd.viz.titlebox(a[k], r'$N_C$= %s' % str(Nc[i]), size=8, 
                            color=colors['black'], bgcolor='white', pad=0.05,
                            boxsize='12%')
        a[j].text(-0.38, 1.05, fig_labels[i][j], transform=a[j].transAxes,
                  fontsize=8)

# Add legend
leg = ax[1][0].legend(title='binding energy \n' + r'$\Delta \varepsilon_{AI}\ [k_BT]$', fontsize=6, loc='lower left',
                      labelspacing=0.1, handlelength=1.5)
plt.setp(leg.get_title(), multialignment='center')
leg.get_title().set_fontsize(6)
plt.tight_layout()
plt.savefig('../figs/ch6_figS7.pdf', bbox_inches='tight')
plt.savefig('../figs/ch6_figS7.png', bbox_inches='tight')


# %%
