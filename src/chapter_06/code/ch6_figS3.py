#%%
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import fsolve
import phd.viz
colors, palette = phd.viz.phd_style()

# Define functions to be used in figure
def pact(IPTG, K_A, K_I, e_AI):
    '''
    Computes the probability that a repressor is active
    Parameters
    ----------
    IPTG : array-like
        Array of IPTG concentrations in uM
    K_A : float
        Dissociation constant for active repressor
    K_I : float
        Dissociation constant for inactive repressor
    e_AI : float
        Energetic difference between the active and inactive state
    Returns
    -------
    probability that repressor is active
    '''
    pact = (1 + IPTG * 1 / K_A)**2 / \
        (((1 + IPTG * 1 / K_A))**2 + np.exp(-e_AI) * (1 + IPTG * 1 / K_I)**2)
    return pact


def fugacity(IPTG, R, Ns, e_s, K_A=139E-6, K_I=0.53E-6, e_AI=4.5, Nc=0, e_c=0):
    '''
    Solves for the fugacity of simple repression with
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
    fugacity at each IPTG concentration
    '''
    NNS = 4.6E6
    lam = []

    def func(x): return -Reff + Ns * (x * np.exp(-e_s)) / (1 + x * np.exp(-e_s)) +\
        NNS * (x) / (1 + x) + \
        Nc * (x * np.exp(-e_c)) / (1 + x * np.exp(-e_c))
    for c in IPTG:
        Reff = R * pact(c, K_A, K_I, e_AI)
        lam.append(fsolve(func, 0))
    return np.array(lam)


def occupancy(lam, e_s):
    '''
    Computes fold-change for simple repression using the fugacity (lam).
    Parameters
    ----------
    lam : float
        fugacity of system as calculated by fugacity()
    e_s : float
        binding energy of specific operator
    Returns
    -------
    fold-change (occupancy)
    '''
    return 1 / (1 + lam * np.exp(-(e_s)))


# Define parameter values
ops = [-15.3, -13.9, -9.7]
op_names = ['O1', 'O2', 'O3']
fig_labels = [['(A)', '(B)', '(C)'], ['(D)', '(E)', '(F)']]
reps = [1740, 1220, 260, 124, 60, 22]
rep_colors = {22: colors['red'], 
              60: colors['brown'],  
              124: colors['green'], 
              260: colors['orange'], 
              1220: colors['purple'], 
              1740: colors['blue']}
Ns = [10, 100]
IPTG = np.logspace(-8, -2, 100)

# Plot figure
fig, ax = plt.subplots(2, 3, sharey=False, figsize=(6, 4))
phd.viz.despine(ax.ravel())
for i, a in enumerate(ax[0]):
    for rep in reps:
        lam_array = fugacity(IPTG, rep, Ns=Ns[0], e_s=ops[i])
        fc = occupancy(lam_array, ops[i])
        _ = a.plot(IPTG*1E6, fc, label=rep, color=rep_colors[rep])
    a.set_xscale('log')
    a.set_ylabel('fold-change')
    a.set_xlabel('IPTG [µM]')
    a.set_ylim(-0.01, 1.1)
    a.set_xlim(1E-2, 1E4)
    a.tick_params(labelsize=6)

    # Define figure text
    phd.viz.titlebox(a, r'%s $\Delta \varepsilon_{RA}= %0.1f\ k_BT$' % (op_names[i], ops[i]),
                     bgcolor='white', color=colors['black'], boxsize='12%', size=8,
                     pad=0.04)



    a.text(-0.28, 1.1, fig_labels[0][i], ha='center', va='center', fontsize=8,
           transform=a.transAxes)

for i, a in enumerate(ax[1]):
    for rep in reps:
        lam_array = fugacity(IPTG, rep, Ns=Ns[1], e_s=ops[i])
        fc = occupancy(lam_array, ops[i])
        _ = a.plot(IPTG*1E6, fc, label=rep, color=rep_colors[rep])
    a.set_xscale('log')
    a.set_ylabel('fold-change', fontsize=8)
    a.set_xlabel('IPTG [µM]', fontsize=8)
    a.set_ylim(-0.01, 1.1)
    a.set_xlim(1E-2, 1E4)
    a.tick_params(labelsize=6)

    # Define figure text
    phd.viz.titlebox(a, r'%s $\Delta \varepsilon_{RA}= %0.1f\ k_BT$' % (op_names[i], ops[i]),
                     bgcolor='white', color=colors['black'], boxsize='12%', size=8,
                     pad=0.04)

    a.text(-0.28, 1.1, fig_labels[1][i], ha='center', va='center', fontsize=8,
           transform=a.transAxes)

# Add legend
leg1 = ax[0][0].legend(title='rep. per cell', loc='upper left', fontsize=6)
leg1.get_title().set_fontsize(6)
plt.subplots_adjust(wspace=0.35, hspace=0.5)
plt.tight_layout()
plt.savefig('../figs/ch6_figS3.pdf', bbox_inches='tight')
plt.savefig('../figs/ch6_figS3.png', bbox_inches='tight')



# %%
