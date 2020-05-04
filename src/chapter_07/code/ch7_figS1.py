#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import phd.viz
import phd.thermo
colors, palette = phd.viz.phd_style()

# %%
# Define the dependent variables
c_range = np.logspace(-2, 4, 500)
kaki_range = np.logspace(-3, 3, 500)
# Compute the components of the delta F as a function of c
ka = 200
ki = 1
ep_ai = 4.5
theta1 = 0.1
theta2 = 4
theta3 = 4
theta4 = 2
theta5 = 0.01
theta6=0.1

def dF_dc(ka, ki, ep_ai, c):
    numer = (ka - ki) * (c + ki)
    denom = (c + ka) * ((c + ka)**2 * ki**2 + np.exp(-ep_ai) * ka**2 * (c + ki)**2)
    return numer / denom

def F(ka, ki, ep_ai, c):
    numer = (1 + c/ka)**2
    denom = numer + np.exp(-ep_ai) * (1 + c/ki)**2
    return - np.log(numer/denom)

# %% 
fig, ax = plt.subplots(2, 2, sharex=True, figsize=(6, 4))
phd.viz.despine(ax.ravel())
for a in ax.ravel():
    a.set_xscale('log')
    a.set_xlabel('$c \, / \, K_A^\mathrm{(wt)}$', fontsize=8)

wt = 2 * np.exp(-ep_ai) * ka**2 * dF_dc(ka ,ki, ep_ai, c_range)
mut1 = 2 * np.exp(-ep_ai) * ka**2 * dF_dc(ka*theta1 ,ki*theta1, ep_ai, c_range)
mut2 = 2 * np.exp(-ep_ai) * ka**2 * dF_dc(ka*theta2 ,ki*theta2, ep_ai, c_range)
mut3 = 2 * np.exp(-ep_ai) * ka**2 * dF_dc(ka*theta3 ,ki*theta4, ep_ai, c_range)
mut4 = 2 * np.exp(-ep_ai) * ka**2 * dF_dc(ka*theta5 ,ki*theta6, ep_ai, c_range)

Fwt = F(ka, ki, ep_ai, c_range)
Fmut1 = F(ka*theta1, ki*theta1, ep_ai, c_range)
Fmut2 = F(ka*theta2, ki*theta2, ep_ai, c_range)
Fmut3 = F(ka*theta3, ki*theta4, ep_ai, c_range)
Fmut4 = F(ka*theta5, ki*theta6, ep_ai, c_range)

ax[0, 0].plot(c_range / ka, Fwt, label='wild type', color=colors['black'])
ax[0, 0].plot(c_range / ka, Fmut1, label=r'$\theta=0.1$', color=colors['green'])
ax[0, 0].plot(c_range / ka, Fmut2, label=r'$\theta=4$', color=colors['blue'])
ax[0, 1].plot(c_range / ka, Fwt - Fwt, color=colors['black'])
ax[0, 1].plot(c_range / ka, Fmut1 - Fwt, color=colors['green'])
ax[0, 1].plot(c_range / ka, Fmut2 - Fwt, color=colors['blue'])

ax[1, 0].plot(c_range / ka, Fwt, label='wild type', color=colors['black'])
ax[1, 0].plot(c_range / ka, Fmut3, label=r'2 $K_A^\mathrm{(wt)} / K_I^\mathrm{(wt)}$', color=colors['green'])
ax[1, 0].plot(c_range / ka, Fmut4, label=r'0.1 $K_A^\mathrm{(wt)} / K_I^\mathrm{(wt)}$', color=colors['blue'])

ax[1, 1].plot(c_range / ka, Fwt - Fwt, color=colors['black'])
ax[1, 1].plot(c_range / ka, Fmut3 - Fwt, color=colors['green'])
ax[1, 1].plot(c_range / ka, Fmut4 - Fwt, color=colors['blue'])

for i in range(2):
    ax[i, 0].set_ylabel('$F\, / \, k_BT$')
    ax[i, 1].set_ylabel('$\Delta F \, / \, k_BT$')

phd.viz.titlebox(ax[1, 0], 'different $K_A/K_I$', bgcolor='white', color=colors['black'], pad=0.05, boxsize='12%')
phd.viz.titlebox(ax[1, 1], 'different $K_A/K_I$', bgcolor='white', color=colors['black'], pad=0.05, boxsize='12%')
phd.viz.titlebox(ax[0, 0], 'same $K_A/K_I$', bgcolor='white', color=colors['black'], pad=0.05, boxsize='12%')
phd.viz.titlebox(ax[0, 1], 'same $K_A/K_I$', bgcolor='white', color=colors['black'], pad=0.05, boxsize='12%')
ax[0, 0].legend(handlelength=0.5, fontsize=6)
ax[1, 0].legend(handlelength=0.5, fontsize=6)

# Add panel labels
fig.text(0.05, 0.87, '(A)', fontsize=8)
fig.text(0.05, 0.45, '(B)', fontsize=8)
plt.subplots_adjust(wspace=0.4)
plt.savefig('../figs/ch7_figS1.pdf', bbbox_inches='tight')
plt.savefig('../figs/ch7_figS1.png', bbbox_inches='tight')

#%%
