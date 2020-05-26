#%%
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import phd.viz
import phd.thermo
import seaborn as sns

# Define the modifier
REP = True
EP = True
C = True

# Define reference parameters
R_ref = 100 # in copy number
ep_ref = -12 # in kT
c_ref = 20 # in micromolar
rep_range = np.logspace(0, 4, 200)
ep_range = np.linspace(0, -20, 200)
c_range = np.logspace(-2, 3, 200)
R_shift = rep_range[0]
ep_shift = ep_range[0]
c_shift = c_range[0]
pact = phd.thermo.MWC(ka=200, ki=1, effector_conc=c_range, ep_ai=5).pact()
pact_ref = phd.thermo.MWC(ka=200, ki=1, effector_conc=c_ref, ep_ai=5).pact()
shift_pact = phd.thermo.MWC(ka=200, ki=1, effector_conc=c_shift, ep_ai=5).pact()

# compute the reference bohr. 
ref_state = phd.thermo.SimpleRepression(R_ref, ep_ref, effector_conc=c_ref,
                                ka=200, ki=1, ep_ai=5)

# Compute the initial perturbed bohr
shift_state = phd.thermo.SimpleRepression(R_shift, ep_shift, effector_conc=c_shift,
                                          ka=200, ki=1, ep_ai=5)
# define the collapse function
bohr_range = np.linspace(-15, 15, 200)
scaling_fn = (1 + np.exp(-bohr_range))**-1


# Set up the figure canvas
fig, ax = plt.subplots(2, 3, figsize=(7, 4))
for a in ax.ravel():
    phd.viz.despine(a)

# Format axees and add labels
for i in range(3):
    ax[0, i].set_xlabel('free energy [$k_BT$]')
    ax[0, i].set_xlim([-15, 15])
    ax[0, i].set_ylabel('fold-change')
    ax[1, i].set_ylabel('free energy shift [$k_BT$]')

ax[1, 0].set_ylim([-4.5, 4.5])
ax[1, 0].set_xlabel(r'$\log(c/c_{ref})$')
ax[1, 1].set_xlabel(r'$\log(R/R_{ref})$')
ax[1, 2].set_xlabel(r'$\Delta\varepsilon_{RA}$ - $\Delta\varepsilon_{RA}^\mathrm{(ref)}$ [$k_BT$]')

# Plot the theory curves. 
for i in range(3):
    ax[0, i].plot(bohr_range, scaling_fn, 'k-', lw=0.75)
phd.viz.titlebox(ax[0, 0], 'change in $c$', color=colors['black'])
phd.viz.titlebox(ax[0, 1], 'change in $R$', color=colors['black'])
phd.viz.titlebox(ax[0, 2], r'change in $\Delta\varepsilon_{RA}$', color=colors['black'])
ax[1,  0].plot(np.log(c_range / c_ref), np.log(pact / pact_ref), 'k-', lw=0.75)
ax[1,  1].plot(np.log(rep_range / R_ref), -np.log(rep_range / R_ref), 'k-', lw=0.75)
ax[1,  2].plot(ep_range - ep_ref, ep_range - ep_ref, 'k-', lw=0.75)


# Plot the reference points
shift_bohrs = []
for i in range(3):
    ax[0, i].plot(ref_state.bohr_parameter(), ref_state.fold_change(), 'ko', ms=5, 
                  markeredgecolor=colors['grey'], markeredgewidth=0.5)
    ax[1, i].plot([0], [0], 'ko', ms=5, 
                  markeredgecolor=colors['grey'], markeredgewidth=0.5)


c_bohr, = ax[0, 0].plot(shift_state.bohr_parameter(), shift_state.fold_change(), 'o', ms=5,
              color=colors['blue'], markeredgecolor=colors['grey'])
r_bohr, = ax[0, 1].plot(shift_state.bohr_parameter(), shift_state.fold_change(), 'o', ms=5,
              color=colors['red'], markeredgecolor=colors['grey'])
ep_bohr, = ax[0, 2].plot(shift_state.bohr_parameter(), shift_state.fold_change(),'o',  ms=5,
              color=colors['orange'], markeredgecolor=colors['grey'])

# Define the shift lines
c_bohr_line, = ax[0, 0].plot([ref_state.bohr_parameter(), shift_state.bohr_parameter()], 
                             [ref_state.fold_change(), ref_state.fold_change()], 
                             color=colors['blue'], lw=1, 
                             zorder=-1)

r_bohr_line, = ax[0, 1].plot([ref_state.bohr_parameter(), shift_state.bohr_parameter()], 
                             [ref_state.fold_change(), ref_state.fold_change()], 
                             color=colors['red'], lw=1,
                             zorder=-1)

ep_bohr_line, = ax[0, 2].plot([ref_state.bohr_parameter(), shift_state.bohr_parameter()], 
                             [ref_state.fold_change(), ref_state.fold_change()], 
                             color=colors['orange'], lw=1,
                             zorder=-1)

bohr_lines = [c_bohr_line, r_bohr_line, ep_bohr_line]

c_shift_line, = ax[1, 0].plot([0, 0],  [0, np.log(shift_pact/pact_ref)], color=colors['blue'],
                          lw=1, zorder=-1)
r_shift_line, = ax[1, 1].plot([0, 0],  [0, -np.log(R_shift/R_ref)], color=colors['red'],
                          lw=1, zorder=-1)
ep_shift_line, = ax[1, 2].plot([0, 0],  [0, ep_shift - ep_ref], color=colors['orange'],
                          lw=1, zorder=-1)

shift_lines = [c_shift_line, r_shift_line, ep_shift_line]


# Plot the perturbed points
c_point, = ax[1, 0].plot(np.log(c_shift/c_ref), np.log(shift_pact/pact_ref), 'o', ms=5, 
              color=colors['blue'], markeredgecolor=colors['grey'])
r_point, = ax[1, 1].plot(np.log(R_shift/R_ref), -np.log(R_shift/R_ref), 'o', ms=5, 
              color=colors['red'], markeredgecolor=colors['grey'])
ep_point, = ax[1, 2].plot(ep_shift - ep_ref, ep_shift - ep_ref, 'o', ms=5, 
              color=colors['orange'], markeredgecolor=colors['grey'])

shift_points = [c_point, r_point, ep_point]

# Define the animation functions
def anim_rshift(i):
    # Update the shift plots
    R_shift = rep_range[i]
    new_bohr = -np.log(R_shift/R_ref)
    r_point.set_data([np.log(R_shift/R_ref)], [-np.log(R_shift/R_ref)])

    # Update the bohr plot
    new_arch = phd.thermo.SimpleRepression(R_shift, ep_ref, effector_conc=c_ref,
                                           ka=200, ki=1, ep_ai=500)
    r_bohr.set_data([new_arch.bohr_parameter()], [new_arch.fold_change()])

    # Update the shift indicators
    r_shift_line.set_ydata([0, -np.log(R_shift/R_ref)])
    r_bohr_line.set_xdata([0, new_arch.bohr_parameter()])

def anim_cshift(i):

    # Update the shift plot
    c_shift = c_range[i]
    new_pact = phd.thermo.MWC(ka=200, ki=1, ep_ai=5, effector_conc=c_shift).pact()
    c_point.set_data([np.log(c_shift/c_ref)], [np.log(new_pact/pact_ref)])

    # Update the bohr plot
    new_arch = phd.thermo.SimpleRepression(R_ref, ep_ref, effector_conc=c_shift,
                                           ka=200, ki=1, ep_ai=5)
    c_bohr.set_data([new_arch.bohr_parameter()], [new_arch.fold_change()])

    # Update the shift indicators
    c_shift_line.set_ydata([0, -np.log(new_pact/pact_ref)])
    c_bohr_line.set_xdata([0, new_arch.bohr_parameter()])


def anim_epshift(i):

    # Update the shift plot
    ep_shift = ep_range[i]
    ep_point.set_data([ep_shift - ep_ref], [ep_shift - ep_ref])

    # Update the bohr plot
    new_arch = phd.thermo.SimpleRepression(R_ref, ep_shift, effector_conc=c_shift,
                                           ka=200, ki=1, ep_ai=500)
    ep_bohr.set_data([new_arch.bohr_parameter()], [new_arch.fold_change()])

    # Update the shift indicators
    ep_shift_line.set_ydata([0, ep_shift - ep_ref])
    ep_bohr_line.set_xdata([0, new_arch.bohr_parameter()])

def all_anim(i):
    anim_rshift(i)
    anim_cshift(i)
    anim_epshift(i)

plt.tight_layout()

if (REP == True) & (EP==True) & (C == True):
    ani = FuncAnimation(fig, all_anim, interval=1, frames=np.arange(0, len(rep_range)))
    name = 'all_shift'
else:
    if REP:
        ani = FuncAnimation(fig, anim_rshift, interval=1, frames=np.arange(0, len(rep_range)))
        name = 'r_shift'
    if EP:
        ani = FuncAnimation(fig, anim_epshift, interval=1, frames=np.arange(0, len(rep_range)))
        name = 'ep_shift'
    if C:
        ani = FuncAnimation(fig, anim_cshift, interval=1, frames=np.arange(0, len(rep_range)))
        name = 'c_shift'

ani.save(f'../figs/slide13_{name}.mp4', fps=50)
# %%
