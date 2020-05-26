
#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import phd.viz
import phd.thermo
import seaborn as sns
colors, palette = phd.viz.phd_style()

# Define the modulators
KAKI = True
EPAI = True
R = True
EPR = True

# Define the reference states. 
R_ref = 200
epr_ref = -14
epai_ref = 1
Ka_ref = 200
Ki_ref = 1


# DEfine the perturbed states and properly roll for animation purposes
def rollarray(x, loc, palette=False):
    _rolled = np.roll(x, loc)
    largev1 = _rolled[:loc]
    largev2 = _rolled[loc:len(x) - loc]
    tot = np.array([largev1, largev1[::-1], largev2[::-1], largev2])
    return np.concatenate(tot).ravel()

def genpalette():
    _palette =  sns.color_palette('GnBu_d', n_colors=200)
    mid_bright = _palette[100:]
    bright_mid = mid_bright[::-1]
    mid_dark = _palette[:101]
    dark_mid = mid_dark[::-1]
    return np.concatenate([mid_bright, bright_mid, mid_dark, dark_mid])

_palette = genpalette()

Ka_mut = np.logspace(0, 4, 200)
loc = np.where(np.isclose(Ka_mut, Ka_ref, atol=5))[0][0]
Ka_mut = rollarray(Ka_mut, -loc)
Ki_mut = Ka_mut * Ki_ref / Ka_ref


# EpAI
epai_mut = np.linspace(-9, 9, 200)
loc = np.where(np.isclose(epai_mut, epai_ref, atol=0.05))[0][0]
epai_mut = rollarray(epai_mut, -loc)


# R 
r_mut = np.logspace(0, 4, 200)
loc = np.where(np.isclose(r_mut, R_ref, atol=5))[0][0]
r_mut = rollarray(r_mut, -loc)

# epR
epr_mut = np.linspace(-18, -10, 200)
loc = np.where(np.isclose(epr_mut, epr_ref, atol=0.1))[0][0]
epr_mut = rollarray(epr_mut, -loc)


# Define the reference architecture
c_range = np.logspace(-4, 4, 200)
ref_arch = phd.thermo.SimpleRepression(R_ref, epr_ref, ka=Ka_ref, ki=Ki_ref,
                                    ep_ai=epai_ref, effector_conc=c_range).fold_change()
ref_mwc = phd.thermo.MWC(ka=Ka_ref, ki=Ki_ref, ep_ai=epai_ref, effector_conc=c_range).pact()

#%%
fig, ax = plt.subplots(2, 4, figsize=(8,4), sharex=True)
for i in range(4):
    # ax[0, i].axis('off')
    ax[1, i].set_xlabel('$c / K_A^\mathrm{(wt)}$')
    ax[0, i].set_ylabel('fold-change')
    ax[1, i].set_ylabel('∆F [$k_BT$]')

for a in ax.ravel():
    phd.viz.despine(a)
    a.set_xscale('log')

# Set the scales
for i in range(4):
    ax[1, i].set_ylim([-15, 15])

# Add titleboxes
# phd.viz.titlebox(ax[0, 0], 'change in $K_A$ and $K_I$', color=colors['black'],
#                  boxsize='12%')
# phd.viz.titlebox(ax[0, 1], r'change in $\Delta\varepsilon_{AI}$', color=colors['black'],
#                   boxsize='12%')
# phd.viz.titlebox(ax[0, 2], r'change in $R$', color=colors['black'],
#                  boxsize='12%')
# phd.viz.titlebox(ax[0, 3], r'change in $\Delta\varepsilon_{RA}$', color=colors['black'],
#                  boxsize='12%')


# Add the wild-type features
for i in range(4):
    ax[1, i].hlines(0, 1E-4/Ka_ref, 1E4/Ka_ref, color=colors['black'], lw=0.75, zorder=1000) 
    ax[0, i].plot(c_range/Ka_ref, ref_arch, '-', color=colors['black'], label='__nolegend__',
                lw=0.75, zorder=1000)


# Instantiate the mutant architectures
kaki_arch = phd.thermo.SimpleRepression(R_ref, epr_ref, ka=Ka_mut[0], ki=Ki_mut[0],
                                        ep_ai=epai_ref, effector_conc=c_range).fold_change()
kaki_mwc = phd.thermo.MWC(ka=Ka_mut[0], ki=Ki_mut[0], ep_ai=epai_ref, effector_conc=c_range).pact()

epai_arch = phd.thermo.SimpleRepression(R_ref, epr_ref, ka=Ka_ref, ki=Ki_ref,
                                        ep_ai=epai_mut[0], effector_conc=c_range).fold_change()
epai_mwc = phd.thermo.MWC(ka=Ka_ref, ki=Ki_ref, ep_ai=epai_mut[0], effector_conc=c_range).pact()

rep_arch = phd.thermo.SimpleRepression(r_mut[0], epr_ref, ka=Ka_ref, ki=Ki_ref,
                                        ep_ai=epai_ref, effector_conc=c_range).fold_change()

epr_arch = phd.thermo.SimpleRepression(R_ref, epr_mut[0], ka=Ka_ref, ki=Ki_ref,
                                        ep_ai=epai_ref, effector_conc=c_range).fold_change()


# Instantiate the mutant curves
kaki_line, = ax[0, 0].plot(c_range/Ka_ref, kaki_arch, color=colors['blue'])
kaki_delf, = ax[1, 0].plot(c_range/Ka_ref, -np.log(kaki_mwc/ref_mwc), color=colors['blue'])


epai_line, = ax[0, 1].plot(c_range/Ka_ref, epai_arch, color=colors['blue'])
epai_delf, = ax[1, 1].plot(c_range/Ka_ref, -np.log(epai_mwc/ref_mwc), color=colors['blue'])



r_line, = ax[0, 2].plot(c_range/Ka_ref, rep_arch, color=colors['red'], label=r'$R^\mathrm{(mut)} = R^\mathrm{(wt)}$')
r_delf, = ax[1, 2].plot([c_range[0]/Ka_ref, c_range[-1]/Ka_ref], [-np.log(r_mut[0]/R_ref), -np.log(r_mut[0]/R_ref)],
                    color=colors['red'], label=r'$R^\mathrm{(mut)} = R^\mathrm{(wt)}$')


epr_line, = ax[0, 3].plot(c_range/Ka_ref, epr_arch, color=colors['orange'], label=r'$\Delta\varepsilon_{RA}^\mathrm{(mut)} > \Delta\varepsilon_{RA}^\mathrm{(wt)}$')
epr_delf, = ax[1, 3].plot([c_range[0]/Ka_ref, c_range[-1]/Ka_ref], [epr_mut[0] - epr_ref, epr_mut[0] - epr_ref],
                    color=colors['orange'], label=r'$\Delta\varepsilon_{RA}^\mathrm{(mut)} > \Delta\varepsilon_{RA}^\mathrm{(wt)}$')

def anim_kaki(i):
    kaki_arch = phd.thermo.SimpleRepression(R_ref, epr_ref, ka=Ka_mut[i],
                                            ki=Ki_mut[i], ep_ai=epai_ref,
                                            effector_conc=c_range).fold_change()
    kaki_mwc = phd.thermo.MWC(ka=Ka_mut[i], ki=Ki_mut[i], ep_ai=epai_ref,
                              effector_conc=c_range).pact()
    if Ka_mut[i] < Ka_ref:
        color = colors['blue']
    else:
        color = colors['green']
    kaki_line.set_color(color)
    kaki_line.set_ydata(kaki_arch)
    kaki_delf.set_ydata(-np.log(kaki_mwc/ref_mwc))
    kaki_delf.set_color(color)


def anim_epai(i):
    epai_arch = phd.thermo.SimpleRepression(R_ref, epr_ref, ka=Ka_ref,
                                            ki=Ki_ref, ep_ai=epai_mut[i],
                                            effector_conc=c_range).fold_change()
    epai_mwc = phd.thermo.MWC(ka=Ka_ref, ki=Ki_ref, ep_ai=epai_mut[i],
                              effector_conc=c_range).pact()
    if epai_mut[i] < epai_ref:
        color=colors['blue']
    else:
        color = colors['green']
    epai_line.set_ydata(epai_arch)
    epai_delf.set_ydata(-np.log(epai_mwc/ref_mwc))
    epai_line.set_color(color)
    epai_delf.set_color(color)

def anim_r(i):
    r_arch = phd.thermo.SimpleRepression(r_mut[i], epr_ref, ka=Ka_ref,
                                            ki=Ki_ref, ep_ai=epai_ref,
                                            effector_conc=c_range).fold_change()
    r_line.set_ydata(r_arch) 
    r_delf.set_ydata([-np.log(r_mut[i]/R_ref), -np.log(r_mut[i]/R_ref)])
    if r_mut[i] < R_ref:
        color = colors['blue']
    else:
        color = colors['green']
    r_line.set_color(color)
    r_delf.set_color(color)


def anim_epr(i):
    epr_arch = phd.thermo.SimpleRepression(R_ref, epr_mut[i], ka=Ka_ref,
                                            ki=Ki_ref, ep_ai=epai_ref,
                                            effector_conc=c_range).fold_change()
 
    epr_line.set_ydata(epr_arch)
    epr_delf.set_ydata([epr_mut[i] - epr_ref, epr_mut[i] - epr_ref])
    if epr_mut[i] < epr_ref:
        color = colors['blue']
    else:
        color = colors['green']
    epr_line.set_color(color)
    epr_delf.set_color(color)


plt.tight_layout()
def anim_all(i):
    if KAKI == True:
        anim_kaki(i)
    if EPAI == True:
        anim_epai(i)
    if R == True:
        anim_r(i)
    if EPR == True:
        anim_epr(i)

ani = FuncAnimation(fig, anim_all, frames=np.arange(0, len(Ka_mut)))

name = ''
if KAKI == True:
    name += '_kaki'
if EPAI == True:
    name += '_epai'
if R == True:
    name += '_R'
if EPR == True:
    name += '_epr'

ani.save(f'../figs/mutant_delF{name}.mp4', fps=100)
# %%


# %%

# Set up small gifs of just the delF for demonstrative purposes

#%% EpRA
fig, ax = plt.subplots(1, 1, figsize=(2, 2))
phd.viz.despine(ax)
ax.set_xscale('log')
ax.set_ylabel(r'∆F [$k_BT$]')
ax.set_xlabel(r'$c / K_A^\mathrm{(wt)}$')
phd.viz.titlebox(ax, 'DNA binding energy', color=colors['black'])
ax.hlines(0, c_range[0]/Ka_ref, c_range[-1]/Ka_ref, color=colors['black'], lw=0.75)
ax.set_ylim([-5, 5])
epr_delf, = ax.plot(c_range/Ka_ref, [epr_mut[0] - epr_ref] * np.ones_like(c_range), 
                    color=colors['blue'])

def anim(i):
    epr_delf.set_ydata([epr_mut[i] - epr_ref] * np.ones_like(c_range)) 
    if epr_mut[i] < epr_ref:
        color = colors['blue']
    else:
        color = colors['green']
    epr_delf.set_color(color)

ani = FuncAnimation(fig, anim, frames=np.arange(0, len(epr_mut)))
plt.tight_layout()
ani.save('../figs/epr_delf_example.mp4', fps=80)

# %%
fig, ax = plt.subplots(2, 1, figsize=(2, 4), sharex=True)
phd.viz.despine(ax)
for a in ax.ravel():
    a.set_xscale('log')
    a.set_ylabel(r'∆F [$k_BT$]')
    a.set_ylim([-9, 9])

ax[1].set_xlabel(r'$c / K_A^\mathrm{(wt)}$')
phd.viz.titlebox(ax[0], 'inducer binding constants', color=colors['black'])
phd.viz.titlebox(ax[1], 'allosteric energy difference', color=colors['black'])
ax[0].hlines(0, c_range[0]/Ka_ref, c_range[-1] / Ka_ref, color=colors['black'], lw=0.75)
ax[1].hlines(0, c_range[0]/Ka_ref, c_range[-1]/Ka_ref , color=colors['black'], lw=0.75)
kaki_delf, = ax[0].plot(c_range/Ka_ref, -np.log(kaki_mwc/ ref_arch),
                    color=colors['blue'])
epai_delf, = ax[1].plot(c_range/Ka_ref, -np.log(epai_mwc/ ref_arch),
                    color=colors['blue'])


def anim(i):
    _kaki_mwc = phd.thermo.MWC(ka=Ka_mut[i], ki=Ki_mut[i], ep_ai=epai_ref, 
                                effector_conc=c_range).pact()
    _epai_mwc = phd.thermo.MWC(ka=Ka_ref, ki=Ki_ref, ep_ai=epai_mut[i],
                            effector_conc=c_range).pact()
    if Ka_mut[i] < Ka_ref:
        kaki_color = colors['blue']
    else:
        kaki_color = colors['green']
    if epai_mut[i] < epai_ref:
        epai_color = colors['blue']
    else:
        epai_color=colors['green']
    kaki_delf.set_ydata(-np.log(_kaki_mwc / ref_mwc))
    epai_delf.set_ydata(-np.log(_epai_mwc / ref_mwc))
    kaki_delf.set_color(kaki_color)
    epai_delf.set_color(epai_color)

ani = FuncAnimation(fig, anim, frames=np.arange(0, len(epai_mut)))
plt.tight_layout()
ani.save('../figs/allo_delf_example.mp4', fps=80)



# %%
