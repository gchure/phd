#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import phd.viz
import phd.thermo
import seaborn as sns
colors, palette = phd.viz.phd_style()

magma = sns.color_palette('viridis', n_colors=6)

# Define some curves
rep_range = np.logspace(0, 3, 5)

# Define the marker styles
markers = ['o', 'v', 'X', '^', 'd']
glyphs = {r:m for r, m in zip(rep_range, markers)}

ep_r = -12
c_range = np.logspace(-2, 3.5, 200)

r, e, c  = np.meshgrid(rep_range, ep_r, c_range)
arch = phd.thermo.SimpleRepression(r, e, effector_conc=c, ka=200, ki=1, ep_ai=5)
fc_curves = arch.fold_change()
bohr_curves = arch.bohr_parameter()

# Define the theory curves.
bohr_range = np.linspace(-10, 10, 200)
collapse = (1 + np.exp(-bohr_range))**-1

# %%

# Instantiate and format axis
fig, ax = plt.subplots(2, 1, figsize=(3, 4))
for a in ax:
    phd.viz.despine(a)
ax[0].set_xscale('log')
ax[0].set_xlabel('inducer concentration')
ax[1].set_xlabel('free energy [$k_BT$]')
ax[0].set_ylabel('fold-change')
ax[1].set_ylabel('fold-change')
ax[1].set_xlim([-10, 10])
phd.viz.titlebox(ax[0], 'PARAMETRIC DETAILS', fontsize=6, color=colors['black'],
                bgcolor=colors['grey'], pad=0.05, boxsize='10%')
phd.viz.titlebox(ax[1], 'PHENOTYPIC OUTPUT', fontsize=6, color=colors['black'],
                bgcolor=colors['grey'], pad=0.05, boxsize='10%')

ax[1].plot(bohr_range, collapse, 'k-')
fc_points = []
bohr_points = []
for i, r in enumerate(rep_range):
    ax[0].plot(c_range, fc_curves[0, i,  :], '-', color = magma[i])
    ax[0].plot(c_range[::10], fc_curves[0, i, ::10], linestyle='none', marker=glyphs[r],
        ms=5, color=magma[i], markeredgecolor=colors['grey'], markeredgewidth=0.5)
    ax[1].plot(bohr_curves[0, i,  ::10], fc_curves[0, i, ::10], color='lightgrey',
              ms=5, linestyle='none', marker=glyphs[r], zorder=-1)
    fc_line, = ax[0].plot([], [], linestyle='none', marker=glyphs[r],
                       ms=5, markeredgecolor=magma[i], markerfacecolor=colors['grey'], 
                       markeredgewidth=0.5)
    bohr_line, = ax[1].plot([], [], linestyle='none', marker=glyphs[r],
                       ms=5, markerfacecolor=magma[i], markeredgecolor=colors['grey'], 
                       markeredgewidth=0.5)

    fc_points.append(fc_line)
    bohr_points.append(bohr_line)

ax[0].set_xlim([1E-2, 1.5E3])
hline1,  = ax[0].plot([1E-4, 3E4], [1.0, 1.0], color=colors['black'], lw=8, alpha=0.25)
hline2,  = ax[1].plot([-15, 15], [1.0, 1.0], color=colors['black'], lw=8, alpha=0.25)
window = 0.01
N = 100
selections = np.linspace(0, 1, N)
fudge = 0.02

def animate(i):
    hline1.set_ydata([selections[i], selections[i]])
    hline2.set_ydata([selections[i], selections[i]])
    for j, p in enumerate(fc_points):
            _fc = fc_curves[0, j, ::10]
            _idx = (_fc >= selections[i] - fudge) & (_fc <= selections[i] + fudge) 
            _bohr = bohr_curves[0, j, ::10]
            # print(_idx)
            _fc = _fc[_idx]
            _bohr = _bohr[_idx]

            if len(_fc) > 0:
                p.set_data(c_range[::10][_idx],  _fc)
                bohr_points[j].set_data(_bohr, _fc)
            else:
                p.set_data([], [])
                bohr_points[j].set_data([], [])

ani = FuncAnimation(fig, animate, interval=1, frames=np.arange(0, N))
# ani.save('./pedagogical_collapse.mp4', fps=10)
plt.tight_layout()
ani.save('../figs/pedagogical_collapse.gif', writer='imagemagick', fps=10)

    # %%
