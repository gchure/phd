#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import  matplotlib.gridspec as gridspec
import phd.viz
colors, palette = phd.viz.phd_style()

# Load the data.
data = pd.read_csv('../../data/ch9_mscl_si/MLG910_electrophysiology.csv')
data.columns = ['time', 'pa', 'mmHg']

# Instantiate the figure.
fig = plt.figure(figsize=(4, 4))
gs = gridspec.GridSpec(3, 1)
ax1 = fig.add_subplot(gs[:2, 0])
ax2 = fig.add_subplot(gs[2, 0])

# Format the axes
phd.viz.despine([ax1, ax2])
ax1.set_xlim([19, 23.8])
ax2.set_xlim([19, 23.8])
ax1.set_ylim([0, 525])
ax2.set_ylim([-350, 0])

# Format the axes
ax1.xaxis.set_ticklabels([])

# Add the appropriate labels
ax1.set_ylabel('current [pA]')
ax2.set_ylabel('pressure\n [mmHg]')
ax2.set_xlabel('time [s]')

# Add marker labels
ax1.text(0.08, 0.93, 'MscS', fontsize=8, transform=ax1.transAxes)
        
ax2.text(0.46, 0.93, 'MscL-sfGFP', fontsize=8, transform=ax1.transAxes)

# Plot the traces and color red
_ = ax1.plot(data['time'], data['pa'], '-', color=colors['purple'], lw=0.5)
_ = ax2.plot(data['time'], data['mmHg'], '-', color=colors['purple'], lw=0.75)

# Label the MscS points
_ = ax1.vlines(19.6, -1.5, 525, lw=31, color=colors['light_blue'], zorder=-1, alpha=0.5)
_ = ax2.vlines(19.6, 4, -350, lw=31, color=colors['light_blue'], zorder=-1, alpha=0.5)

# Label the MscL points
_ = ax1.vlines(21.7, -1.5, 525, lw=100, color=colors['light_orange'], zorder=-1, alpha=0.5)
_ = ax2.vlines(21.7, 4, -350, lw=100, color=colors['light_orange'], zorder=-1, alpha=0.5)
plt.savefig('../figs/ch9_figS1.pdf', bbox_inches='tight')
plt.savefig('../figs/ch9_figS1.png', bbox_inches='tight') 
  

# %%
