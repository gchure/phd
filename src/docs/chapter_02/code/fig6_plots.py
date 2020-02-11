#%%
import numpy as np
import matplotlib.pyplot as plt
import  matplotlib.gridspec as gridspec 
import pandas as pd
import phd.viz
import phd.stats
import pickle 
colors, palette = phd.viz.phd_style()
constants = phd.thermo.load_constants()

# Load the data set
data = pd.read_csv('../../data/ch2_induction/RazoMejia_2018.csv', comment='#')

# Load the flatchains for the prediction measurements. 
with open('../../data/ch2_induction/mcmc/SI_I_O2_R260.pkl', 'rb') as file:
    unpickler = pickle.Unpickler(file)
    gauss_flatchain = unpickler.load()
    gauss_flatlnprobability = unpickler.load()
ka_fc = np.exp(-gauss_flatchain[:, 0])[::100]
ki_fc = np.exp(-gauss_flatchain[:, 1])[::100]


#%%
# Compute the theoretical property curves. 
rep_range = np.logspace(0, 4, 200)
prop_df = pd.DataFrame([])
ops = {'O1':constants['O1'], 'O2':constants['O2'], 'O3':constants['O3']}
for op, ep in ops.items():
    for i, r in enumerate(rep_range):
        arch = phd.thermo.SimpleRepression(r, ep, ka=ka_fc, ki=ki_fc, 
                                          ep_ai=constants['ep_AI'], effector_conc=0)
        props = arch.compute_properties()
        for prop, val in props.items():
            if prop == 'leakiness':
                hpd_min, hpd_max = val, val
            else:
                hpd_min, hpd_max = phd.stats.compute_hpd(val, 0.95)
            prop_df = prop_df.append({'operator':op, 'binding_energy':ep,
                                      prop:val, 'repressors':r,
                                      'hpd_min':hpd_min, 'hpd_max':hpd_max},
                                      ignore_index=True)

# %%

#Instantiate the figure canvas. 
fig, ax = plt.subplots(2, 3, figsize=(6, 3.5), dpi=100)

# Format axes as needed. 
for i, a in enumerate(ax.ravel()):
    a.set_xscale('log')

ax[0, 0].set_yscale('log')
ax[0, 0].set_ylabel('leakiness')
ax[0, 1].set_ylabel('saturation')
ax[0, 2].set_ylabel('dynamic range')
ax[1, 0].set_ylabel(r'$EC_{50}$ [µM]')
ax[1, 0].set_yscale('log')
ax[1, 1].set_ylabel('effective Hill coefficient')
ax[1, 2].axis('off')

# Define the color palette
viridis = 


plt.tight_layout()

# %%
