#%%
import numpy as np
import pandas as pd
import phd.stats
import pickle
import glob

# %%
MASS_FRAC = 0.95
files = glob.glob('../../data/ch2_induction/mcmc/SI*.pkl')
df = pd.DataFrame([])
for f in files:
    _, _, op, rep = f.split('/')[-1].split('_')
    rep = int(rep.split('.')[0][1:])
    with open(f, 'rb') as _f:
        unpickler = pickle.Unpickler(_f)
        gauss_flatchain = unpickler.load()
        gauss_flatlnprob = unpickler.load()
    max_idx = np.argmax(gauss_flatlnprob, axis=0)
    ea, ei, sigma = gauss_flatchain[max_idx]
    ka_mode = np.exp(-ea)
    ki_mode = np.exp(-ei)
    ka = np.exp(-gauss_flatchain[:, 0])
    ki = np.exp(-gauss_flatchain[:, 1])
    ka_min, ka_max = phd.stats.compute_hpd(ka, MASS_FRAC)
    ki_min, ki_max = phd.stats.compute_hpd(ki, MASS_FRAC)
    df = df.append({'parameter': 'ka',
                    'hpd_min': ka_min,
                    'hpd_max': ka_max,
                    'mode':ka_mode,
                    'mass_frac':MASS_FRAC,
                    'operator': op,
                    'repressors': rep},
                    ignore_index=True)
    df = df.append({'parameter': 'ki',
                    'hpd_min': ki_min,
                    'hpd_max': ki_max,
                    'mode':ki_mode,
                    'mass_frac':MASS_FRAC,
                    'operator': op,
                    'repressors': rep},
                    ignore_index=True)
df.to_csv('../../data/ch2_induction/RazoMejia_KaKi_estimates.csv', index=False)


# %%
