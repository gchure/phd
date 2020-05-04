#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import phd.viz
import scipy.io
import phd.image
import scipy.ndimage
import scipy.signal
import skimage.morphology
import skimage.filters
import skimage.measure
import skimage.io
import glob
colors, color_list = phd.viz.phd_style()
IP_DIST = 0.065 # in nm per pixel

#%%
# Find the root folders for each snapshot
folders = glob.glob('../../data/ch8_growth_si/example_mats/*')
fig, ax = plt.subplots(1, 2, figsize=(5, 2.5), dpi=100)
phd.viz.despine(ax.ravel())
for a in ax:
    a.set_ylabel('length [µm]')
    a.xaxis.grid(False)

locs = {32:0, 37: 20, 42: 40, 'acetate': 0, 'glycerol': 20, 'glucose': 40}
temp_colors = {32: colors['blue'], 37:colors['purple'], 42: colors['red']}
carb_colors = {'acetate':colors['dark_brown'], 'glucose':colors['purple'], 
               'glycerol':colors['green']}
lengths = {'glucose':[], 'glycerol':[], 'acetate':[], 37:[], 42:[], 32:[]}
widths = {'glucose':[], 'glycerol':[], 'acetate':[], 37:[], 42:[], 32:[]}
for f in folders:
    # Determine the temperature
    carbon, temp, _ = f.split('/')[-1].split('_')
    temp = int(temp[:-1]) 
    # Determine the axis
    if carbon != 'glucose':
        _ax = [ax[0]]
        color = carb_colors[carbon]
        loc = [locs[carbon]]
    if temp != 37:
        _ax = [ax[1]]
        color = temp_colors[temp]
        loc = [locs[temp]]
    if (carbon == 'glucose') & (temp==37):
        _ax = ax
        color = colors['purple']
        loc = [locs[carbon], locs[temp]]

    # Get the mat files
    mats = glob.glob(f'{f}/*.mat')

    # Instantiate a data frame for easy calculation of the mean. 
    dfs = [] 
    _lengths = []
    _widths = []
    for j, m in enumerate(mats):
        mat = scipy.io.loadmat(m)
        cellA = mat['CellA'][0][0][0][0]
        seg = cellA[3]
        length, width = cellA[9][0]
        _lengths.append(length)
        _widths.append(width)
        lab = skimage.measure.label(seg)
        props = skimage.measure.regionprops(lab)
        prop = [p for p in props]
        o = props[0].orientation
        rot = scipy.ndimage.rotate(seg, -o * 180/np.pi)
        conts = skimage.measure.find_contours(rot, 0)[0]
        _conts_x = conts[:,0] - np.min(conts[:, 0])
        _conts_y = conts[:, 1] - np.min(conts[:, 1])

        if (np.sum(seg) < 1000) & (np.sum(seg) > 50):
            # Circularly permute the contour from the minimum x
            ind = np.where(_conts_y == 0)[0][0]
            perm_conts_y = np.roll(_conts_y, -ind)
            filt_conts_y = scipy.signal.savgol_filter(perm_conts_y, polyorder=2, window_length=15)
            perm_conts_x = np.roll(_conts_x, -ind)
            filt_conts_x = scipy.signal.savgol_filter(perm_conts_x, polyorder=2, window_length=15)
            _df = pd.DataFrame().from_dict({'x':perm_conts_x, 'idx':np.arange(len(conts)),
                                          'y':perm_conts_y * IP_DIST})
            _df['cell_id'] = j
            dfs.append(_df)
            for i, a in enumerate(_ax):
                a.plot(filt_conts_y + loc[i], filt_conts_x * IP_DIST, lw=0.1, 
                        alpha=0.75, color=color)
    # Compute the average lengths and widths
    avg_length = np.mean(np.array(_lengths))
    avg_width = np.mean(np.array(_widths)) 
    for k, a in enumerate(_ax):
        a.vlines(loc[k], IP_DIST * avg_width/2,  (avg_length -  avg_width) * IP_DIST, 'k', lw=2, zorder=1000)
        a.vlines(loc[k] + avg_width, IP_DIST * avg_width/2,  (avg_length - avg_width) * IP_DIST, 'k', lw=2, zorder=1000)
        bottom_arc = matplotlib.patches.Arc((loc[k] + avg_width/2, IP_DIST * avg_width/2), 
                      width=avg_width, height=IP_DIST * avg_width, theta1=180, 
                      theta2=0, color='k', zorder=1000, lw=2) 
        top_arc = matplotlib.patches.Arc((loc[k] + avg_width/2, IP_DIST * (avg_length - avg_width)), 
                      width=avg_width, height=IP_DIST * avg_width, theta1=0, 
                      theta2=180, color='k', zorder=1000, lw=2) 
        a.add_patch(bottom_arc)
        a.add_patch(top_arc)


# Compute the average shape parameters

ax[0].set_xlim([-5, 65])    
ax[0].set_xticks([10, 30, 50])
ax[0].set_xticklabels(['acetate', 'glycerol', 'glucose'])
ax[1].set_xlim([-5, 65])    
ax[1].set_xticks([10, 30, 50])
ax[1].set_xticklabels(['32° C', '37° C', '42° C'])
plt.tight_layout()
fig.text(0, 0.9, '(A)', fontsize=8)
fig.text(0.5, 0.9, '(B)', fontsize=8)
plt.savefig('../figs/ch8_figS2.pdf', bbox_inches='tight')
plt.savefig('../figs/ch8_figS2.png', bbox_inches='tight')

# %%
