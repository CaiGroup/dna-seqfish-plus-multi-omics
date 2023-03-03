import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np

mean_chrm_pw_dists = np.loadtxt(snakemake.input[0], delimiter=",", dtype=float)
mean_chrm_pw_dists /= 1000
fig, ax = plt.subplots(figsize=(20,20))

mpl_im = ax.imshow(mean_chrm_pw_dists, cmap = 'Reds')

ax.tick_params(bottom=False, top=False, left=False, labelleft=False, labelbottom=False, labeltop=False)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)

cbar = plt.colorbar(mpl_im, cax=cax)

cbar.set_label("Distance (microns)")

plt.savefig(snakemake.output[0])