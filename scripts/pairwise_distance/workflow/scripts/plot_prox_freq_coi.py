import skimage.io
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np

prox_counts = skimage.io.imread(snakemake.input[0])
cooccurances= skimage.io.imread(snakemake.input[1])
fig, ax = plt.subplots(figsize=(20,20))
normed = prox_counts/cooccurances

mpl_im = ax.imshow(np.log10(normed+0.0000000001), cmap = 'Reds')

ax.tick_params(bottom=False, top=False, left=False, labelleft=False, labelbottom=False, labeltop=False)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)

cbar = plt.colorbar(mpl_im, cax=cax)

cbar.set_label("log10(Proximity Frequency)")

plt.savefig(snakemake.output[0])