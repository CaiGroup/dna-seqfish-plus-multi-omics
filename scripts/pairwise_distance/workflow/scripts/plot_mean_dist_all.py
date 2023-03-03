import numpy as np
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import skimage.io
import PIL

PIL.Image.MAX_IMAGE_PIXELS = 272449424


img_nbrs =  skimage.io.imread(snakemake.input[0])

missing_pairs = img_nbrs == 2**16-1
img_nbrs[missing_pairs] = 0
img_nbrs_microns = img_nbrs/1000

fig, ax = plt.subplots(figsize=(15,15))

q_im = ax.imshow(img_nbrs_microns, cmap="Reds_r")
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = plt.colorbar(q_im, cax=cax)

img_nbrs_normed = img_nbrs/np.max(img_nbrs)

reds_r = mpl.colormaps['Reds_r']
cimg = reds_r(img_nbrs_normed)

cimg[missing_pairs] = np.array([0.5, 0.5, 0.5, 1])

q_im = ax.imshow(cimg)
ax.tick_params(bottom=False, top=False, left=False, labelleft=False, labelbottom=False, labeltop=False)


cbar.set_label("Mean Distance (microns)")

plt.savefig(snakemake.output[0])
