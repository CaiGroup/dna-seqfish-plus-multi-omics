import pandas as pd
from pathlib import Path
import sys
sys.path.insert(0, "./")
import numpy as np
from scipy.stats import zscore
from util import pil_imread
from mayavi import mlab
from matplotlib import colors
from tvtk.api import tvtk
import configparser

# read in cell information and needed data for plotting
config = configparser.ConfigParser()
config.read('H4K20me3_MajSat_volume_config.ini')
# cell information
rep = int(config["CellInfo"]["rep"])
pos = int(config["CellInfo"]["pos"])
cell_id = int(config["CellInfo"]["cell_id"])
# data storage pathway
dot_file = config["DataStorage"]["dot_file"]
IF_file = config["DataStorage"]["IF_file"]
mask_file = config["DataStorage"]["mask_file"]
annot_file = config["DataStorage"]["annot_file"]
input_path = Path(config["DataStorage"]["input_path"])

# preset color for different chromosomes
color_name = ['#9370db', '#0000c8', '#3c64e6', '#789bf2', '#b0e0e6', '#20b2aa', '#9acd32', '#2e8b57', '#f5e6be', '#deb887',
              '#ffe100', '#ffa500', '#ff4500', '#b22222', '#ffb6c1', '#ff1493', '#e21f26', '#723b7a', '#295f8a', '#777777']
# transfer the Hex code to RGB code since Mayavi only take RGB color
rgba_name = [ np.array(colors.to_rgb(c)) for c in color_name]

# preset color for different gene family
color_v = (198 / 255, 148 / 255, 50 / 255) # Vmn
color_o = (129 / 255, 151 / 255, 244 / 255)# Olfr

# get the whole field of view IF voxel intensity, the data formated is MajSat intensity in first channel, and H4K20me3 intensity in second channel
IF_arr = pil_imread(str(input_path / IF_file),swapaxes=True)
ma_arr = pil_imread(str(input_path / mask_file),swapaxes=True)[0]

# read bin annotation and bin location file
bin_annot = pd.read_csv(str(input_path / annot_file))
df = pd.read_csv(str(input_path / dot_file))[["name", "cellID", "rep", "x", "y", "z", "x_um", "y_um", "z_um", "chrom"]]

# format data

# crop out 3D cube of the selected cell
z_span, y_span, x_span = np.where(ma_arr == cell_id)
pad = 2 # set a pad to avoid half dots
z_min, z_max = z_span.min() , z_span.max()
y_min, y_max = y_span.min() - pad, y_span.max() + pad
x_min, x_max = x_span.min() - pad, x_span.max() + pad
zoom = np.s_[z_span.min() :z_span.max(), y_span.min() - pad :y_span.max() + pad , x_span.min() - pad: x_span.max() + pad]

# crop out the IF_arr for a quick test
sub_arr = np.copy(IF_arr[0][zoom]).astype(float)
sub_mask = np.copy(ma_arr[zoom])
sub_arr[sub_mask != cell_id] = np.NaN
zscore_arr = zscore(sub_arr, axis=None, ddof=0, nan_policy='omit')
zscore_arr = np.nan_to_num(zscore_arr)

sub_arr1 = np.copy(IF_arr[1][zoom]).astype(float)
sub_mask = np.copy(ma_arr[zoom])
sub_arr1[sub_mask != cell_id] = np.NaN
zscore_arr1 = zscore(sub_arr1, axis=None, ddof=0, nan_policy='omit')
zscore_arr1 = np.nan_to_num(zscore_arr1)

# prepare dot location information
df = df[~df["chrom"].str.contains("control")]
df["chrom_id"] = df["chrom"].str[3:].replace("X", 20).astype(int)

# add dot annotation
df = df.merge(bin_annot)

df_show = df[df["H4K20me3_cat5"] == 2] # H4K20me3 strong bins
df_o = df[(df["Olfr"] == 1) & (df["H4K20me3_cat5"] != 0) ] # olfr family
df_v = df[(df["Vmn"] == 1)  & (df["H4K20me3_cat5"] != 0) ] # vmn family

# H4K20me3 dots, color by chromosome
x_show = df_show["x"].values
y_show = df_show["y"].values
z_show = df_show["z"].values
c_show = df_show["chrom_id"].values

# all detected dots, show as background dots
x = df["x"].values
y = df["y"].values
z = df["z"].values
c = df["chrom_id"].values

# olfr gene family dots
x_o = df_o["x"].values
y_o = df_o["y"].values
z_o = df_o["z"].values

# Vmn gene family dots
x_v = df_v["x"].values
y_v = df_v["y"].values
z_v = df_v["z"].values

# color for H4K20me3 bin chromosome
rgba_show = np.array([rgba_name[i - 1] for i in c_show]) * 255

# Plotting starts here
# final plot
fig1 = mlab.figure(f"MajSat", size = (2048,2048),  bgcolor=(1, 1, 1), fgcolor = (0.5, 0.5, 0.5))

src = mlab.pipeline.scalar_field(zscore_arr)
surface = mlab.pipeline.iso_surface(src, contours=[2], opacity=0.35, color = (124/255, 208/255, 103/255))
surface.actor.actor.scale = (1.92, 1.0, 1.0)

# plot all bg dots
p3d2 = mlab.points3d( (z-z_min) * 1.92, y-y_min, x-x_min,
                      scale_mode='vector', scale_factor = 0.7,
                      opacity = 0.13, mode='sphere', resolution = 30)


#change the color of pts
pts=mlab.points3d((z_show - z_min) * 1.92, y_show - y_min, x_show - x_min, resolution = 60)
sc=tvtk.UnsignedCharArray()
sc.from_array(rgba_show)
pts.mlab_source.dataset.point_data.scalars=sc
pts.mlab_source.dataset.modified()
pts.glyph.scale_mode = 'data_scaling_off'
pts.glyph.glyph.scale_factor = 2.5


fig2 = mlab.figure(f"H4K20me3", size = (2048,2048),  bgcolor=(1, 1, 1), fgcolor = (0.5, 0.5, 0.5))

src1 = mlab.pipeline.scalar_field(zscore_arr1)
surface1 = mlab.pipeline.iso_surface(src1, contours=[2], opacity=0.35, color = (207/255, 131/255, 194/255))
surface1.actor.actor.scale = (1.92, 1.0, 1.0)

# plot all bg dots
p3d2 = mlab.points3d( (z-z_min) * 1.92, y-y_min, x-x_min,
                      scale_mode='vector', scale_factor = 0.7,
                      opacity = 0.13, mode='sphere', resolution = 30)


#change the color of pts
pts=mlab.points3d((z_show - z_min) * 1.92, y_show - y_min, x_show - x_min, resolution = 60)
sc=tvtk.UnsignedCharArray()
sc.from_array(rgba_show)
pts.mlab_source.dataset.point_data.scalars=sc
pts.mlab_source.dataset.modified()
pts.glyph.scale_mode = 'data_scaling_off'
pts.glyph.glyph.scale_factor = 2.5

fig3 = mlab.figure(f"interchrom", size = (2048,2048),  bgcolor=(1, 1, 1), fgcolor = (0.5, 0.5, 0.5))

src1 = mlab.pipeline.scalar_field(zscore_arr1)
surface1 = mlab.pipeline.iso_surface(src1, contours=[2], opacity=0.35, color = (207/255, 131/255, 194/255))
surface1.actor.actor.scale = (1.92, 1.0, 1.0)

#change the color of pts
pts=mlab.points3d((z_show - z_min) * 1.92, y_show - y_min, x_show - x_min, resolution = 60)
sc=tvtk.UnsignedCharArray()
sc.from_array(rgba_show)
pts.mlab_source.dataset.point_data.scalars=sc
pts.mlab_source.dataset.modified()
pts.glyph.scale_mode = 'data_scaling_off'
pts.glyph.glyph.scale_factor = 2.5


fig4 = mlab.figure(f"gene_family", size = (2048,2048),  bgcolor=(1, 1, 1), fgcolor = (0.5, 0.5, 0.5))

src1 = mlab.pipeline.scalar_field(zscore_arr1)
surface1 = mlab.pipeline.iso_surface(src1, contours=[2], opacity=0.35, color = (207/255, 131/255, 194/255))
surface1.actor.actor.scale = (1.92, 1.0, 1.0)

# test 1
# plot all bg dots
p3d2 = mlab.points3d( (z-z_min) * 1.92, y-y_min, x-x_min,
                      scale_mode='vector', scale_factor = 0.7,
                      opacity = 0.13, mode='sphere', resolution = 30)
pts_o = mlab.points3d((z_o - z_min) * 1.92, y_o - y_min, x_o - x_min, resolution = 60, opacity = 1, scale_factor = 2.5, color = color_o)
pts_v = mlab.points3d((z_v - z_min) * 1.92, y_v - y_min, x_v - x_min, resolution = 60, opacity = 1, scale_factor = 2.5, color = color_v)

# mlab.show()
mlab.sync_camera(fig2, fig1)
mlab.sync_camera(fig2, fig3)
mlab.sync_camera(fig2, fig4)

# set lighting
for fig in [fig1, fig2, fig3, fig4]:
    fig.scene.light_manager.light_mode = "vtk"
    camera_light0 = fig.scene.light_manager.lights[0]
    camera_light0.intensity = 1
    camera_light1 = fig.scene.light_manager.lights[1]
    camera_light1.intensity = 0.3
    camera_light0.activate = True
    camera_light1.activate = True