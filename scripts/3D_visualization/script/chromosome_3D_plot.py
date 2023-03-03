import pandas as pd
import numpy as np
from pathlib import Path
from mayavi import mlab
import configparser

# read in cell information and needed data for plotting
config = configparser.ConfigParser()
config.read('chromosome_config.ini')

experiment = config["CellInfo"]["experiment"]
rep = int(config["CellInfo"]["rep"])
pos = int(config["CellInfo"]["pos"])
cell_id = int(config["CellInfo"]["cell_id"])
chroms = config["CellInfo"]["chroms"].split(",")

dot_file = config["DataStorage"]["dot_file"]
input_path = Path(config["DataStorage"]["input_path"])

dbscan = "dbscan_ldp_nbr_allele"
fname = str(input_path / f"{experiment}_rep_{rep}_pos_{pos}_cell_{cell_id}_leiden_0.csv")

df = pd.read_csv(fname)
df = df[~df["name"].str.contains("control")]
# sort df according to allele
df["g"] = df["name"].str.split("-").str[1].astype(int)
df = df.sort_values(by = ["chrom", dbscan, "g"])

figs = []
for chrom_show in chroms:
    fig = mlab.figure(f"{chrom_show}", size = (2048,2048),  bgcolor=(1, 1, 1), fgcolor = (0.5, 0.5, 0.5))
    df_show = df[df["chrom"] == chrom_show]

    allele_id = df_show[dbscan].values
    alleles = df_show[dbscan].unique()

    for al in alleles:
        if al != -1:
            idx = allele_id == al

            al_df = df_show[df_show[dbscan]== al]
            al_df = al_df.sort_values(by = "g")

            x = al_df["x_um"].values
            y = -al_df["y_um"].values
            z = al_df["z_um"].values

            s = []
            connections = []
            for i in range(len(x) - 1):
                s.append(i)
                connections.append(np.array([i, i+ 1]))
            s.append(i + 1)
            connections = np.array(connections).astype(float)

            pts = mlab.points3d(x, y, z, s, scale_mode='vector', scale_factor = 0.35, colormap = "viridis", opacity = 1, mode='sphere', resolution = 60)
            pts.mlab_source.dataset.lines = np.array(connections)
            pts.mlab_source.on_ratio = 60
            tube = mlab.pipeline.tube(pts, tube_radius=0.1)
            tube.filter.radius_factor = 1.
            tube.filter.vary_radius = 'vary_radius_by_scalar'
            mlab.pipeline.surface(tube, colormap = "viridis", opacity = 1)
    mlab.view(0, 0)
    figs.append(fig)

# sync camera
if len(figs) > 1:
    for sfig in figs[1:]:
        mlab.sync_camera(figs[0], sfig)

# set light of scene
for fig in figs:
    fig.scene.light_manager.light_mode = "vtk"
    camera_light0 = fig.scene.light_manager.lights[0]
    camera_light0.intensity = 0.8
    camera_light1 = fig.scene.light_manager.lights[1]
    camera_light1.intensity = 0.4
    camera_light0.activate = True
    camera_light1.activate = True
