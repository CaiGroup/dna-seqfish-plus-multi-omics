import os
from pathlib import Path

import pandas as pd
import numpy as np
import dask.dataframe as dd
from utils import *


# accept inputs from the user
def read_input():
    import argparse

    # get the input variables from the command line using argparse
    parser = argparse.ArgumentParser(description='Calculate the distance for each bin pair')
    parser.add_argument('--chrom', type=str, help='the chromosome number you want to check')
    parser.add_argument('--input_dir', type=str, help='the directory of the input files')
    parser.add_argument('--celltype', type=str, help='the directory of the input files')
    # parser.add_argument('--block_size', type=int, help='chunk size when breaking down the dataframe')
    parser.add_argument('--output_dir', type=str, help='the directory of the output files as median final result format')
    parser.add_argument('--rep', type=str, help='which replicate to choose')
    # parser.add_argument('--max_bin', type=str, help='the maximum bin number')

    # parse the arguments
    args = parser.parse_args()
    
    input_dir = args.input_dir
    chrom = args.chrom
    # block_size = args.block_size
    output_dir = args.output_dir
    rep = args.rep
    celltype = args.celltype
    # max_bin = args.max_bin

    return input_dir, chrom, output_dir, rep, celltype

# input_dir, output_dir = Path("/groups/CaiLab/personal/yujing/Yodai/100k/data/raw/E14/profile_per_dot/radical_organization/"), Path("/home/yy4ng/group_dir/personal/yujing/Yodai/100k/data/distance/25kb/E14/tmp/")

# chrom = "chr19"
# rep = "1"
# save intermidiate result in tmp folder

# clean temporary files
def clean(dir):
    import shutil
    for files in os.listdir(dir):
        path = os.path.join(dir, files)
        try:
            shutil.rmtree(path)
        except OSError:
            os.remove(path)

def main():

    df = dd.read_csv( str(input_dir / f"*rep{rep}*_update.csv"),  usecols= ["name", "chrom", "g", "x_um", "y_um", "z_um", dbscan_col, "fov", "cellID"],  dtype = dtype_dict)
    # select the chromorosome
    df = df[(df.chrom == chrom) & (df[dbscan_col] != -1)]
    dist_df = df.groupby(["fov", "cellID", dbscan_col]).apply(calc_intra_chrom_distance, meta={'g1': 'int32', 'g2': 'int32', 'dist_um': 'float16'}).reset_index(drop = True)
    tmp_dist = dist_df.groupby(["g1", "g2"])["dist_um"].apply(lambda x: x.values, meta=('dist_um', 'f8')).compute()
    clean(python_temp_dir)
    tmp_dist = tmp_dist.reset_index()

    # look for paints
    tmp_dist["p1"] = tmp_dist["g1"].map(lambda x: g2p[x])
    tmp_dist["p2"] = tmp_dist["g2"].map(lambda x: g2p[x])
    tmp_dist["p_dist"] = tmp_dist["p2"] - tmp_dist["p1"]
    tmp_dist["g_dist"] = tmp_dist["g2"] - tmp_dist["g1"]

    # norm_dfs = []
    # for g_dist in tmp_dist["g_dist"].unique():
    #     norm_dfs.append(gdist_correct(tmp_dist, g_dist))

    # # concat the norm_dfs
    # norm_df = pd.concat(norm_dfs)
    #save the norm_df
    pd.to_pickle(tmp_dist, output_dir / f"{chrom}_{chrom}_rep{rep}_25kb.pkl")

if __name__ == "__main__":
    input_dir, chrom, output_dir, rep, celltype = read_input()

    # set the TMPDIR environment variable to the path of the new temporary directory to try to avoid OS 39 error
    python_temp_dir = f'/central/scratch/yy4ng/tmp/{celltype}_{chrom}_{rep}'
    Path(python_temp_dir).mkdir(parents=True, exist_ok=True)
    os.environ['TMPDIR'] = python_temp_dir

    # get the celltype
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    # default values that stays constant
    dbscan_col = "dummy_dbscan"
    # default values that stays constant
    dtype_dict = {'name': str, 'chrom': str, 'g' : int, 
                  'x_um': float, 'y_um': float, 'z_um': float,
                  dbscan_col: int, 'fov': int, 'cellID': int}

    data_root = Path("/home/yy4ng/group_dir/personal/yujing/Yodai/100k/data/distance")
    coords = pd.read_csv(os.path.join(data_root, "annot", "LC1-100k-0516-mm10-25kb-GC-repeat.csv"))[["name", "chrom", "start", "end"]]
    paints = pd.read_csv(os.path.join(data_root, "annot", "2020-09-20-NewBalance-loci-paint-barcodes.csv"))[["name", "paint"]]
    coords = coords.merge(paints, how = "left")
    # replace row_idx to real coordinates
    # replace g with real coordinates
    coord_df = coords[coords["chrom"] == chrom].dropna()
    coord_df.loc[:, "g"] = coords["name"].str.split("-").str[-1].astype(int)
    g2p = {g:p for g, p in zip(coord_df['g'], coord_df['paint'])}
    # run the main function
    main()
    print ("Done! XD")