import pandas as pd
import numpy as np
import dask.dataframe as dd
import dask.config
from pathlib import Path
from utils import *
import os


# clean temporary files
def clean(dir):
    import shutil
    for files in os.listdir(dir):
        path = os.path.join(dir, files)
        try:
            shutil.rmtree(path)
        except OSError:
            os.remove(path)
            
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

# # try the longest chromosome chr1 with all replicates chromosome 1


def main():
    
    data_root = Path("/home/yy4ng/group_dir/personal/yujing/Yodai/100k/data/distance")
    
    coords = pd.read_csv(os.path.join(data_root, "annot", "LC1-100k-0516-mm10-25kb-GC-repeat.csv"))[["name", "chrom", "start", "end"]]
    paints = pd.read_csv(os.path.join(data_root, "annot", "2020-09-20-NewBalance-loci-paint-barcodes.csv"))[["name", "paint"]]
    coords = coords.merge(paints, how = "left")
    coords.loc[:, "g"] = coords["name"].str.split("-").str[-1].astype(int)
    # replace row_idx to real coordinates
    # replace g with real coordinates
    coord_df1 = coords[coords["chrom"] == chrom1].dropna()
    coord_df2 = coords[coords["chrom"] == chrom2].dropna()
    g2p_chrom1 = {g:p for g, p in zip(coord_df1['g'], coord_df1['paint'])}
    g2p_chrom2 = {g:p for g, p in zip(coord_df2['g'], coord_df2['paint'])}
    g2s_chrom1 = {g:s for g, s in zip(coord_df1['g'], coord_df1['start'])}
    g2s_chrom2 = {g:s for g, s in zip(coord_df2['g'], coord_df2['start'])}


    df = dd.read_csv(str(input_dir / f"*rep{rep}*_update.csv"), usecols= ["name", "chrom", "g", "x_um", "y_um", "z_um", dbscan_col, "fov", "cellID"], dtype = dtype_dict)
    # select the chromorosome
    df = df[(df.chrom.isin([chrom1, chrom2])) & (df[dbscan_col] != -1)]
    df["chrom_id"] = df["chrom"].str.slice(start=3).astype(int)
    df = df.sort_values(by = "chrom_id")
    dist_df = df.groupby(["fov", "cellID"]).apply(calc_inter_chrom_distance, meta={chrom1:  'int32', chrom2 :  'int32', 'dist_um': 'float16'}).reset_index(drop = True)
    tmp_dist = dist_df.groupby([chrom1, chrom2])["dist_um"].apply(lambda x: x.values, meta=('dist_um', 'float')).compute().reset_index()
    tmp_dist["p1"] = tmp_dist[chrom1].map(g2p_chrom1)
    tmp_dist["p2"] = tmp_dist[chrom2].map(g2p_chrom2)
    tmp_dist["s1"] = tmp_dist[chrom1].map(g2s_chrom1)
    tmp_dist["s2"] = tmp_dist[chrom2].map(g2s_chrom2)

    pd.to_pickle(tmp_dist, output_dir / f"{chrom1}_{chrom2}_rep{rep}_25kb.pkl")
    
if __name__ == "__main__":
    
    # default values that stays constant
    dbscan_col = "dummy_dbscan"
    dtype_dict = {'name': str, 'chrom': str, 'g' : int, 'x_um': float, 'y_um': float, 'z_um': float,
                  dbscan_col: int, 'fov': int, 'cellID': int}
    
    input_dir, chrom, output_dir, rep, celltype = read_input()
    
    print (chrom)
    print (chrom.split("-"))
    chrom1, chrom2 = chrom.split("-")
    
    # save intermidiate result in tmp folder
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    
    # set the TMPDIR environment variable to the path of the new temporary directory to try to avoid OS 39 error
    python_temp_dir = f'/central/scratch/yy4ng/tmp/{celltype}_{chrom}_{rep}'
    Path(python_temp_dir).mkdir(parents=True, exist_ok=True)
    os.environ['TMPDIR'] = python_temp_dir

    main()
    
    clean(python_temp_dir)
    print ("Done, XD")
    
