import pandas as pd
import numpy as np
from pathlib import Path
import os
from scipy.sparse import csr_matrix, save_npz

# read in line

# def bin_distance(df, col, type = "median"):
#     if type == "median":
#         return df[col].apply(lambda x: np.median(x)).values
#     elif type == "mean":
#         return df[col].apply(lambda x: np.mean(x)).values
#     elif isinstance(type, (int, float)) and 0 <= type <= 100:
#         return df[col].apply(lambda x: np.percentile(x, type)).values
#     else:
#         return None

# def bin_probability(df, col, thresh = 0.3):
#     return df[col].apply(lambda x: np.sum(x < thresh) / len(x)).values

def format_mtx_interchrom(df):
    arr = df.values
    rows, cols, vals = arr[:, 0].astype(int), arr[:, 1].astype(int), arr[:, 2]
    s1, s2 = np.max(rows) + 1, np.max(cols) + 1
    matrix = csr_matrix( (vals, (rows, cols)), shape = (s1, s2))
    return matrix

# accept inputs from the user
def read_input():
    import argparse

    # get the input variables from the command line using argparse
    parser = argparse.ArgumentParser(description='Calculate the distance for each bin pair')
    parser.add_argument('--chrom', type=str, help='the chromosome number you want to check')
    parser.add_argument('--input_dir', type=str, help='the directory of the input files')
    parser.add_argument('--celltype', type=str, help='the directory of the input files')
    # parser.add_argument('--block_size', type=int, help='chunk size when breaking down the dataframe')
    # parser.add_argument('--output_dir', type=str, help='the directory of the output files as median final result format')
    parser.add_argument('--rep', type=str, help='which replicate to choose')
    # parser.add_argument('--max_bin', type=str, help='the maximum bin number')

    # parse the arguments
    args = parser.parse_args()
    
    input_dir = args.input_dir
    chrom = args.chrom
    # block_size = args.block_size
    # output_dir = args.output_dir
    rep = args.rep
    celltype = args.celltype
    # max_bin = args.max_bin

    return input_dir, chrom, rep, celltype


def main():
    rs = reps.split("-")
    chrom1, chrom2 = chrom.split("-")
    print (chrom1, chrom2)
    dfs = []
    ress = [ "paint", 1000000, 500000, 200000]

    for rep in rs:
        df = pd.read_pickle(str(input_dir / "tmp"/ celltype /  f"{chrom1}_{chrom2}_rep{rep}_25kb.pkl")).sort_values(by=[chrom1, chrom2])
        print ("done reading")
        dfs.append(df)
        # calculte, using an example
        # use the small chromosome

        for res in ress:

            if res == "paint":
                df["b1"] = df["p1"] - df["p1"].min()
                df["b2"] = df["p2"] - df["p2"].min()
            else:
                df["b1"] = df["s1"] // res
                df["b2"] = df["s2"] // res

            tdf = df.groupby(["b1", "b2"])["dist_um"].apply(lambda x: np.hstack(x)).reset_index()
            tdf["median"] = tdf["dist_um"].apply(lambda x : np.median(x))
            tdf["quantile"] = tdf["dist_um"].apply(lambda x : np.quantile(x, 0.25) ) 
            mmtx = format_mtx_interchrom(tdf[["b1", "b2", "median"]])
            qmtx = format_mtx_interchrom(tdf[["b1", "b2", "quantile"]])
            # save 
            save_npz(orig_output_dir / f'{chrom1}_{chrom2}_rep{rep}_{res}_median.npz', mmtx)
            save_npz(orig_output_dir / f'{chrom1}_{chrom2}_rep{rep}_{res}_quantile.npz', qmtx)
            print (f"Finish {res} {rep} ")

    # end to pool all reps
    df = pd.concat(dfs)
    del dfs
    for res in ress:
        if res == "paint":
            df["b1"] = df["p1"] - df["p1"].min()
            df["b2"] = df["p2"] - df["p2"].min()
        else:
            df["b1"] = df["s1"] // res
            df["b2"] = df["s2"] // res

        tdf = df.groupby(["b1", "b2"])["dist_um"].apply(lambda x: np.hstack(x)).reset_index()
        tdf["median"] = tdf["dist_um"].apply(lambda x : np.median(x))
        tdf["quantile"] = tdf["dist_um"].apply(lambda x : np.quantile(x, 0.25) ) 
        mmtx = format_mtx_interchrom(tdf[["b1", "b2", "median"]])
        qmtx = format_mtx_interchrom(tdf[["b1", "b2", "quantile"]])
        # save 
        save_npz(orig_output_dir / f'{chrom1}_{chrom2}_rep{reps}_{res}_median.npz', mmtx)
        save_npz(orig_output_dir / f'{chrom1}_{chrom2}_rep{reps}_{res}_quantile.npz', qmtx)
        print (f"Finish {res} {reps} ")
        
if __name__ == "__main__":  
    # get the input variables from the command line using argparse
    input_dir, chrom, reps, celltype = read_input()
    input_dir = Path(input_dir)
    # check if the output dir exists
    orig_output_dir = input_dir / "original" / celltype
#     corr_output_dir = input_dir / "corrected" / celltype

    # make the directory if not exists
    orig_output_dir.mkdir(parents=True, exist_ok=True)
#     corr_output_dir.mkdir(parents=True, exist_ok=True)
    rs = reps.split("-")
    main()




