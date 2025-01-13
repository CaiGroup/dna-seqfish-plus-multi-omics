import pandas as pd
import numpy as np
from pathlib import Path
import os
from scipy.sparse import csr_matrix, save_npz


def format_mtx(df):
    arr = df.values
    rows, cols, vals = arr[:, 0].astype(int), arr[:, 1].astype(int), arr[:, 2]
    s = max(np.max(rows), np.max(cols)) + 1
    matrix = csr_matrix( (vals, (rows, cols)), shape = (s, s))
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
    parser.add_argument('--output_dir', type=str, help='the directory of the output files as median final result format')
    parser.add_argument('--rep', type=str, help='which replicate to choose')
    # parser.add_argument('--max_bin', type=str, help='the maximum bin number')

    # parse the arguments
    args = parser.parse_args()
    
    input_dir = args.input_dir
    chrom = args.chrom
    # block_size = args.block_size
    output_dir = args.output_dir
    # rep = args.rep
    celltype = args.celltype
    # max_bin = args.max_bin

    return input_dir,  output_dir, chrom, celltype

def main():

    reps = [1, 2]
    ress = [25000, 50000, 100000, 200000, 500000, 1000000]

    # make the matrix
    data_root = Path("/home/yy4ng/group_dir/personal/yujing/Yodai/100k/data/distance")
    coords = pd.read_csv(os.path.join(data_root, "annot", "LC1-100k-0516-mm10-25kb-GC-repeat.csv"))[["name", "chrom", "start", "end"]]
    paints = pd.read_csv(os.path.join(data_root, "annot", "2020-09-20-NewBalance-loci-paint-barcodes.csv"))[["name", "paint"]]
    coords = coords.merge(paints, how = "left")
    coord_df = coords[coords["chrom"] == chrom]
    coord_df.loc[:, "g"] = coords["name"].str.split("-").str[-1].astype(int)
    translation_table = coord_df[["g", "start"]].set_index("g").to_dict()["start"]

    # calculte, using an example
    # use the small chromosome
    df1 = pd.read_pickle(str(input_dir  /  f"{chrom}_{chrom}_rep{reps[0]}_25kb.pkl")).sort_values(by=['g1', 'g2'])
    df2 = pd.read_pickle(str(input_dir  /  f"{chrom}_{chrom}_rep{reps[1]}_25kb.pkl")).sort_values(by=['g1', 'g2'])

    for df, rep in zip([df1, df2], reps):
        df["s1"] = df["g1"].map(translation_table)
        df["s2"] = df["g2"].map(translation_table)
        for res in ress:
            df["b1"] = df["s1"] // res
            df["b2"] = df["s2"] // res
            norm  = format_mtx(df.groupby(["b1", "b2"])["norm_um"].apply(lambda x: np.median(np.hstack(x))).reset_index())
            # save 
            save_npz(corr_output_dir / f'{chrom}_rep{rep}_{res}.npz', norm)
    # grab all reps
    df = pd.concat([df1, df2])
    del df1
    del df2
    df["s1"] = df["g1"].map(translation_table)
    df["s2"] = df["g2"].map(translation_table)

    for res in ress:
        df["b1"] = df["s1"] // res
        df["b2"] = df["s2"] // res
        norm  = format_mtx(df.groupby(["b1", "b2"])["norm_um"].apply(lambda x: np.median(np.hstack(x))).reset_index())
        # save 
        save_npz(corr_output_dir / f'{chrom}_repall_{res}.npz', norm)
    
if __name__ == "__main__":  
    input_dir,  output_dir, chrom, celltype = read_input()
    input_dir = Path(input_dir)
    corr_output_dir = Path(output_dir)
    main()
    print ("done")
    #     chrom = "chr6", celltype = "E14", input_dir = "/groups/CaiLab/personal/yujing/Yodai/100k/data/distance_final/tmp_gamma/E14"
    # # corr_output_dir = "/groups/CaiLab/personal/yujing/Yodai/100k/data/distance_final/corrected/E14"