import pandas as pd
import numpy as np
import dask.dataframe as dd
import dask.config
from pathlib import Path
import fastparquet as fp
import concurrent.futures
from concurrent.futures import ThreadPoolExecutor
import multiprocessing
from utils import *

# default values that stays constant
dbscan_col = "dbscan_ldp_allele"
dtype_dict = {'name': str, 'chrom': str, 'g' : int, 'x_um': float, 'y_um': float, 'z_um': float,
              dbscan_col: int, 'fov': int, 'cellID': int}


# accept inputs from the user
def read_input():
    import argparse

    # get the input variables from the command line using argparse
    parser = argparse.ArgumentParser(description='Calculate the distance for each bin pair')
    parser.add_argument('--chrom', type=str, help='the chromosome number you want to check')
    parser.add_argument('--input_dir', type=str, help='the directory of the input files')
    parser.add_argument('--block_size', type=int, help='chunk size when breaking down the dataframe')
    parser.add_argument('--output_dir', type=str, help='the directory of the output files as median final result format')
    parser.add_argument('--rep', type=str, help='which replicate to choose')

    # parse the arguments
    args = parser.parse_args()
    
    input_dir = args.input_dir
    chrom = args.chrom
    block_size = args.block_size
    output_dir = args.output_dir
    rep = args.rep

    return input_dir, chrom, int(block_size), output_dir, rep

def main():
    input_dir, chrom, block_size, output_dir, reps = read_input()
    
    print ("chrom", chrom)
    print ("rep", reps)
    print ("block size", block_size)
    print ("output", output_dir)
    print ("input_dir", input_dir)

    # convert the input_dir and output_dir to Path object
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)

    chrom1, chrom2 = chrom.split("&")
    print (f" start calculating distances between {chrom1}, {chrom2}")

    # start computing distance
    dist_df_list = []
    # compute rep1 and rep2 separately while retaining distance information for combined result
    for rep in reps.split("&"):
        print ("Start calculating distance matrix of replicate", rep)
        df = dd.read_csv(str(input_dir / f"*replicate{rep}*.csv"), usecols= ["name", "chrom", "g", "x_um", "y_um", "z_um", dbscan_col, "fov", "cellID"], dtype = dtype_dict)
        # select the chromorosome
        df = df[(df.chrom.isin([chrom1, chrom2])) & (df[dbscan_col] != -1)]
        dist_df = df.groupby(["fov", "cellID"]).apply(calc_inter_chrom_distance, meta={chrom1:  'int32', chrom2 :  'int32', 'dist_um': 'float16'}).reset_index(drop = True).compute()
        dist_df_list.append(dist_df)
        print (f"distance dataframe for {rep} is created, delete the original dataframe")
        # remove useless data and release memory
        del df

    # after save the result for further reading, first break down the dataframe into chunks
    # get g_blocks
    g_max1, g_max2 = dist_df_list[0][chrom1].max(), dist_df_list[0][chrom2].max()
    g1_block_start = np.arange(1, g_max1 + 1, block_size)
    g1_block_end = g1_block_start + block_size - 1
    g2_block_start = np.arange(1, g_max2 + 1, block_size)
    g2_block_end = g2_block_start + block_size - 1

    combinations = [(g1_block_start[i], g1_block_end[i], g2_block_start[j], g2_block_end[j]) for i in range(len(g1_block_start)) for j in range(len(g2_block_start))]

    #
    # parallelize the dataframe breakdown and median calculation

    NUM_WORKERS = multiprocessing.cpu_count() - 1

    for i, dist_df in enumerate(dist_df_list):
        print (f"start processing {i}")
        result = []
        with ThreadPoolExecutor(max_workers=NUM_WORKERS) as executor:
            futures = [executor.submit(chunk_interchrom_df_parallel, dist_df, comb , chrom1, chrom2) for comb in combinations]
            # Wait for all tasks to complete
            for future in concurrent.futures.as_completed(futures):
                # Handle any potential exceptions from the tasks
                try:
                    res = future.result() # get the result from the future
                    if res is not None:
                        result.append(res)
                except Exception as e:
                    print(f"An error occurred: {e}")


        # remove useless data and release memory
        print ("done calculating median pairwise distance, keep the distance dataframe for all reps")
        # concatenate results and sort, then save it
        final_df = pd.concat(result).sort_values([chrom1, chrom2])
        final_df[chrom1] = final_df[chrom1].astype(np.int32)
        final_df[chrom2] = final_df[chrom2].astype(np.int32)
        final_df["dist_um"] = final_df["dist_um"].astype(np.float16)
        fp.write(str(output_dir / f"{chrom1}_{chrom2}_rep{i+1}_25kb.parquet"), final_df, compression = "GZIP")

    # calculate all reps
    dist_df = pd.concat(dist_df_list)
    print ("start processing all reps, delete the list container and release memory")
    del dist_df_list
    result = []
    with ThreadPoolExecutor(max_workers=NUM_WORKERS) as executor:
        futures = [executor.submit(chunk_interchrom_df_parallel, dist_df, comb , chrom1, chrom2) for comb in combinations]
        # Wait for all tasks to complete
        for future in concurrent.futures.as_completed(futures):
            # Handle any potential exceptions from the tasks
            try:
                future.result()
                if res is not None:
                    result.append(res)
            except Exception as e:
                print(f"An error occurred: {e}")

    final_df = pd.concat(result).sort_values([chrom1, chrom2])
    final_df[chrom1] = final_df[chrom1].astype(np.int32)
    final_df[chrom2] = final_df[chrom2].astype(np.int32)
    final_df["dist_um"] = final_df["dist_um"].astype(np.float16)
    fp.write(str(output_dir / f"{chrom1}_{chrom2}_rep{reps}_25kb.parquet"), final_df, compression = "GZIP")
    print ("Complete calculation is done :)")

# start main function
if __name__ == "__main__":
    main()



