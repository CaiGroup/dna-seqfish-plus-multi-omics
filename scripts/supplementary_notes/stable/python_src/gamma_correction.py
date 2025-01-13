import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
import dask.dataframe as dd
from lmfit.models import Model
from scipy.stats import gamma, percentileofscore

# This will suppress all warnings
warnings.simplefilter(action='ignore', category=Warning)

import os
from concurrent.futures import ProcessPoolExecutor

cpu_count = os.cpu_count()
print(cpu_count)

# try another way
# wrap up gamma correction
# correct for each g_dist

def split_into_groups(pdf, col):
    size = len(pdf)
    print(f"Processing partition of size: {size}")
    
    # Create a dictionary where keys are the unique values of 'col' and values are the sub-dataframes
    groups = {key: sub_df for key, sub_df in pdf.groupby(col)}
    return groups

    
def concat_and_return_with_index(group):
    """
    Concatenate 'dist_um' values in the provided group.
    
    Parameters:
    -----------
    group : DataFrame
        Data group with a 'dist_um' column.
        
    Returns:
    --------
    ndarray
        Concatenated values of 'dist_um' column from the group.
    """
    concatenated_values = np.hstack(group["dist_um"])
    return concatenated_values

def retrieve_distribution(dist_container):
    """
    Retrieve distribution from the provided container.
    
    Parameters:
    -----------
    dist_container : Series/DataFrame
        Container holding distributions.
        
    Returns:
    --------
    tuple
        Index of the container (p_dist values) and the values of the distribution.
    """
    return dist_container.index, dist_container.values

def prepare_distance(sub_df):
    """
    Prepare distance information by grouping based on 'p_dist' and concatenating 'dist_um'.
    
    Parameters:
    -----------
    sub_df : DataFrame
        DataFrame with 'p_dist' column for grouping and 'dist_um' column for data extraction.
        
    Returns:
    --------
    tuple
        Indexes and concatenated values based on 'p_dist' grouping.
    """
    return retrieve_distribution(sub_df.groupby("p_dist").apply(concat_and_return_with_index))

def fit_contatenate_gamma(*samples):
    s = np.concatenate(*samples)
    # Convert raw data to histogram
    y, edges = np.histogram(s, bins=100, density=True)
    x = (edges[:-1] + edges[1:]) / 2  # bin centers

    # Get initial guess from scipy's gamma fitting
    alpha_init, loc_init, scale_init = gamma.fit(s, floc=0)

    # Define the gamma distribution model
    def gamma_model(x, alpha, loc, scale):
        return gamma.pdf(x, a=alpha, loc=loc, scale=scale)

    gmodel = Model(gamma_model)

    # Initial parameter values and bounds using scipy's results
    params = gmodel.make_params(alpha=alpha_init, loc=loc_init, scale=scale_init)
    params['alpha'].set(min=0.01, max=10)
    params['loc'].set(min=-10, max=10)
    params['scale'].set(min=0.01, max=10)

    # Perform the fitting
    result = gmodel.fit(y, x=x, params=params)
    return [result.params['alpha'].value, result.params['loc'].value, result.params['scale'].value]

def calc_corr_df(sub_df):
    # get the distance and corresponding p_dist
    ps, ss = prepare_distance(sub_df)

    # fit the general gamma distribution
    params = fit_contatenate_gamma(ss)

    # grab the sub_df and correct for median
    # calculate median first
    sub_df["dist_median"] = sub_df["dist_um"].apply(lambda x : np.median(x))

    # map medians back to general gamma distribution
    # to do this more efficiently, try to use vectorize then mapping back
    # split the dist_median into groups by p_dist
    corr_dfs = []
    for p, orig_dist in zip(ps,ss) :
        ss_df = sub_df[sub_df["p_dist"] == p].copy()
        # calculate percent
        orig_median = ss_df["dist_median"].values
        prct = percentileofscore(orig_dist, orig_median)
        ss_df["gamma_median"] = gamma.ppf(prct / 100, a=params[0], loc = params[1], scale=params[2])
        corr_dfs.append(ss_df)
    return pd.concat(corr_dfs)

def map_value_to_gamma(x, params_orig, params_new):
    # Compute the CDF value (percentile) of x in the original gamma distribution
    p = gamma.cdf(x, a=params_orig[0], loc=params_orig[1], scale=params_orig[2])
    
    # Find the corresponding value in the new gamma distribution
    mapped_value = gamma.ppf(p, a=params_new[0], loc=params_new[1], scale=params_new[2])
    
    return mapped_value

def calc_corr_cdf_single_allele_df(sub_df):
    # get the distance and corresponding p_dist
    ps, ss = prepare_distance(sub_df)

    # fit the general gamma distribution
    params_all = fit_contatenate_gamma(ss)
    # grab the sub_df and correct for median
    # calculate corresponding values for each column
    corr_dfs = []
    for p, orig_dist in zip(ps,ss) :
        ss_df = sub_df[sub_df["p_dist"] == p].copy()
        # calculate the original params
        params_orig = fit_contatenate_gamma([orig_dist])
        # calculate percent
        ss_df["norm_um"] = ss_df["dist_um"].apply(lambda x: map_value_to_gamma(x, params_orig, params_all))
        corr_dfs.append(ss_df)
    return pd.concat(corr_dfs)

def driver_function(sub_df):
    g = sub_df["g_dist"].iloc[0]
    print (f"start {g}")
    final_df =  calc_corr_cdf_single_allele_df(sub_df) 
    print (f"end {g}")
    return final_df

# read in parameters
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

# read in parameters
def main():
    input_dir, chrom, output_dir, rep, celltype = read_input()

    # input_dir = "/home/yy4ng/group_dir/personal/yujing/Yodai/100k/data/distance_final/tmp/E14"
    # output_dir = "/home/yy4ng/group_dir/personal/yujing/Yodai/100k/data/distance_final/tmp_gamma/E14"
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    # chrom = "chr19"
    # rep = 1
    # celltype = "E14"

    # get the tmp dataframe
    df = pd.read_pickle(str(input_dir/f"{chrom}_{chrom}_rep{rep}_25kb.pkl"))

    ddf = dd.from_pandas(df[['g1', 'g2', 'dist_um', 'p1', 'p2', 'p_dist', 'g_dist']], npartitions=50)  # Adjust npartitions based on your needs

    # release memory
    del df

    # Apply the function to every partition
    result = ddf.map_partitions(split_into_groups, col='g_dist', meta=object).compute()

    # Merge sub-dataframes having the same key
    merged_groups = {}
    for d in result:
        for key, sub_df in d.items():
            if key in merged_groups:
                # Stitch dataframes having the same key
                merged_groups[key] = pd.concat([merged_groups[key], sub_df], axis=0)
            else:
                merged_groups[key] = sub_df

    # Convert dictionary values to a list to get the final list of sub-dataframes
    sub_dfs = list(merged_groups.values())

    # release memory
    del ddf

    # try multiplexing
    with ProcessPoolExecutor(max_workers=cpu_count - 1) as executor:
        # The map method blocks until all results are available
        tdf_list = list(executor.map(driver_function, sub_dfs))

    del sub_dfs
    # save df
    final_df = pd.concat(tdf_list)
    final_df.to_pickle(str(output_dir/f"{chrom}_{chrom}_rep{rep}_25kb.pkl") )
    print ("Done XD")

if __name__ == "__main__":
    main()