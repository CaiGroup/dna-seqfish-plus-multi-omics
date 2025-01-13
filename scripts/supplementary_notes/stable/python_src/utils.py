
import pandas as pd
import numpy as np


def calc_intra_chrom_distance(sub_df):
    """
    This function calculates the intra-chromosomal distances based on the 3D genomic physical coordinates.
    
    Parameters:
    sub_df (DataFrame): A pandas DataFrame containing the 3D genomic physical coordinates, along with the bin values.
    
    Returns:
    DataFrame: A new DataFrame with the calculated distances between different genomic bins.
    """
    
    # import nessecary packages
    from scipy.spatial.distance import pdist
    sub_df = sub_df.sort_values(by = "g")
    # get the 3D genomic physical coordinates and the bin value in the sub_df
    xyz = sub_df[["x_um", "y_um", "z_um"]].values
    xyz = (xyz - np.mean(xyz, axis=0)).astype(np.float16)
    g = sub_df["g"].values
    
    # calculate the 1D pdist distance 
    dist = pdist(xyz).astype(np.float16)
    
    # calculate the indexing pair of bins
    n = len(g)
    rows, cols = np.triu_indices(n, k=1)  # k=1 to exclude the diagonal
    # Replace each index in index_pairs with the corresponding value from g
    # print
    if np.all(g[rows] == g[cols]):
        print (sub_df)

    return pd.DataFrame( {"g1" : g[rows], "g2" : g[cols], "dist_um": dist}) 


def calc_inter_chrom_distance(sub_df):
    from scipy.spatial.distance import pdist

    sub_df = sub_df.sort_values(by = ["chrom_id", "g"])

    c = sub_df["chrom"].values
    g = sub_df["g"].values

    # mean center xyz
    xyz = sub_df[["x_um", "y_um", "z_um"]].values
    xyz = (xyz - np.mean(xyz, axis=0)).astype(np.float16)
    # calculate pdist
    # calculate the 1D pdist distance 
    dist = pdist(xyz).astype(np.float16)
    # calculate the indexing pair of bins
    n = len(g)
    rows, cols = np.triu_indices(n, k=1)  # k=1 to exclude the diagonal

    # get the coordinate system
    c_rows = c[rows]
    c_cols = c[cols]

    g_rows = g[rows]
    g_cols = g[cols]

    mask = c_rows != c_cols
    # filter out intra chromosome interaction
    c_rows_filtered = c_rows[mask]
    c_cols_filtered = c_cols[mask]
    g_rows_filtered = g_rows[mask]
    g_cols_filtered = g_cols[mask]
    dist_filtered = dist[mask]
    if mask.sum() != 0:
        # get the chrom
        fdf = pd.DataFrame({
            c_rows_filtered[0] : g_rows_filtered,
            c_cols_filtered[0] : g_cols_filtered,        
            "dist_um" : dist[mask]
        })
        # return a dataframe
        return fdf
    else:
        print ("An error occured that this cell does not have one or all chromosomes")
        

# correct for distance
def correct_samples(*samples):
    # Calculate the size-weighted mean and standard deviation
    samples = [s.astype("float64") for s in samples]
    pool_mean = np.average([np.mean(sample) for sample in samples], weights=[len(sample) for sample in samples])
    pool_std = np.sqrt(np.average([np.var(sample, ddof=1) for sample in samples], weights=[len(sample) for sample in samples]))

    # Create a list to store the corrected samples
    corrected_samples = []

    # Standardize, shift, and scale each sample
    for sample in samples:
        standardized = (sample - np.mean(sample)) / np.std(sample, ddof=1)
        corrected = standardized * pool_std + pool_mean
        corrected_samples.append(corrected.astype("float16"))
    return corrected_samples

# apply correction to each g_dist
def gdist_correct(df, select_g):
    print (f"start {select_g}")
    sub_df = df[df["g_dist"] == select_g].copy()
    corrected_dist = []
    grouped_sub = sub_df.groupby(["p_dist"])["dist_um"].apply(lambda x : np.hstack(x)).reset_index()
    ss = grouped_sub["dist_um"].to_list()
    ps = grouped_sub["p_dist"].to_list()
    cc = correct_samples(*ss)

    look_up_dist = {}
    for p, c, s in zip(ps, cc, ss): 
        look_up_table = dict(zip(s, c))
        look_up_dist[p] = look_up_table  

    for i, r in sub_df.iterrows():
        lookup_table = look_up_dist[r["p_dist"]]
        def lookup_func(input_val):
            return lookup_table[input_val]
        vectorized_lookup_func = np.vectorize(lookup_func)
        x = r["dist_um"]
        corrected_dist.append(vectorized_lookup_func(x))
    sub_df["norm_um"] = corrected_dist
    print (f"end {select_g}")
    return sub_df

