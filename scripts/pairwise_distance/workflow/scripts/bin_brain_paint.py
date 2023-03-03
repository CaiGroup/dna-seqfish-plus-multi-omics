import pandas as pd

loci_fname = snakemake.input[0]
bin_key_fname = snakemake.input[1]

loci = pd.read_csv(loci_fname)
loci = loci.loc[:,['fov', 'cellID', 'x', 'y','z','name','chrom','dbscan_ldp_nbr_allele']]

bin_key = pd.read_csv(bin_key_fname)
bin_key_chrm_min = bin_key.groupby("chrom").min()
bin_key = bin_key.loc[:, ["name", "paint"]]

loci.set_index("name", inplace=True)
bin_key.set_index("name", inplace=True)

joined = loci.join(bin_key)

g_paint = [r.paint - bin_key_chrm_min.loc[r.chrom, "paint"] + 1 for (i,r) in joined.iterrows()]

joined["g"] = g_paint

joined.to_csv(snakemake.output[0])
