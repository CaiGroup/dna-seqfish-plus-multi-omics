import pandas as pd

loci_fname = snakemake.input[0]
bin_key_fname = snakemake.input[1]

loci = pd.read_csv(loci_fname)
loci = loci.loc[:,['fov', 'cellID', 'x', 'y','z','name','chrom','dbscan_ldp_nbr_allele']]

bin_key = pd.read_csv(bin_key_fname)
bin_key = bin_key.loc[:, ["name", "200kb name", "200kb bin", "start", "end"]]

loci.set_index("name", inplace=True)
bin_key.set_index("name", inplace=True)

joined = loci.join(bin_key)

joined["g"] = joined["200kb bin"]

joined.to_csv(snakemake.output[0])
