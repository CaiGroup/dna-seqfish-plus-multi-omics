using DataFrames
using CSV
using DNASeqFISHChromosomeAssignment
#try
    #using CPLEX
    #opt = CPLEX.Optimizer
#catch e
using GLPK
opt = GLPK.Optimizer
#end

#using CPLEX
#opt = CPLEX.Optimizer

prm = ChromSepParams()


set_r_dbscan(prm, snakemake.params["r_dbscan"])
set_r_ldp(prm, snakemake.params["r_ldp"])
set_r_ldp_nbr(prm, snakemake.params["r_dbscan"])
set_sigma(prm, snakemake.params["s"])
#max_strands = snakemake.params["max_strands"]
#max_x_strands = snakemake.params["xmax_strands"]
set_min_size(prm, snakemake.params["min_size"])
set_minProp_unique(prm, snakemake.params["unique_prop_thresh"])
set_dbscan_min_nbrs(prm, snakemake.params["dbscan_min_nbrs"])


dbscan_r_min = 400
dbscan_r_max = 800
dbscan_r_inc = 50
overlap_thresh = 0.05

filename = snakemake.input[1]
pnts = DataFrame(CSV.File(filename))
filter!(:chrom => !=("control"), pnts)

pnts[!,"x"] .*= 103 #nm/pixel
pnts[!,"y"] .*= 103 #nm/pixel
pnts[!,"z"] .*= 250 #nm/slice
pnts[!,"locusID"] = Array(1:nrow(pnts))


pnts[!,"g"] = [parse(Int64, split(nm, "-")[2]) for nm in pnts[!,"name"]]

res = assign_chromosomes(pnts, prm, opt)

#res[!,"x"] .= xp #1/103 #pixel/nm
#res[!,"y"] .= yp #1/103 #pixel/nm
#res[!,"z"] .= zs #1/250 #slice/nm

CSV.write(snakemake.output[1], res)
