using DataFrames
using CSV
using FileIO
using Images
using Statistics
using DelimitedFiles

get_chr_pw_dist_mat_from_cell(df) = [sqrt(sum((Vector(df[i, [:x_mean, :y_mean, :z_mean]]) .- Vector(df[j, [:x_mean, :y_mean, :z_mean]])).^2)) for i in 1:20, j in 1:20]


fov_pw_dist_mats = []
for fname in values(snakemake.input)
    println("fname: $fname")
    pnts = DataFrame(CSV.File(fname))
    pnts[!, "chrom"] = map(elem -> elem[4] == 'X' ? 20 : parse(Int64, elem[4:end]), pnts.chrom)
    pnts[!,"x"] .*= 103 #nm/pixel
    pnts[!,"y"] .*= 103 #nm/pixel
    pnts[!,"z"] .*= 250 #nm/slice
    sort!(pnts, :chrom)
    chrm_centroids = combine(groupby(pnts, [:cellID, :chrom]), [:x => mean, :y => mean, :z => mean])
    sort!(chrm_centroids, :chrom)
    cell_grps = groupby(chrm_centroids, :cellID)
    filtered_cell_grps = filter(cell_grp -> nrow(cell_grp) == 20, cell_grps)
    pw_dist_mats = map(get_chr_pw_dist_mat_from_cell, collect(filtered_cell_grps))
    push!(fov_pw_dist_mats, pw_dist_mats)
end

pw_dist_mats = vcat(fov_pw_dist_mats...)

mean_pw_dist_mat = mean(pw_dist_mats)

open(snakemake.output[1], "w") do io
    writedlm(io, mean_pw_dist_mat, ',')
end
