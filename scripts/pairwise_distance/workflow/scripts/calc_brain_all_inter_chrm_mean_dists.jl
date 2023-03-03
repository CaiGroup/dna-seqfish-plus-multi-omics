using DataFrames
using CSV
using FileIO
using Images
using Statistics
using DelimitedFiles

get_chr_pw_dist_mat_from_cell(df) = [sqrt(sum((Vector(df[i, [:x_mean, :y_mean, :z_mean]]) .- Vector(df[j, [:x_mean, :y_mean, :z_mean]])).^2)) for i in 1:20, j in 1:20]

cell_info = DataFrame(CSV.File(snakemake.input[1]))
cell_info = cell_info[:,Not([:x, :y, :z])]

clusters_of_interest = [[2],[4],[7],[6,11]]

fnames = [e[2] for e in snakemake.input if e[1] != 1 ]

for (i, clusters) in enumerate(clusters_of_interest)
    println("clusters: ", clusters)
    fov_pw_dist_mats = []
    for fname in fnames
        println("fname: $fname")
        pnts = DataFrame(CSV.File(fname))
        pnts[!,"x"] .*= 103 #nm/pixel
        pnts[!,"y"] .*= 103 #nm/pixel
        pnts[!,"z"] .*= 250 #nm/slice
        pnts[!, "chrom"] = map(elem -> elem[4] == 'X' ? 20 : parse(Int64, elem[4:end]), pnts.chrom)
        sort!(pnts, :chrom)

        if occursin("rep", fname)
            rep_substr_pos = findfirst("rep", fname)
            rep = parse(Int64, fname[rep_substr_pos[end]+1])
            pnts[!,"rep"] .= rep
        end
        pnts = innerjoin(pnts, cell_info, on=[:rep, :fov, :cellID])
        filter!(row -> row.doublet != 1, pnts)
        cluster_pnts = filter(locus -> locus.leiden in clusters, pnts)
        chrm_centroids = combine(groupby(cluster_pnts, [:cellID, :chrom]), [:x => mean, :y => mean, :z => mean])
        sort!(chrm_centroids, :chrom)
        cell_grps = groupby(chrm_centroids, :cellID)
        filtered_cell_grps = filter(cell_grp -> nrow(cell_grp) == 20, cell_grps)
        pw_dist_mats = map(get_chr_pw_dist_mat_from_cell, collect(filtered_cell_grps))
        push!(fov_pw_dist_mats, pw_dist_mats)
    end

    pw_dist_mats = vcat(fov_pw_dist_mats...)

    mean_pw_dist_mat = mean(pw_dist_mats)

    open(snakemake.output[i], "w") do io
        writedlm(io, mean_pw_dist_mat, ',')
    end

end