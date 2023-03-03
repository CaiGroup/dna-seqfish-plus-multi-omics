using DataFrames
using CSV
using FileIO
using Images

"""
Reads points, throws out unused columns to save memory
"""
function read_pnts(path, allele_col)
    pnts_ = DataFrame(CSV.File(path))
    pnts = pnts_[:, ["cellID", "x","y","z","chrom", "g", allele_col]]
    pnts.cellID = Int16.(pnts.cellID)
    rename!(pnts, allele_col => "allele")
    pnts.allele = Int8.(pnts.allele)
    filter!(locus -> locus.allele != -1, pnts)
    pnts[!,"x"] .*= 103 #nm/pixel
    pnts[!,"y"] .*= 103 #nm/pixel
    pnts[!,"z"] .*= 250 #nm/slice
    if occursin("rep", path)
        rep_substr_pos = findfirst("rep", path)
        rep = parse(Int64, path[rep_substr_pos[end]+1])
        pnts[!,"rep"] .= rep
    end
    return pnts
end

allele_col = "dbscan_ldp_nbr_allele"

function get_pnts(file_path, allele_col)
    fov_pnts = read_pnts(file_path, allele_col)
    csv_filename = split(file_path, "/")[3]
    pos = split(split(csv_filename, "pos_")[2], "_")[1]
    fov = parse(Int16, pos)
    fov_pnts[:,"fov"] .= fov
    return fov_pnts
end

println("read in")
pnts = get_pnts(snakemake.input[2], allele_col)
#pnts[!,"rep"] .= parse(Int64, split(snakemake.input[2],"/")[3][10])
println(first(pnts, 5))
for path_num in 3:length(snakemake.input)
    new_pnts = get_pnts(snakemake.input[path_num], allele_col)
#    new_pnts[!, "rep"] .= parse(Int64, split(snakemake.input[path_num], "/")[3][10])
    global pnts = vcat(pnts,new_pnts)
#    println("nrow(pnts): ", nrow(pnts))
end
println("nrow(pnts): ", nrow(pnts))
println(first(pnts, 5))

clusters_of_interest = [[2],[4],[7],[6,11]]


filter!(locus -> locus.chrom == snakemake.wildcards["coi"], pnts)

cell_info = DataFrame(CSV.File(snakemake.input[1]))
cell_info = cell_info[:,Not([:x, :y, :z])]
println("cell info:")
println(first(cell_info, 5))
println("pnts: ")
println(first(pnts, 5))
pnts = innerjoin(pnts, cell_info, on=[:rep, :fov, :cellID])
filter!(row -> row.doublet != 1, pnts)

pnts[:,"finalcellID"] .= Int16(0)
pnts_cell_grps = groupby(pnts, [:fov, :cellID, :rep])
ncells = length(pnts_cell_grps)
for (finalcellid, cell_pnts) in enumerate(pnts_cell_grps)
    cell_pnts[:, "finalcellID"] .= Int16(finalcellid)
end
pnts = DataFrame(pnts_cell_grps)
pnts = pnts[:,Not(:cellID)]
println("added finalcellID")

println(first(pnts, 5))

sort!(pnts, [:chrom, :finalcellID, :allele, :g])
pnts_grpd = groupby(pnts, :chrom)
chrom_names = unique(pnts.chrom)
nchroms = length(chrom_names)

function get_chrm_chrm_dists(chr1, chr2)
    max_chr1 = maximum(chr1.g)
    max_chr2 = maximum(chr2.g)
    dists = zeros(Float64, max_chr1, max_chr2)
    nsums = zeros(Float64, max_chr1, max_chr2)
    strands = groupby(chr1, [:finalcellID, :allele])
    for strand in strands
        cell_dists, cell_nsums = get_cell_chrom_dists(strand, strand, max_chr1, max_chr2)
        dists .+= cell_dists
        nsums .+= cell_nsums
    end
    dists[nsums .== 0] .= 0xffff
    nsums[nsums .== 0] .= 1
    return dists ./ nsums, nsums
end

function get_cell_chrom_dists(chr1, chr2, max1, max2)
    dists = zeros(Float64, max1, max2)
    nsums= zeros(Float64, max1, max2)
    temp = Float64[0,0,0]
    dist = 0.0
    for locus1 in eachrow(chr1), locus2 in eachrow(chr2)
        @inbounds temp[1] = locus1["x"] - locus2["x"]
        @inbounds temp[2] = locus1["y"] - locus2["y"]
        @inbounds temp[3] = locus1["z"] - locus2["z"]
        temp .= temp .^2
        dist = sum(temp)
        dist = sqrt(dist)
        dists[locus1.g, locus2.g] += dist
        nsums[locus1.g, locus2.g] += 1
    end
    return dists, nsums
end

for (i, clusters) in enumerate(clusters_of_interest)
    println("starting computation for ", snakemake.output[i]) 
    println("cluster: ", clusters)
    cluster_pnts = filter(locus -> locus.leiden in clusters, pnts)
    println(first(cluster_pnts, 5))
    dists, nsums = get_chrm_chrm_dists(cluster_pnts, cluster_pnts)

    if maximum(dists) < 2^16
        dists = round.(UInt16, dists)
    else
        dists = round.(UInt32, dists)
    end

    println("saving..." , snakemake.output[i])
    save(snakemake.output[i], dists)
    println("done!")
end
