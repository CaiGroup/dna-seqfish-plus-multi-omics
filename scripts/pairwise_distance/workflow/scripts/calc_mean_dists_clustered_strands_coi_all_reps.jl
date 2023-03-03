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
    pnts.cellID = Int8.(pnts.cellID)
    rename!(pnts, allele_col => "allele")
    pnts.allele = Int8.(pnts.allele)
    filter!(locus -> locus.allele != -1, pnts)
    pnts[!,"x"] .*= 103 #nm/pixel
    pnts[!,"y"] .*= 103 #nm/pixel
    pnts[!,"z"] .*= 250 #nm/slice
    if occursin("rep", path)
        rep_substr_pos = findfirst("rep", path)
        rep = path[rep_substr_pos[end]+1]
        pnts[!,"rep"] .= rep
    else
        pnts[!,"rep"] .= 1
    end
    return pnts
end

allele_col = "dbscan_ldp_nbr_allele_sorted"

function get_pnts(file_path, allele_col)
    fov_pnts = read_pnts(file_path, allele_col)
    csv_filename = split(file_path, "/")[3]
    pos = split(split(csv_filename, "pos_")[2], "_")[1]
    fov = parse(Int8, pos)
    fov_pnts[:,"fov"] .= fov
    return fov_pnts
end

println("read in")
pnts = get_pnts(snakemake.input[1], allele_col)
for path_num in 2:(length(snakemake.input)-1)
    global pnts = vcat(pnts,get_pnts(snakemake.input[path_num], allele_col))
end
println(first(pnts, 5))
filter!(locus -> locus.chrom == snakemake.wildcards["coi"], pnts)
println(first(pnts, 5))
pnts[:,"finalcellID"] .= Int16(0)
pnts_cell_grps = groupby(pnts, [:fov, :cellID], :rep])
ncells = length(pnts_cell_grps)
for (finalcellid, cell_pnts) in enumerate(pnts_cell_grps)
    cell_pnts[:, "finalcellID"] .= Int16(finalcellid)
end
pnts = DataFrame(pnts_cell_grps)
println("added finalcellID")

# read clusters
clusters = DataFrame(CSV.File(snakemake.input[length(snakemake.input)]))
batch_2_exp_dict = Dict(1 => "E14_mDuxCa_wo", 2 => "E14_mDuxCA_24hr_rep1", 3 => "E14_mDuxCA_24hr_rep2")
println("filter experiment/batch")
println(first(clusters, 5))
filter!(row -> batch_2_exp_dict[row.batch] == snakemake.params["experiment"], clusters)
println(first(clusters, 5))
println("join")
println(first(pnts,5))
pnts = innerjoin(pnts, clusters, on=[:fov, :cellID])
println(first(pnts,5))
println("filter cluster")
cluster = parse(Int64,snakemake.wildcards["cluster"])
println("typeof(cluster): ", typeof(cluster))
println("cluster: ", cluster)
println("typeof(pnts.leiden): ", typeof(pnts.leiden))
println(unique(pnts.leiden))
filter!(row -> row.leiden == cluster, pnts)
println(first(pnts, 5))
println("filtered batch cluster")

pnts = pnts[:,Not(:cellID)]
sort!(pnts, [:chrom, :finalcellID, :allele, :g])
pnts_grpd = groupby(pnts, :chrom)
chrom_names = unique(pnts.chrom)
println("chrom_names: ", chrom_names)
nchroms = length(chrom_names)

function get_chrm_chrm_dists(chr1, chr2, ncells)
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

function get_dist_mat(inputs)
    println(inputs[1].chrom[1], ", ", inputs[2].chrom[1])
    dists, nsums = get_chrm_chrm_dists(inputs[1], inputs[2], inputs[3])
    dists = round.(UInt16, dists)
    return dists, inputs[4], inputs[5]#, nsums
end

println("do pmap")
#do pmap
inputs = []
chrm_dfs = collect(pnts_grpd)
for i in 1:nchroms, j in i:nchroms
    push!(inputs, (copy(chrm_dfs[i]),  copy(chrm_dfs[j]), ncells, i, j))
end
#println("inputs: ", inputs)
mats = map(get_dist_mat, inputs)
println("typeof(mats): ", typeof(mats))
println(mats)
println("typeof(mats[1]): ", typeof(mats[1]))

println("saving...")
save(snakemake.output[1], mats[1][1])
println("done!")
