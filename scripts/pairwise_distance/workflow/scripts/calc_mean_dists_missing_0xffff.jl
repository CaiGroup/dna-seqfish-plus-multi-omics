using Distributed
@everywhere using DataFrames
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
for path_num in 2:length(snakemake.input)
    global pnts = vcat(pnts,get_pnts(snakemake.input[path_num], allele_col))
end

pnts[:,"finalcellID"] .= Int16(0)
pnts_cell_grps = groupby(pnts, [:fov, :cellID])
ncells = length(pnts_cell_grps)
for (finalcellid, cell_pnts) in enumerate(pnts_cell_grps)
    cell_pnts[:, "finalcellID"] .= Int16(finalcellid)
end
pnts = DataFrame(pnts_cell_grps)
pnts = pnts[:,Not(:cellID)]
println("added finalcellID")

# for testing
#filter!(locus -> locus.finalcellID in [1,2], pnts)
#ncells=2
pnts[!, "chrom"] = map(elem -> elem[4] == 'X' ? 20 : parse(Int64, elem[4:end]), pnts.chrom)

sort!(pnts, [:chrom, :finalcellID, :allele, :g])
pnts_grpd = groupby(pnts, :chrom)
chrom_names = unique(pnts.chrom)
nchroms = length(chrom_names)

println(snakemake.threads)
extra_procs = snakemake.threads
println("given threads: $extra_procs")
println("nprocs init: ", nprocs())
addprocs(extra_procs)
println("nprocs after add: ", nprocs())

@everywhere begin
    using DataFrames
    using InlineStrings

    function get_chrm_chrm_dists(chr1, chr2, ncells)
        max_chr1 = maximum(chr1.g)
        max_chr2 = maximum(chr2.g)

        dists = zeros(Float64, max_chr1, max_chr2)
        nsums = zeros(Float64, max_chr1, max_chr2)
        for cell in 1:ncells
            println("cell: ", cell)
            cell_dists, cell_nsums = get_cell_chrom_dists(chr1[chr1.finalcellID .== cell,:], chr2[chr2.finalcellID .== cell,:], max_chr1, max_chr2)
            dists .+= cell_dists
            nsums .+= cell_nsums
        end
        dists[nsums .== 0] .= 0xffff
        nsums[nsums .== 0] .= 1
        return dists ./ nsums
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
        dists = get_chrm_chrm_dists(inputs[1], inputs[2], inputs[3])
	    dists = round.(UInt16, dists)
        return dists, inputs[4], inputs[5]
    end
end

#do pmap
inputs = []
chrm_dfs = collect(pnts_grpd)
for i in 1:nchroms, j in i:nchroms
    push!(inputs, (copy(chrm_dfs[i]),  copy(chrm_dfs[j]), ncells, i, j))
end
println("starting pmap")
mats = pmap(get_dist_mat, inputs)
println("pmapped")
rmprocs(collect(2:(extra_procs+1))...)

#piece back together
ng_chroms = DataFrame(map(chrm -> (chrm.chrom[1], maximum(chrm.g)), collect(pnts_grpd)))
ng_chroms[:,"end_ind"] = cumsum(ng_chroms[:,2])
ng_chroms[:,"start_ind"] = ng_chroms[:, "end_ind"] .- ng_chroms[:,2] .+ 1
println(ng_chroms)

img = zeros(UInt16, ng_chroms.end_ind[end], ng_chroms.end_ind[end])
#im_nsums = zeros(UInt32, ng_chroms.end_ind[end], ng_chroms.end_ind[end])
for pair_mat in mats
    dists, i, j = pair_mat
    println("Chrms: ", i, ", ", j)
    global img[ng_chroms.start_ind[i]:ng_chroms.end_ind[i], ng_chroms.start_ind[j]:ng_chroms.end_ind[j]] .= round.(UInt16, dists)
    #im_nsums[ng_chroms.start_ind[i]:ng_chroms.end_ind[i], ng_chroms.start_ind[j]:ng_chroms.end_ind[j]] .= round.(UInt32, nsums)
    if i != j
        global img[ng_chroms.start_ind[j]:ng_chroms.end_ind[j], ng_chroms.start_ind[i]:ng_chroms.end_ind[i]] .= round.(UInt16, dists')
        #im_nsums[ng_chroms.start_ind[j]:ng_chroms.end_ind[j], ng_chroms.start_ind[i]:ng_chroms.end_ind[i]] .= round.(UInt32, nsums')
    end
end


println("saving...")
save(snakemake.output[1], img)
#save(snakemake.output[2], im_nsums)
println("done!")
