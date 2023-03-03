using Distributed
@everywhere using DataFrames
using CSV
using FileIO
using Images
using SparseArrays

"""
Reads points, throws out unused columns to save memory
"""
function read_pnts(path, allele_col)
    pnts_ = DataFrame(CSV.File(path))
    pnts = pnts_[:, ["cellID", "x","y","z","chrom", "g", "name", allele_col]]
    pnts.cellID = Int8.(pnts.cellID)
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
pnts_cell_grps = groupby(pnts, [:rep, :fov, :cellID])
ncells = length(pnts_cell_grps)
for (finalcellid, cell_pnts) in enumerate(pnts_cell_grps)
    cell_pnts[:, "finalcellID"] .= Int16(finalcellid)
end
pnts = DataFrame(pnts_cell_grps)
pnts = pnts[:,Not(:cellID)]
println("added finalcellID")

sort!(pnts, [:chrom, :finalcellID, :allele, :g])
cell_grps = groupby(pnts, :finalcellID)

chrm_gmax = combine(groupby(pnts, :chrom), :g => maximum => :g_max)
chrm_gmax[!, "cum_g"] .= cumsum(chrm_gmax[!, "g_max"])
chrm_start_m1_dict = Dict(collect(zip(chrm_gmax[:, :chrom], [[0]; chrm_gmax[1:end-1, :cum_g]])))
ind_max = chrm_gmax[end, :cum_g]

println(snakemake.threads)
extra_procs = snakemake.threads
println("given threads: $extra_procs")
println("nprocs init: ", nprocs())
addprocs(extra_procs)
println("nprocs after add: ", nprocs())

@everywhere begin
    using DataFrames
    using InlineStrings
    using NearestNeighbors
    using SparseArrays

    get_strand_pw(x :: Tuple{DataFrame, Int64, Dict{String7, Int64}, Int64}) = get_strand_pw(x...)

    function get_strand_pw(points, r,  chrm_start_m1_dict, ind_max)
        println("cell: ", points.finalcellID[1])
        freq_norm = spzeros(UInt16, ind_max, ind_max)
        locus_freq = combine(groupby(points, [:chrom, :g]), nrow)
        for l1 in eachrow(locus_freq), l2 in eachrow(locus_freq)
            l1_ind = chrm_start_m1_dict[l1.chrom] + l1.g
            l2_ind = chrm_start_m1_dict[l2.chrom] + l2.g
            if l1.g == l2.g 
                freq_norm[l1_ind, l2_ind] = l1.nrow*l2.nrow - l1.nrow
            else
                freq_norm[l1_ind, l2_ind] = l1.nrow*l2.nrow
            end
        end
        return freq_norm
    end

end

#do pmap
inputs = []
cell_dfs = collect(cell_grps)
ncells = length(cell_dfs)
for i in 1:ncells
    push!(inputs, (copy(cell_dfs[i]), 500, chrm_start_m1_dict, ind_max))
end

co_occurances = pmap(get_strand_pw, inputs)
all_co_occurances = zeros(UInt16, ind_max, ind_max)
for co_occurance in co_occurances
    all_co_occurances .+= co_occurance
end

#rmprocs(collect(2:51)...)

println("saving...")
save(snakemake.output[1], all_co_occurances)
#println("saving down sampled version...")
#save("contacts_plot_500nm_"*source_dir*"_downsampled_ncells_$ncells"*".png", ceil.(UInt16, imresize(all_contacts, (30000, 30000))))
println("done!")

