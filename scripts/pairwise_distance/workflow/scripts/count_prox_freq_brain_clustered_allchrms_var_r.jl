using DataFrames
using CSV
using FileIO
using Images
using SparseArrays
using NearestNeighbors

"""
Reads points, throws out unused columns to save memory
"""
function read_pnts(path, allele_col)
    pnts_ = DataFrame(CSV.File(path))
    pnts = pnts_[:, ["cellID", "x","y","z","chrom", "g", allele_col]]
    #pnts.cellID = Int8.(pnts.cellID)
    rename!(pnts, allele_col => "allele")
    #pnts.allele = Int8.(pnts.allele)
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

search_radius=parse(Int64, snakemake.wildcards["r"])

allele_col = "dbscan_ldp_nbr_allele"

function get_pnts(file_path, allele_col)
    fov_pnts = read_pnts(file_path, allele_col)
    csv_filename = split(file_path, "/")[3]
    pos = split(split(csv_filename, "pos_")[2], "_")[1]
    fov = parse(Int8, pos)
    fov_pnts[:,"fov"] .= fov
    return fov_pnts
end

println("read in")
pnts = get_pnts(snakemake.input[2], allele_col)
for path_num in 3:length(snakemake.input)
    global pnts = vcat(pnts,get_pnts(snakemake.input[path_num], allele_col))
end

clusters_of_interest = [[2],[4],[7],[6,11]]


filter!(locus -> locus.chrom == snakemake.wildcards["coi"], pnts)

cell_info = DataFrame(CSV.File(snakemake.input[1]))
cell_info = cell_info[:,Not([:x, :y, :z])]

pnts = innerjoin(pnts, cell_info, on=[:rep, :fov, :cellID])
filter!(row -> row.doublet != 1, pnts)


pnts[:,"finalcellID"] .= Int16(0)
pnts_cell_grps = groupby(pnts, [:fov, :cellID])
ncells = length(pnts_cell_grps)
for (finalcellid, cell_pnts) in enumerate(pnts_cell_grps)
    cell_pnts[:, "finalcellID"] .= Int16(finalcellid)
end
pnts = DataFrame(pnts_cell_grps)
pnts = pnts[:,Not(:cellID)]
println("added finalcellID")

sort!(pnts, [:chrom, :finalcellID, :allele, :g])

max_g = maximum(pnts.g)

@everywhere begin
    using DataFrames
    using InlineStrings
    using NearestNeighbors
    using SparseArrays

    get_pnts_mat(points) = Matrix(Matrix(points[:, ["x", "y", "z"]])')

    get_cell_contacts(x :: Tuple{DataFrame, Int64, Dict{String15, Int64}}) = get_cell_contacts(x...)
    
    function get_cell_contacts(points, radius, names_inds_dict)
        println("cell: ", points.finalcellID[1])
        points_mat = get_pnts_mat(points)
        tree = KDTree(points_mat)
        pnts_nbrs = inrange(tree, points_mat, radius)
        nloci = length(names_inds_dict)
        contacts = spzeros(UInt16, nloci, nloci)

        for (pnt_ind, pnt_nbrs) in enumerate(pnts_nbrs)
            g_pnt = points.name[pnt_ind]
            g_pnt_ind = names_inds_dict[g_pnt]
            for nbr_ind in pnt_nbrs
                if pnt_ind != nbr_ind
                    g_nbr = points.name[nbr_ind]
                    g_nbr_ind = names_inds_dict[g_nbr]
                    contacts[g_pnt_ind, g_nbr_ind] += 0x0001
                end
            end
        end

        return contacts
    end
end

#do pmap


for (i, clusters) in enumerate(clusters_of_interest)
    println("starting computation for ", snakemake.output[i]) 
    println("cluster: ", clusters)
    cluster_pnts = filter(locus -> locus.leiden in clusters, pnts)
    println(first(cluster_pnts, 5))
    cell_grps = groupby(cluster_pnts, :finalcellID)

    inputs = []
    cell_dfs = collect(cell_grps)
    ncells = length(cell_dfs)
    for i in 1:ncells
        push!(inputs, (copy(cell_dfs[i]), 500, name_2_ind_dict))
    end

    mats = pmap(get_cell_contacts, inputs)
    all_contacts = zeros(UInt16, nloci, nloci)
    for mat in mats
        all_contacts .+= mat
    end
    mats = []

    co_occurances = pmap(get_strand_pw, inputs)
    all_co_occurances = zeros(UInt16, max_g, max_g)
    for co_occurance in co_occurances
        all_co_occurances .+= co_occurance
    end


    save(snakemake.output[i], all_contacts)

    save(snakemake.output[i+length(clusters_of_interest)], all_norm)

end
