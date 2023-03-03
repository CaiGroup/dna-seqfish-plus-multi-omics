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
for path_num in 2:(length(snakemake.input)-1)
    global pnts = vcat(pnts,get_pnts(snakemake.input[path_num], allele_col))
end

filter!(locus -> locus.chrom == snakemake.wildcards["coi"], pnts)

pnts[:,"finalcellID"] .= Int16(0)
pnts_cell_grps = groupby(pnts, [:fov, :cellID])
ncells = length(pnts_cell_grps)
for (finalcellid, cell_pnts) in enumerate(pnts_cell_grps)
    cell_pnts[:, "finalcellID"] .= Int16(finalcellid)
end

pnts = DataFrame(pnts_cell_grps)
pnts = pnts[:,Not(:cellID)]
println("added finalcellID")

# read clusters
clusters = DataFrames(CSV.File(snakemake.input[end]))
batch_2_exp_dict = Dict(1 => "E14_mDuxCa_wo", 2 => "E14_mDuxCA_24hr_rep1", 3 => "E14_mDuxCA_24hr_rep2")
filter!(row -> batch_2_exp_dict[row.batch] == "E14_mDuxCA_24hr_rep1", clusters)
pnts = innerjoin(pnts, clusters, on=[:fov, :cellID])
filter!(row -> row.leiden == snakemake.wildcards["cluster"], pnts)

sort!(pnts, [:chrom, :finalcellID, :allele, :g])
strand_grps = groupby(pnts, [:finalcellID, :allele])

max_g = maximum(pnts.g)

get_pnts_mat(points) = Matrix(Matrix(points[:, ["x", "y", "z"]])')

get_cell_contacts(x :: Tuple{DataFrame, Int64, Int64}) = get_cell_contacts(x...)
    
function get_cell_contacts(points, radius, max_g)
    println("cell: ", points.finalcellID[1])
    points_mat = get_pnts_mat(points)
    tree = KDTree(points_mat)
    pnts_nbrs = inrange(tree, points_mat, radius)
    contacts = spzeros(UInt16, max_g, max_g)
    for (pnt_ind, pnt_nbrs) in enumerate(pnts_nbrs)
        g_pnt = points.g[pnt_ind]
        for nbr_ind in pnt_nbrs
            if pnt_ind != nbr_ind
                g_nbr = points.g[nbr_ind]
                contacts[g_pnt, g_nbr] += 0x0001
            end
        end
    end
    return contacts
end

get_strand_pw(x :: Tuple{DataFrame, Int64, Int64}) = get_strand_pw(x...)

function get_strand_pw(points, r,  max_g)
    println("cell: ", points.finalcellID[1])
    freq_norm = zeros(UInt16, max_g, max_g)
    locus_freq = combine(groupby(points, :g), nrow)
    for l1 in eachrow(locus_freq), l2 in eachrow(locus_freq)
        if l1.g == l2.g 
            freq_norm[l1.g, l2.g] = l1.nrow*l2.nrow - l1.nrow
        else
            freq_norm[l1.g, l2.g] = l1.nrow*l2.nrow
        end
    end
    return freq_norm
end

#do pmap
inputs = []
strand_dfs = collect(strand_grps)
nstrands = length(strand_dfs)
all_contacts = zeros(UInt16, max_g, max_g)
all_norm = zeros(UInt16, max_g, max_g)
for i in 1:nstrands
    all_contacts .+= get_cell_contacts(strand_dfs[i], 500, max_g)
    all_norm .+= get_strand_pw(strand_dfs[i], 500, max_g)
end

mats = map(get_cell_contacts, inputs)
norms = map(get_strand_pw, inputs)


println("saving...")
save(snakemake.output[1], all_contacts)

save(snakemake.output[2], all_norm)

println("done!")
