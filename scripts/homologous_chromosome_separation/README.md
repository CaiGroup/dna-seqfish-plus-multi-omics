# DNA SeqFISH Chromosome Allele Separation Workflow

This repository separates different copies of the same chromosome from DNASeqFISH data into alleles using [DNASeqFISHChromsomeAssignment.jl](https://github.com/CaiGroup/DNASeqFISHChromosomeAssignment).

## Installation

Install the latest version of DNASeqFISHChromsomeAssignment.jl and other necessary julia packages using the commands

```
julia
julia> ]
pkg> add https://github.com/CaiGroup/DNASeqFISHChromosomeAssignment.jl
pkg> add DataFrames
pkg> add CSV
pkg> add GLPK
```


If you have not already installed snakemake, you can do so by following commands from the snakemake full installation [instructions](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

Make sure your snakemake conda environment is active

```
conda activate snakemake
```

## Set up

If using on the Caltech HPC, first run the `make_softlinks.sh` script to make softlinks to the datafiles to chromosome copy numbers to using the command

```
bash make_softlinks.sh
```

In the snakefile, you can adjust the parameters for different datasets.

* r_dbscan - the [DBSCAN](https://en.wikipedia.org/wiki/DBSCAN) search radius
* r_ldp - LDP search radius
* s - the LDP spatial distance penalty scaling factor
* min_size - the minimum allowed size of of a DBSCAN cluster or LDP
* r_ldp_nbr - the search radius within which to add unclassified neighboring points to LDP clusters
* unique_prop_thresh - the genomic loci uniqueness threshold below which DBSCAN clusters are split into two LDPs
* dbscan_min_nbrs - the minimum number of neighbors within a [DBSCAN](https://en.wikipedia.org/wiki/DBSCAN) search radius a point must have to be a core point

## Running the Workflow

To check what files the snakemake needs to process, run

```
snakemake -n
```

To start the workflow on the HPC, run

```
sbatch run-slurm.sh
```

To run locally, run

```
snakemake -c<n>
```
where `<n>` is the number of cores you would like to provide.

plt_cell_chrm_pnts_strnds.py contains a python function to help plot the results.
