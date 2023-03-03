# Two-layer DNASeqFISH+ pairwise_maps

This workflow calculates pairwise distance maps from tables of chromosome separation Two-layer DNASeqFISH+ loci.

## Installing packages

To install [snakemake](https://snakemake.readthedocs.io/en/stable/index.html), use the following commands from the snakemake full installation [instructions](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html): 

<pre> <code>
conda install -n base -c conda-forge mamba
conda activate base
mamba create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake
</code> </pre>

You will need to activate your snakemake conda environment whenever you run the workflow.

Before running the workflow, it is necessary to install the Julia pacakges:

DataFrames.jl

CSV.jl

Images.jl

FileIO.jl

DNASeqFISHChromosomeAssignment.jl

You can install these packages by opening a julia session and typing the commands

```
>>using Pkg
>>Pkg.add("CSV")
>>Pkg.add("DataFrames")
>>Pkg.add("Image")
>>Pkg.add("FileIO")
>>Pkg.add(https://github.com/CaiGroup/DNASeqFISHChromosomeAssignment)
```

You will need to install the python packages in your snakemake anaconda environment:

pandas

numpy

matplotlib

scikit-image

PIL

## Input Data

Data for embryonic stem cells are placed in a folder

```resources/filtered_renumbered_dots```

and data for cerebellum is placed in a folder

```resources/LC1-100k-0516-mm10-25kb-GC-repeat.csv```

Other files containing information on clustering and binning configurations are placed directly in the resources folder.

## Running

From the command line in the main workflow directory, run

```
snakemake --use-conda -c<n>
```

where ```<n>``` is the number of cores you would like to use for running the workflow.
