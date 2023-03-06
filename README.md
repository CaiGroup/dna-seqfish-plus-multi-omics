# Spatial multi-omics data analysis
***
DNA seqFISH+ is a method for high-resolution, multi-omics profiling in cell culture and complex tissues.
<br>
This repository contains experimental materials, scripts for imaging processing and downstream analysis.
***
# Organization to get start
| Directory | Subdirectory               | Contents                                                                                                  |
|-----------|----------------------------|-----------------------------------------------------------------------------------------------------------|
| data      | experimental-resources     | Probe sequences, barcoding scheme and antibody information                                                |
| data      | annotation                 | Gene annotation, sequence features of DNA locus at different resolutions                                  |
| data      | CellCulcure and cerebellum | Processed ensemble and single cell level data together with resource data for figures in this paper       |
| scripts   | image_processing           | Code used for imaging preprocessing, image alignment, dot calling and decoding (Matlab and Python)        |
| scripts   | pairwise_distance          | Scripts for calculating physical distance between detected DNA locus (Julia and Python)                   |
| scripts   | 3D_visualization           | Code and test dataset used for 3D visualization of chromosomes and subcompartment reconstruction (Python) |
| scripts   | downstream_analysis        | Code for downstream analysis (Python)                                                                     |

***
## Image processing
Process raw images from DNA seqFISH+ experiment, together with RNAseqFISH and sequential immunofluorescence images.
> ### Dependencies:
> 1. Matlab (R2019a)
> 2. Python 3.8
> 3. Fiji installed
> <br>
> This pipeline is built upon previous image processing pipeline (https://github.com/CaiGroup/dna-seqfish-plus-tissue).

***
## Pairwise distance calculation
Input physical coordinates of DNA locus after homologous chromosome separation, output ensemble level pairwise distance of genomic bins in selected resolution.
> ### Dependencies:
> 1. Julia
> 2. Snakemake
> <br>
>Usage refer to readme document under pairwise_distance folder

***
## 3D Visualization:
Input single cell DNA locus coordinates, meta information and immunoflourescence to reconstruct 3D cell visualization
> ### Dependencies
> 1. Python 3.8
> 2. Mayavi 4.8.1
><br>
>Usage refer to readme document uder 3D_visualization folder

## License
Free for non-commercial and academic research. The software is bound by the licensing rules of California Institute of Technology (Caltech)

## Contact
* Contact the corresponding authors: Long Cai (lcai@caltech.edu) or Yodai Takei (ytakei@caltech.edu) for any inquiry.





