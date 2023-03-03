# Two-layer DNASeqFISH+ Downstream Analyze
***
Two-layer DNASeqFISH+ is a method high-resolution, multi-omics profiling 
<br>
This repository contains experimental materials, scripts for imaging processing, barcode calling and processed data used in paper: link
***
# Organization to get start
| Directory | Subdirectory               | Contents                                                                                                  |
|-----------|----------------------------|-----------------------------------------------------------------------------------------------------------|
| data      | experimental-resources     | Probe sequences, barcoding scheme and antibody information                                                |
| data      | annotation                 | Gene annotation, sequence features of DNA locus at different resolutions                                  |
| data      | CellCulcure and cerebellum | Processed ensemble and single cell level data together with resource data for figures in this paper       |
| scripts   | image_processing           | Code used for imaging preprocessing, image alignment, dot calling and decoding (Matlab)                   |
| scripts   | pairwise_distance          | Scripts for calculating physical distance between detected DNA locus (Julia and Python)                   |
| scripts   | 3D_visualization           | Code and test dataset used for 3D visualization of chromosomes and subcompartment reconstruction (Python) |
| scripts   | other                      | Code for replicating downstream analysis (Python)                                                         |

***
## Image processing
Process raw images from Two-layer DNASeqFISH experiment, together with mRNA and DNA probe barcoding scheme to output physical coordinate of each detected individual DNA loci, mRNA or intron dots.
> ### Dependencies:
> 1.Matlab (R2019a)
> <br>
> 2.Fiji installed ()
> ### Requirements:
> 
> 
> <br>
> Usage refer to readme document under image_processing folder

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
* Contact the corresponding author: Long Cai (lcai@caltech.edu) for any inquiry.





