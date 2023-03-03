rule bin_25k_to_200k:
    input: 
        "resources/filtered_renumbered_dots/{experiment}_pos_{pos}_filtered_renumbered.csv",
        "resources/LC1-100k-0516-mm10-25kb-GC-repeat.csv"
    output: "results/filtered_renumbered_dots_binned/{experiment}_pos_{pos}_filtered_renumbered.csv"
    resources: mem_mb=2000
    script: "../scripts/bin_25k_loci_to_200k.py"

rule bin_25k_to_200k_brain:
    input: 
        "resources/cerebellum/Brain_rep{rep}_pos_{pos}_radial_scores_wofilter.csv",
        "resources/LC1-100k-0516-mm10-25kb-GC-repeat.csv"
    output: "results/cerebellum_binned/Brain_rep{rep}_pos_{pos}_radial_scores_wofilter.csv"
    resources: mem_mb=2000
    script: "../scripts/bin_brain_25k_loci_to_200k.py"

rule bin_paint:
    input:
        "resources/filtered_renumbered_dots/{experiment}_pos_{pos}_filtered_renumbered.csv",
        "resources/2020-09-20-NewBalance-loci-paint-barcodes.csv"
    output:
        "results/filtered_renumbered_dots_paint_binned/{experiment}_pos_{pos}_filtered_renumbered.csv"
    resources: mem_mb=2000
    script: "../scripts/bin_paint.py"
rule bin_paint_brain:
    input:
        "resources/cerebellum/Brain_rep{rep}_pos_{pos}_radial_scores_wofilter.csv",
        "resources/2020-09-20-NewBalance-loci-paint-barcodes.csv"
    output:
        "results/cerebellum_paint_binned/Brain_rep{rep}_pos_{pos}_radial_scores_wofilter.csv"
    resources: mem_mb=2000
    script: "../scripts/bin_brain_paint.py"
