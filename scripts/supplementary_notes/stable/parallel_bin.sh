#/bin/bash

for chrom_id in {1..5}; 
do
    for rep_id in {1..2};
    do
    # run bash file
    bash binning_pipeline.sh E14 "/groups/CaiLab/personal/yujing/Yodai/100k/data/distance_final" chr$chrom_id $rep_id  "/groups/CaiLab/personal/yujing/Yodai/100k/data/distance_final/log/E14"
    done
done

#bash binning_pipeline.sh E14 "/groups/CaiLab/personal/yujing/Yodai/100k/data/distance_final" chrX 1  "/groups/CaiLab/personal/yujing/Yodai/100k/data/distance_final/log/E14"
#bash binning_pipeline.sh E14 "/groups/CaiLab/personal/yujing/Yodai/100k/data/distance_final" chrX 2  "/groups/CaiLab/personal/yujing/Yodai/100k/data/distance_final/log/E14"

