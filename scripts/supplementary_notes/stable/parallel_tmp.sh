#/bin/bash

for chrom_id in {1..19}; 
do
    for rep_id in {1..2};
    do
    # run bash file
    bash cerebellum_distance_pipeline.sh MLI1 /groups/CaiLab/personal/yujing/Yodai/100k/data/raw/MLI1/profile_per_dot/radial_organization /groups/CaiLab/personal/yujing/Yodai/100k/data/distance_final chr$chrom_id $rep_id
    done
done

bash cerebellum_distance_pipeline.sh MLI1 /groups/CaiLab/personal/yujing/Yodai/100k/data/raw/MLI1/profile_per_dot/radial_organization /groups/CaiLab/personal/yujing/Yodai/100k/data/distance_final chrX 1

bash cerebellum_distance_pipeline.sh MLI1 /groups/CaiLab/personal/yujing/Yodai/100k/data/raw/MLI1/profile_per_dot/radial_organization /groups/CaiLab/personal/yujing/Yodai/100k/data/distance_final chrX 2

