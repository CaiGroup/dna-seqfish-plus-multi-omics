for chrom_id1 in {2..9};
do
    for chrom_id2 in {10..19};
    do
        comb=chr$chrom_id1"-chr"$chrom_id2
        for rep in {1..2};
        do
             bash cerebellum_interchrom_distance_pipeline.sh Purkinje /groups/CaiLab/personal/yujing/Yodai/100k/data/raw/Purkinje/profile_per_dot/radial_organization /groups/CaiLab/personal/yujing/Yodai/100k/data/distance_final $comb $rep
             bash cerebellum_interchrom_distance_pipeline.sh MLI1 /groups/CaiLab/personal/yujing/Yodai/100k/data/raw/MLI1/profile_per_dot/radial_organization /groups/CaiLab/personal/yujing/Yodai/100k/data/distance_final $comb $rep
             bash cerebellum_interchrom_distance_pipeline.sh MLI2_PLI /groups/CaiLab/personal/yujing/Yodai/100k/data/MLI2_PLI/Purkinje/profile_per_dot/radial_organization /groups/CaiLab/personal/yujing/Yodai/100k/data/distance_final $comb $rep
#             bash cerebellum_interchrom_distance_pipeline.sh Granule /groups/CaiLab/personal/yujing/Yodai/100k/data/raw/Granule/profile_per_dot/radial_organization /groups/CaiLab/personal/yujing/Yodai/100k/data/distance_final $comb $rep
#             bash cerebellum_interchrom_distance_pipeline.sh Bergmann /groups/CaiLab/personal/yujing/Yodai/100k/data/raw/Bergmann/profile_per_dot/radial_organization /groups/CaiLab/personal/yujing/Yodai/100k/data/distance_final $comb $rep
        done
    done
done
