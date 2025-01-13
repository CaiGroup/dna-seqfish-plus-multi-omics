#/bin/bash
for chrom_id1 in {1..18};
do
    for ((chrom_id2 = $chrom_id1+1; chrom_id2 <=19; chrom_id2++));
    do
        comb=chr$chrom_id1"-chr"$chrom_id2
        # run bash file
        bash binning_interchrom_pipeline.sh Granule /groups/CaiLab/personal/yujing/Yodai/100k/data/distance_final $comb 1-2  /groups/CaiLab/personal/yujing/Yodai/100k/data/distance_final/log/Granule
    done
done


