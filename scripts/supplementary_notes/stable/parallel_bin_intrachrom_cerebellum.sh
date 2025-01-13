#/bin/bash
for chrom_id in {1..19};
do
    bash binning_cerebellum_intrachrom_pipeline.sh MLI2_PLI /home/yy4ng/group_dir/personal/yujing/Yodai/100k/data/distance_final/tmp/MLI2_PLI/ chr$chrom_id /home/yy4ng/group_dir/personal/yujing/Yodai/100k/data/distance_final/original/MLI2_PLI/ /home/yy4ng/group_dir/personal/yujing/Yodai/100k/data/distance_final/log/MLI2_PLI
done

bash binning_cerebellum_intrachrom_pipeline.sh MLI2_PLI /home/yy4ng/group_dir/personal/yujing/Yodai/100k/data/distance_final/tmp/MLI2_PLI/ chrX /home/yy4ng/group_dir/personal/yujing/Yodai/100k/data/distance_final/original/MLI2_PLI/ /home/yy4ng/group_dir/personal/yujing/Yodai/100k/data/distance_final/log/MLI2_PLI


