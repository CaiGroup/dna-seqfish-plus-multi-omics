#/bin/bash
# default paramters
python_path=/central/home/yy4ng/anaconda3/bin/python3.8

cell=$1
input_path=$2
output_path=$3
chrom=$4
rep=$5

calc_tmp=./python_src/calc_tmp_dist_cerebellum.py
# calc_orig=./median_original.py
# calc_corr=./median_corrected.py

# make directories for storing file
temp_dir=$output_path/tmp/$cell
orig_dir=$output_path/original/$cell
corr_dir=$output_path/corrected/$cell
log_dir=$output_path/log/$cell

echo $input_path
echo $output_path
# echo $temp_dir
# echo $orig_dir
# make directories
mkdir -p $temp_dir
mkdir -p $orig_dir
mkdir -p $corr_dir
mkdir -p $log_dir

# if calculating intra chrom distance
# Generate the intermediate result, original distance list for each bin pair. Then gaussian correct, save output in tmp_path
# Create a temporary script file
tmp_script=$(mktemp ${log_dir}/tmp-script.XXXXXX)
# save the intermediate result in tmp_path
cat > "$tmp_script" << EOF
#!/bin/bash
#SBATCH --job-name=${chrom}_${rep}_tmp
#SBATCH --output=${log_dir}/slurm-$chrom-rep$rep.out
#SBATCH --error=${log_dir}/slurm-$chrom-rep$rep.err
#SBATCH --time 5:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=20G

${python_path} -u ${calc_tmp} --chrom ${chrom} --input_dir $input_path --output_dir $temp_dir  --rep $rep > ${log_dir}/tmp_${chrom}_${rep}.log
EOF
# Submit the job
sbatch "$tmp_script"






