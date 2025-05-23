#/bin/bash
# default paramters
python_path=/central/home/yy4ng/anaconda3/bin/python3.8
bin_tmp=./python_src/binning_gamma.py

cell=$1
input_path=$2
chrom=$3
output_path=$4
log_dir=$5

tmp_script=$(mktemp ${log_dir}/tmp-script.XXXXXX)
# save the intermediate result in tmp_path
cat > "$tmp_script" << EOF
#!/bin/bash
#SBATCH --job-name=${chrom}_tmp
#SBATCH --output=${log_dir}/slurm-$chrom-rep$rep.out
#SBATCH --error=${log_dir}/slurm-$chrom-rep$rep.err
#SBATCH --time 6:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=40G

${python_path} -u ${bin_tmp} --chrom ${chrom} --input_dir $input_path  --output_dir $output_path --celltype $cell > ${log_dir}/binning-$cell-$chrom-rep$rep.log 2>&1
EOF
# Submit the job
sbatch "$tmp_script"
