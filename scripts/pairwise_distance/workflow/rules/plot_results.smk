
rule plot_pw_proximities:
    input: "results/prox_freq_outputs/{experiment}_prox_counts_{coi}.png", "results/prox_freq_outputs/{experiment}_ncoocur_{coi}.png"
    output: "results/plots/prox_freq_{experiment}_{coi}.png"
    resources: mem_mb = 1000
    conda: "../envs/plot.yaml"
    script: "../scripts/plot_prox_freq_coi.py"

rule plot_pw_proximities_binned:
    input: "results/prox_freq_outputs_binned/{experiment}_prox_counts_{coi}.png", "results/prox_freq_outputs_binned/{experiment}_ncoocur_{coi}.png"
    output: "results/plots/prox_freq_binned_{experiment}_{coi}.png"
    resources: mem_mb = 1000
    conda: "../envs/plot.yaml"
    script: "../scripts/plot_prox_freq_coi.py"

rule plot_pw_proximities_paint_binned:
    input: "results/prox_freq_outputs_paint_binned/{experiment}_prox_counts_{coi}.png", "results/prox_freq_outputs_binned/{experiment}_ncoocur_{coi}.png"
    output: "results/plots/prox_freq_paint_binned_{experiment}_{coi}.png"
    resources: mem_mb = 1000
    conda: "../envs/plot.yaml"
    script: "../scripts/plot_prox_freq_coi.py"

rule plot_pw_proximities_brain:
    input: "results/prox_freq_outputs/exp_{experiment}_prox_counts_{coi}_r{r}_cluster{cluster}.png",
           "results/prox_freq_outputs/exp_{experiment}_ncoocur_{coi}_r{r}_cluster{cluster}.png"
    output: "results/plots/prox_freq_{experiment}_{coi}_r{r}_cluster{cluster}.png"
    resources: mem_mb = 2000
    conda: "../envs/plot.yaml"
    wildcard_constraints: experiment = "Brain_rep\d+&*\d*"
    script: "../scripts/plot_prox_freq_coi.py"

rule plot_pw_proximities_brain_binned:
    input: "results/prox_freq_outputs_binned/exp_{experiment}_prox_counts_{coi}_r{r}_cluster{cluster}.png",
           "results/prox_freq_outputs_binned/exp_{experiment}_ncoocur_{coi}_r{r}_cluster{cluster}.png"
    output: "results/plots/prox_freq_binned_{experiment}_{coi}_r{r}_cluster{cluster}.png"
    resources: mem_mb = 1000
    conda: "../envs/plot.yaml"
    wildcard_constraints: experiment = "Brain_rep\d+&*\d*"
    script: "../scripts/plot_prox_freq_coi.py"

rule plot_pw_mean_dist:
    input: "results/mean_outputs/exp_{experiment}_mean_dists_{coi}.png"
    output: "results/plots/mean_dists_{experiment}_{coi}.png"
    resources: mem_mb = 10000
    conda: "../envs/plot.yaml"
    script: "../scripts/plot_mean_dist_coi.py"

rule plot_pw_mean_dist_binned:
    input: "results/mean_outputs_binned/exp_{experiment}_mean_dists_{coi}.png"
    output: "results/plots/mean_dists_binned_{experiment}_{coi}.png"
    resources: mem_mb = 10000
    conda: "../envs/plot.yaml"
    script: "../scripts/plot_mean_dist_coi.py"

rule plot_pw_mean_dist_paint_binned:
    input: "results/mean_outputs_paint_binned/exp_{experiment}_mean_dists_{coi}.png"
    output: "results/plots/mean_dists_paint_binned_{experiment}_{coi}.png"
    resources: mem_mb = 10000
    conda: "../envs/plot.yaml"
    script: "../scripts/plot_mean_dist_coi.py"

rule plot_pw_mean_strand_clusterdist:
    input: "results/mean_outputs/exp_{experiment}_mean_dists_{coi}_cluster{cluster}.png" #"results/mean_outputs/exp_{experiment}_mean_dists_{coi}.png"
    output: "results/plots/mean_dists_{experiment}_{coi}_cluster{cluster}.png"
    resources: mem_mb = 10000
    conda: "../envs/plot.yaml"
    script: "../scripts/plot_mean_dist_coi.py"

rule plot_pw_mean_strand_clusterdist_binned:
    input: "results/mean_outputs_binned/exp_{experiment}_mean_dists_{coi}_cluster{cluster}.png" #"results/mean_outputs/exp_{experiment}_mean_dists_{coi}.png"
    output: "results/plots/mean_dists_binned_{experiment}_{coi}_cluster{cluster}.png"
    resources: mem_mb = 10000
    conda: "../envs/plot.yaml"
    script: "../scripts/plot_mean_dist_coi.py"

rule plot_all_pw_mean_strand_clusterdist_binned:
    input: "results/mean_outputs_binned/exp_{experiment}_mean_dists_all_chrm_cluster{cluster}.png" #"results/mean_outputs/exp_{experiment}_mean_dists_{coi}.png"
    output: "results/plots/mean_dists_binned_{experiment}_all_chrm_cluster{cluster}.png"
    resources: mem_mb = 30000
    conda: "../envs/plot.yaml"
    script: "../scripts/plot_mean_dist_coi.py"

rule plot_pw_mean_strand_clusterdist_paint_binned:
    input: "results/mean_outputs_paint_binned/exp_{experiment}_mean_dists_{coi}_cluster{cluster}.png" #"results/mean_outputs/exp_{experiment}_mean_dists_{coi}.png"
    output: "results/plots/mean_dists_paint_binned_{experiment}_{coi}_cluster{cluster}.png"
    resources: mem_mb = 10000
    conda: "../envs/plot.yaml"
    script: "../scripts/plot_mean_dist_coi.py"

rule plot_all_pw_mean_strand_clusterdist_paint_binned:
    input: "results/mean_outputs_paint_binned/exp_{experiment}_mean_dists_all_chrm_cluster{cluster}.png" #"results/mean_outputs/exp_{experiment}_mean_dists_{coi}.png"
    output: "results/plots/mean_dists_paint_binned_{experiment}_all_chrm_cluster{cluster}.png"
    resources: mem_mb = 30000
    conda: "../envs/plot.yaml"
    script: "../scripts/plot_mean_dist_coi.py"


rule plot_pw_median_dist:
    input: "results/median_outputs/exp_{experiment}_median_dists_{coi}.png"
    output: "results/plots/median_dists_{experiment}_{coi}.png"
    resources: mem_mb = 1000
    conda: "../envs/plot.yaml"
    script: "../scripts/plot_mean_dist_coi.py"

    
rule plot_pw_all_mean_dist_binned:
    input: "results/mean_outputs_binned/exp_{experiment}_mean_dists_all_chrm.png"
    output: "results/plots/mean_dists_binned_{experiment}_all_chrm.png"
    resources: mem_mb = 30000
    conda: "../envs/plot.yaml"
    script: "../scripts/plot_mean_dist_all.py"

rule plot_pw_all_mean_dist_paint_binned:
    input: "results/mean_outputs_paint_binned/exp_{experiment}_mean_dists_all_chrm.png"
    output: "results/plots/mean_dists_paint_binned_{experiment}_all_chrm.png"
    resources: mem_mb = 30000
    conda: "../envs/plot.yaml"
    script: "../scripts/plot_mean_dist_all.py"