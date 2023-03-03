rule get_pw_mean_dists_Brain_rep1_rep2:
    input:
        "resources/cerebellum_nuc_vol_filtered_rep1+2_stats.csv",
        expand("resources/cerebellum/Brain_rep1_pos_{pos}_radial_scores_wofilter.csv", pos=config["experiment_positions"]["Brain_rep1"]),
        expand("resources/cerebellum/Brain_rep2_pos_{pos}_radial_scores_wofilter.csv", pos=config["experiment_positions"]["Brain_rep2"]),
    output: 
        "results/mean_outputs/exp_Brain_rep1&2_mean_dists_{coi}_cluster2.png",
        "results/mean_outputs/exp_Brain_rep1&2_mean_dists_{coi}_cluster4.png",
        "results/mean_outputs/exp_Brain_rep1&2_mean_dists_{coi}_cluster7.png",
        "results/mean_outputs/exp_Brain_rep1&2_mean_dists_{coi}_cluster6&11.png"
    resources: mem_mb = 30000
    script: "../scripts/calc_mean_dists_brain_clust_strands_coi.jl"

rule get_pw_proximities_Brain_rep1_rep2:
    input:
        "resources/cerebellum_nuc_vol_filtered_rep1+2_stats.csv",
        expand("resources/cerebellum/Brain_rep1_pos_{pos}_radial_scores_wofilter.csv", pos=config["experiment_positions"]["Brain_rep1"]),
        expand("resources/cerebellum/Brain_rep2_pos_{pos}_radial_scores_wofilter.csv", pos=config["experiment_positions"]["Brain_rep2"]),
    output: 
        "results/prox_freq_outputs/exp_Brain_rep1&2_prox_counts_{coi}_r{r}_cluster2.png",
        "results/prox_freq_outputs/exp_Brain_rep1&2_prox_counts_{coi}_r{r}_cluster4.png",
        "results/prox_freq_outputs/exp_Brain_rep1&2_prox_counts_{coi}_r{r}_cluster7.png",
        "results/prox_freq_outputs/exp_Brain_rep1&2_prox_counts_{coi}_r{r}_cluster6&11.png",
        "results/prox_freq_outputs/exp_Brain_rep1&2_ncoocur_{coi}_r{r}_cluster2.png",
        "results/prox_freq_outputs/exp_Brain_rep1&2_ncoocur_{coi}_r{r}_cluster4.png",
        "results/prox_freq_outputs/exp_Brain_rep1&2_ncoocur_{coi}_r{r}_cluster7.png",
        "results/prox_freq_outputs/exp_Brain_rep1&2_ncoocur_{coi}_r{r}_cluster6&11.png",
    resources: mem_mb = 30000
    script: "../scripts/count_prox_freq_brain_clustered_strands_coi_var_r.jl"



rule get_pw_mean_dists_Brain_rep1_rep2_binned:
    input:
        "resources/cerebellum_nuc_vol_filtered_rep1+2_stats.csv",
        expand("results/cerebellum_binned/Brain_rep1_pos_{pos}_radial_scores_wofilter.csv", pos=config["experiment_positions"]["Brain_rep1"]),
        expand("results/cerebellum_binned/Brain_rep2_pos_{pos}_radial_scores_wofilter.csv", pos=config["experiment_positions"]["Brain_rep2"]),
    output: 
        "results/mean_outputs_binned/exp_Brain_rep1&2_mean_dists_{coi}_cluster2.png",
        "results/mean_outputs_binned/exp_Brain_rep1&2_mean_dists_{coi}_cluster4.png",
        "results/mean_outputs_binned/exp_Brain_rep1&2_mean_dists_{coi}_cluster7.png",
        "results/mean_outputs_binned/exp_Brain_rep1&2_mean_dists_{coi}_cluster6&11.png",
    resources: mem_mb = 30000
    script: "../scripts/calc_mean_dists_brain_clust_strands_coi.jl"

rule get_pw_mean_dists_Brain_rep1_rep2_paint_binned:
    input:
        "resources/cerebellum_nuc_vol_filtered_rep1+2_stats.csv",
        expand("results/cerebellum_paint_binned/Brain_rep1_pos_{pos}_radial_scores_wofilter.csv", pos=config["experiment_positions"]["Brain_rep1"]),
        expand("results/cerebellum_paint_binned/Brain_rep2_pos_{pos}_radial_scores_wofilter.csv", pos=config["experiment_positions"]["Brain_rep2"]),
    output: 
        "results/mean_outputs_paint_binned/exp_Brain_rep1&2_mean_dists_{coi}_cluster2.png",
        "results/mean_outputs_paint_binned/exp_Brain_rep1&2_mean_dists_{coi}_cluster4.png",
        "results/mean_outputs_paint_binned/exp_Brain_rep1&2_mean_dists_{coi}_cluster7.png",
        "results/mean_outputs_paint_binned/exp_Brain_rep1&2_mean_dists_{coi}_cluster6&11.png",
    resources: mem_mb = 30000
    script: "../scripts/calc_mean_dists_brain_clust_strands_coi.jl"


rule get_all_pw_mean_dists_Brain_rep1_rep2_binned:
    input:
        "resources/cerebellum_nuc_vol_filtered_rep1+2_stats.csv",
        expand("results/cerebellum_binned/Brain_rep1_pos_{pos}_radial_scores_wofilter.csv", pos=config["experiment_positions"]["Brain_rep1"]),
        expand("results/cerebellum_binned/Brain_rep2_pos_{pos}_radial_scores_wofilter.csv", pos=config["experiment_positions"]["Brain_rep2"]),
    output: 
        "results/mean_outputs_binned/exp_Brain_rep1&2_mean_dists_all_chrm_cluster2.png",
        "results/mean_outputs_binned/exp_Brain_rep1&2_mean_dists_all_chrm_cluster4.png",
        "results/mean_outputs_binned/exp_Brain_rep1&2_mean_dists_all_chrm_cluster7.png",
        "results/mean_outputs_binned/exp_Brain_rep1&2_mean_dists_all_chrm_cluster6&11.png"
    threads: 20
    resources: mem_mb = 4000
    script: "../scripts/calc_mean_dists_brain_clust_all_chrms.jl"

rule get_all_pw_mean_dists_Brain_rep1_rep2_paint_binned:
    input:
        "resources/cerebellum_nuc_vol_filtered_rep1+2_stats.csv",
        expand("results/cerebellum_paint_binned/Brain_rep1_pos_{pos}_radial_scores_wofilter.csv", pos=config["experiment_positions"]["Brain_rep1"]),
        expand("results/cerebellum_paint_binned/Brain_rep2_pos_{pos}_radial_scores_wofilter.csv", pos=config["experiment_positions"]["Brain_rep2"]),
    output: 
        "results/mean_outputs_paint_binned/exp_Brain_rep1&2_mean_dists_all_chrm_cluster2.png",
        "results/mean_outputs_paint_binned/exp_Brain_rep1&2_mean_dists_all_chrm_cluster4.png",
        "results/mean_outputs_paint_binned/exp_Brain_rep1&2_mean_dists_all_chrm_cluster7.png",
        "results/mean_outputs_paint_binned/exp_Brain_rep1&2_mean_dists_all_chrm_cluster6&11.png"
    threads: 20
    resources: mem_mb = 4000
    script: "../scripts/calc_mean_dists_brain_clust_all_chrms.jl"


rule get_all_pw_prox_reqs_brain_rep1_rep2:
        input: 
            "resources/cerebellum_nuc_vol_filtered_rep1+2_stats.csv",
            expand("results/cerebellum_binned/Brain_rep1_pos_{pos}_filtered_renumbered.csv", pos=config["experiment_positions"]["Brain_rep1"]),
            expand("results/cerebellum_binned/Brain_rep2_pos_{pos}_filtered_renumbered.csv", pos=config["experiment_positions"]["Brain_rep2"])
        output: 
            "results/prox_freq_outputs_binned/exp_Brain_rep1&2_prox_counts_all_r{r}.png",
            "results/prox_freq_outputs/exp_Brain_rep1&2_ncoocur_all_r{r}.png"
        resources: mem_mb = 3000
        threads: 20
        script: "../scripts/count_prox_freq_all.jl"

rule get_pw_proximities_Brain_rep1_rep2_binned:
    input:
        "resources/cerebellum_nuc_vol_filtered_rep1+2_stats.csv",
        expand("results/cerebellum_binned/Brain_rep1_pos_{pos}_radial_scores_wofilter.csv", pos=config["experiment_positions"]["Brain_rep1"]),
        expand("results/cerebellum_binned/Brain_rep2_pos_{pos}_radial_scores_wofilter.csv", pos=config["experiment_positions"]["Brain_rep2"]),
    output: 
        "results/prox_freq_outputs_binned/exp_Brain_rep1&2_prox_counts_{coi}_r{r}_cluster2.png",
        "results/prox_freq_outputs_binned/exp_Brain_rep1&2_prox_counts_{coi}_r{r}_cluster4.png",
        "results/prox_freq_outputs_binned/exp_Brain_rep1&2_prox_counts_{coi}_r{r}_cluster7.png",
        "results/prox_freq_outputs_binned/exp_Brain_rep1&2_prox_counts_{coi}_r{r}_cluster6&11.png",
        "results/prox_freq_outputs_binned/exp_Brain_rep1&2_ncoocur_{coi}_r{r}_cluster2.png",
        "results/prox_freq_outputs_binned/exp_Brain_rep1&2_ncoocur_{coi}_r{r}_cluster4.png",
        "results/prox_freq_outputs_binned/exp_Brain_rep1&2_ncoocur_{coi}_r{r}_cluster7.png",
        "results/prox_freq_outputs_binned/exp_Brain_rep1&2_ncoocur_{coi}_r{r}_cluster6&11.png"
    resources: mem_mb = 30000
    script: "../scripts/count_prox_freq_brain_clustered_strands_coi_var_r.jl"

rule get_pw_proximities_Brain_rep1_rep2_paint_binned:
    input:
        "resources/cerebellum_nuc_vol_filtered_rep1+2_stats.csv",
        expand("results/cerebellum_paint_binned/Brain_rep1_pos_{pos}_radial_scores_wofilter.csv", pos=config["experiment_positions"]["Brain_rep1"]),
        expand("results/cerebellum_paint_binned/Brain_rep2_pos_{pos}_radial_scores_wofilter.csv", pos=config["experiment_positions"]["Brain_rep2"]),
    output: 
        "results/prox_freq_outputs_paint_binned/exp_Brain_rep1&2_prox_counts_{coi}_r{r}_cluster2.png",
        "results/prox_freq_outputs_paint_binned/exp_Brain_rep1&2_prox_counts_{coi}_r{r}_cluster4.png",
        "results/prox_freq_outputs_paint_binned/exp_Brain_rep1&2_prox_counts_{coi}_r{r}_cluster7.png",
        "results/prox_freq_outputs_paint_binned/exp_Brain_rep1&2_prox_counts_{coi}_r{r}_cluster6&11.png",
        "results/prox_freq_outputs_paint_binned/exp_Brain_rep1&2_ncoocur_{coi}_r{r}_cluster2.png",
        "results/prox_freq_outputs_paint_binned/exp_Brain_rep1&2_ncoocur_{coi}_r{r}_cluster4.png",
        "results/prox_freq_outputs_paint_binned/exp_Brain_rep1&2_ncoocur_{coi}_r{r}_cluster7.png",
        "results/prox_freq_outputs_paint_binned/exp_Brain_rep1&2_ncoocur_{coi}_r{r}_cluster6&11.png"
    resources: mem_mb = 30000
    script: "../scripts/count_prox_freq_brain_clustered_strands_coi_var_r.jl"

