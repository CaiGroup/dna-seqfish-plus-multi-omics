
#brain single reps
for rep in [1,2]:
    # get_pw_mean_dists_Brain_single rep
    rule:
        input:
            "resources/cerebellum_nuc_vol_filtered_rep1+2_stats.csv",
            expand("resources/cerebellum/Brain_rep%d_pos_{pos}_radial_scores_wofilter.csv" % rep, pos=config["experiment_positions"]["Brain_rep%d" % rep]),
        output: 
            "results/mean_outputs/exp_Brain_rep%d_mean_dists_{coi}_cluster2.png" % rep,
            "results/mean_outputs/exp_Brain_rep%d_mean_dists_{coi}_cluster4.png" % rep,
            "results/mean_outputs/exp_Brain_rep%d_mean_dists_{coi}_cluster7.png" % rep,
            "results/mean_outputs/exp_Brain_rep%d_mean_dists_{coi}_cluster6&11.png" % rep
        resources: mem_mb = 20000
        script: "../scripts/calc_mean_dists_brain_clust_strands_coi.jl"

    #get_pw_mean_dists_Brain_single rep binned:
    rule:
        input:
            "resources/cerebellum_nuc_vol_filtered_rep1+2_stats.csv",
            expand("results/cerebellum_binned/Brain_rep%d_pos_{pos}_radial_scores_wofilter.csv" % rep, pos=config["experiment_positions"]["Brain_rep%d" % rep]),
        output: 
            "results/mean_outputs_binned/exp_Brain_rep%d_mean_dists_{coi}_cluster2.png" % rep,
            "results/mean_outputs_binned/exp_Brain_rep%d_mean_dists_{coi}_cluster4.png" % rep,
            "results/mean_outputs_binned/exp_Brain_rep%d_mean_dists_{coi}_cluster7.png" % rep,
            "results/mean_outputs_binned/exp_Brain_rep%d_mean_dists_{coi}_cluster6&11.png" % rep
    	resources: mem_mb = 10000
    	script: "../scripts/calc_mean_dists_brain_clust_strands_coi.jl"

    rule:
        input:
            "resources/cerebellum_nuc_vol_filtered_rep1+2_stats.csv",
            expand("results/cerebellum_paint_binned/Brain_rep%d_pos_{pos}_radial_scores_wofilter.csv" % rep, pos=config["experiment_positions"]["Brain_rep%d" % rep]),
        output: 
            "results/mean_outputs_paint_binned/exp_Brain_rep%d_mean_dists_{coi}_cluster2.png" % rep,
            "results/mean_outputs_paint_binned/exp_Brain_rep%d_mean_dists_{coi}_cluster4.png" % rep,
            "results/mean_outputs_paint_binned/exp_Brain_rep%d_mean_dists_{coi}_cluster7.png" % rep,
            "results/mean_outputs_paint_binned/exp_Brain_rep%d_mean_dists_{coi}_cluster6&11.png" % rep
    	resources: mem_mb = 10000
    	script: "../scripts/calc_mean_dists_brain_clust_strands_coi.jl"

    # get all pw meand dists brain single rep
    rule:
        input:
            "resources/cerebellum_nuc_vol_filtered_rep1+2_stats.csv",
            expand("results/cerebellum_binned/Brain_rep%d_pos_{pos}_radial_scores_wofilter.csv" % rep, pos=config["experiment_positions"]["Brain_rep%d" % rep]),
        output: 
            "results/mean_outputs_binned/exp_Brain_rep%d_mean_dists_all_chrm_cluster2.png" % rep,
            "results/mean_outputs_binned/exp_Brain_rep%d_mean_dists_all_chrm_cluster4.png" % rep,
            "results/mean_outputs_binned/exp_Brain_rep%d_mean_dists_all_chrm_cluster7.png" % rep,
            "results/mean_outputs_binned/exp_Brain_rep%d_mean_dists_all_chrm_cluster6&11.png" % rep
    	threads: 20 
        resources: mem_mb = 3000
    	script: "../scripts/calc_mean_dists_brain_clust_all_chrms.jl"

    rule:
        input:
            "resources/cerebellum_nuc_vol_filtered_rep1+2_stats.csv",
            expand("results/cerebellum_paint_binned/Brain_rep%d_pos_{pos}_radial_scores_wofilter.csv" % rep, pos=config["experiment_positions"]["Brain_rep%d" % rep]),
        output: 
            "results/mean_outputs_paint_binned/exp_Brain_rep%d_mean_dists_all_chrm_cluster2.png" % rep,
            "results/mean_outputs_paint_binned/exp_Brain_rep%d_mean_dists_all_chrm_cluster4.png" % rep,
            "results/mean_outputs_paint_binned/exp_Brain_rep%d_mean_dists_all_chrm_cluster7.png" % rep,
            "results/mean_outputs_paint_binned/exp_Brain_rep%d_mean_dists_all_chrm_cluster6&11.png" % rep
    	threads: 20 
        resources: mem_mb = 3000
    	script: "../scripts/calc_mean_dists_brain_clust_all_chrms.jl"


        # get_pw_proximities_Brain_single rep
    rule:
        input:
            "resources/cerebellum_nuc_vol_filtered_rep1+2_stats.csv",
            expand("resources/cerebellum/Brain_rep%d_pos_{pos}_radial_scores_wofilter.csv" % rep, pos=config["experiment_positions"]["Brain_rep%d" % rep]),
        output: 
            "results/prox_freq_outputs/exp_Brain_rep%d_prox_counts_{coi}_r{r}_cluster2.png" % rep,
            "results/prox_freq_outputs/exp_Brain_rep%d_prox_counts_{coi}_r{r}_cluster4.png" % rep,
            "results/prox_freq_outputs/exp_Brain_rep%d_prox_counts_{coi}_r{r}_cluster7.png" % rep,
            "results/prox_freq_outputs/exp_Brain_rep%d_prox_counts_{coi}_r{r}_cluster6&11.png" % rep,
            "results/prox_freq_outputs/exp_Brain_rep%d_ncoocur_{coi}_r{r}_cluster2.png" % rep,
            "results/prox_freq_outputs/exp_Brain_rep%d_ncoocur_{coi}_r{r}_cluster4.png" % rep,
            "results/prox_freq_outputs/exp_Brain_rep%d_ncoocur_{coi}_r{r}_cluster7.png" % rep,
            "results/prox_freq_outputs/exp_Brain_rep%d_ncoocur_{coi}_r{r}_cluster6&11.png" % rep,
        resources: mem_mb = 20000
        script: "../scripts/count_prox_freq_brain_clustered_strands_coi_var_r.jl"

    #get_pw_proximities_Brain_single rep binned:
    rule:
        input:
            "resources/cerebellum_nuc_vol_filtered_rep1+2_stats.csv",
            expand("results/cerebellum_binned/Brain_rep%d_pos_{pos}_radial_scores_wofilter.csv" % rep, pos=config["experiment_positions"]["Brain_rep%d" % rep]),
        output: 
            "results/prox_freq_outputs_binned/exp_Brain_rep%d_prox_counts_{coi}_r{r}_cluster2.png" % rep,
            "results/prox_freq_outputs_binned/exp_Brain_rep%d_prox_counts_{coi}_r{r}_cluster4.png" % rep,
            "results/prox_freq_outputs_binned/exp_Brain_rep%d_prox_counts_{coi}_r{r}_cluster7.png" % rep,
            "results/prox_freq_outputs_binned/exp_Brain_rep%d_prox_counts_{coi}_r{r}_cluster6&11.png" % rep,
            "results/prox_freq_outputs_binned/exp_Brain_rep%d_ncoocur_{coi}_r{r}_cluster2.png" % rep,
            "results/prox_freq_outputs_binned/exp_Brain_rep%d_ncoocur_{coi}_r{r}_cluster4.png" % rep,
            "results/prox_freq_outputs_binned/exp_Brain_rep%d_ncoocur_{coi}_r{r}_cluster7.png" % rep,
            "results/prox_freq_outputs_binned/exp_Brain_rep%d_ncoocur_{coi}_r{r}_cluster6&11.png" % rep
    	resources: mem_mb = 10000
    	script: "../scripts/count_prox_freq_brain_clustered_strands_coi_var_r.jl"

    rule:
        input:
            "resources/cerebellum_nuc_vol_filtered_rep1+2_stats.csv",
            expand("results/cerebellum_paint_binned/Brain_rep%d_pos_{pos}_radial_scores_wofilter.csv" % rep, pos=config["experiment_positions"]["Brain_rep%d" % rep]),
        output: 
            "results/prox_freq_outputs_paint_binned/exp_Brain_rep%d_prox_counts_{coi}_r{r}_cluster2.png" % rep,
            "results/prox_freq_outputs_paint_binned/exp_Brain_rep%d_prox_counts_{coi}_r{r}_cluster4.png" % rep,
            "results/prox_freq_outputs_paint_binned/exp_Brain_rep%d_prox_counts_{coi}_r{r}_cluster7.png" % rep,
            "results/prox_freq_outputs_paint_binned/exp_Brain_rep%d_prox_counts_{coi}_r{r}_cluster6&11.png" % rep,
            "results/prox_freq_outputs_paint_binned/exp_Brain_rep%d_ncoocur_{coi}_r{r}_cluster2.png" % rep,
            "results/prox_freq_outputs_paint_binned/exp_Brain_rep%d_ncoocur_{coi}_r{r}_cluster4.png" % rep,
            "results/prox_freq_outputs_paint_binned/exp_Brain_rep%d_ncoocur_{coi}_r{r}_cluster7.png" % rep,
            "results/prox_freq_outputs_paint_binned/exp_Brain_rep%d_ncoocur_{coi}_r{r}_cluster6&11.png" % rep
    	resources: mem_mb = 10000
    	script: "../scripts/count_prox_freq_brain_clustered_strands_coi_var_r.jl"