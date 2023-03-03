for experiment in config["experiments"]:       
    # get pairwise mean distances 
    rule:
        input: expand("results/filtered_renumbered_dots_binned/"+experiment+"_pos_{pos}_filtered_renumbered.csv", pos=config["experiment_positions"][experiment])
        output: "results/mean_outputs_binned/exp_"+experiment+"_mean_dists_all_chrm.png"
        resources: mem_mb = 4000
        threads: 20 
        script: "../scripts/calc_mean_dists_missing_0xffff.jl"

    rule:
        input: expand("results/filtered_renumbered_dots_paint_binned/"+experiment+"_pos_{pos}_filtered_renumbered.csv", pos=config["experiment_positions"][experiment])
        output: "results/mean_outputs_paint_binned/exp_"+experiment+"_mean_dists_all_chrm.png"
        resources: mem_mb = 4000
        threads: 20 
        script: "../scripts/calc_mean_dists_missing_0xffff.jl"

    # get pairwise proximity frequencies
    rule:
        input: expand("results/filtered_renumbered_dots_binned/"+experiment+"_pos_{pos}_filtered_renumbered.csv", pos=config["experiment_positions"][experiment])
        output: "results/prox_freq_outputs_binned/exp_"+experiment+"_prox_counts_all.png"
        resources: mem_mb = 3000
        threads: 20
        script: "../scripts/count_prox_freq_all.jl"

    rule:
        input: expand("results/filtered_renumbered_dots_binned/"+experiment+"_pos_{pos}_filtered_renumbered.csv", pos=config["experiment_positions"][experiment])
        output: "results/prox_freq_outputs_binned/exp_"+experiment+"_ncoocur_all.png"
        resources: mem_mb = 10000
        threads: 15
        script: "../scripts/count_cooccurences_all.jl"

    rule:
        input: expand("results/filtered_renumbered_dots_paint_binned/"+experiment+"_pos_{pos}_filtered_renumbered.csv", pos=config["experiment_positions"][experiment])
        output: "results/prox_freq_outputs_paint_binned/exp_"+experiment+"_prox_counts_all.png"
        resources: mem_mb = 3000
        threads: 20
        script: "../scripts/count_prox_freq_all.jl"

    rule:
        input: expand("results/filtered_renumbered_dots_paint_binned/"+experiment+"_pos_{pos}_filtered_renumbered.csv", pos=config["experiment_positions"][experiment])
        output: "results/prox_freq_outputs_paint_binned/exp_"+experiment+"_ncoocur_all.png"
        resources: mem_mb = 10000
        threads: 10
        script: "../scripts/count_cooccurences_all.jl"

    # Process Chrommosomes of interest
    #############################

    # get chromosome of interest pairwise mean distances
    rule:
        input: expand("resources/filtered_renumbered_dots/"+experiment+"_pos_{pos}_filtered_renumbered.csv", pos=config["experiment_positions"][experiment])
        output: "results/mean_outputs/exp_"+experiment+"_mean_dists_{coi}.png"#, "results/mean_outputs/"+experiment+"_nsums_{coi}.png"
        resources: mem_mb = 15000
        script: "../scripts/calc_mean_dists_strands_coi.jl"

    rule:
        input: expand("results/filtered_renumbered_dots_binned/"+experiment+"_pos_{pos}_filtered_renumbered.csv", pos=config["experiment_positions"][experiment])
        output: "results/mean_outputs_binned/exp_"+experiment+"_mean_dists_{coi}.png"#, "results/mean_outputs_binned/"+experiment+"_nsums_{coi}.png"
        resources: mem_mb = 15000
        script: "../scripts/calc_mean_dists_strands_coi.jl"

    rule:
        input: expand("results/filtered_renumbered_dots_paint_binned/"+experiment+"_pos_{pos}_filtered_renumbered.csv", pos=config["experiment_positions"][experiment])
        output: "results/mean_outputs_paint_binned/exp_"+experiment+"_mean_dists_{coi}.png"
        resources: mem_mb = 15000
        script: "../scripts/calc_mean_dists_strands_coi.jl"

    # get chromosome of interest pairwise median distances
    rule:
        input: expand("resources/filtered_renumbered_dots/"+experiment+"_pos_{pos}_filtered_renumbered.csv", pos=config["experiment_positions"][experiment])
        output: "results/median_outputs/exp_"+experiment+"_median_dists_{coi}.png" #, "results/median_outputs/"+experiment+"_nsums_{coi}.png"
        resources: mem_mb = 10000
        script: "../scripts/calc_median_dists_strands_coi.jl"

    # get chromosome of interest pairwise proximity frequencies
    rule:
        input: expand("resources/filtered_renumbered_dots/"+experiment+"_pos_{pos}_filtered_renumbered.csv", pos=config["experiment_positions"][experiment])
        output: "results/prox_freq_outputs/"+experiment+"_prox_counts_{coi}.png", "results/prox_freq_outputs/"+experiment+"_ncoocur_{coi}.png"
        resources: mem_mb = 10000
        script: "../scripts/count_prox_freq_strands_coi.jl"

     # get chromosome of interest pairwise proximity frequencies binned
    rule:
        input: expand("results/filtered_renumbered_dots_binned/"+experiment+"_pos_{pos}_filtered_renumbered.csv", pos=config["experiment_positions"][experiment])
        output: "results/prox_freq_outputs_binned/"+experiment+"_prox_counts_{coi}.png", "results/prox_freq_outputs_binned/"+experiment+"_ncoocur_{coi}.png"
        resources: mem_mb = 10000
        script: "../scripts/count_prox_freq_strands_coi.jl"

    rule:
        input: expand("results/filtered_renumbered_dots_binned/"+experiment+"_pos_{pos}_filtered_renumbered.csv", pos=config["experiment_positions"][experiment])
        output: "results/prox_freq_outputs_binned_binned/"+experiment+"_prox_counts_{coi}.png", "results/prox_freq_outputs_binned/"+experiment+"_ncoocur_{coi}.png"
        resources: mem_mb = 10000
        script: "../scripts/count_prox_freq_strands_coi.jl"

    # for expeiments with clustering, find separate pairwise distance/proximity maps for each cluster
    if experiment in config["clustered_experiments"]:

        #all pairwise mean distances
        rule:
            input: 
                expand("results/filtered_renumbered_dots_binned/"+experiment+"_pos_{pos}_filtered_renumbered.csv", pos=config["experiment_positions"][experiment]),
                "resources/100k-001-027-E14mDuxCA-mRNA-clustering-0430.csv"
            params: experiment = experiment
            output: "results/mean_outputs_binned/exp_"+experiment+"_mean_dists_all_chrm_cluster{cluster}.png" #, "results/mean_outputs/"+experiment+"_nsums_all_chrm.png"
            resources: mem_mb = 4000
            threads: 20 
            script: "../scripts/calc_mean_dists_missing_0xffff_cluster.jl"

        rule:
            input: 
                expand("results/filtered_renumbered_dots_paint_binned/"+experiment+"_pos_{pos}_filtered_renumbered.csv", pos=config["experiment_positions"][experiment]),
                "resources/100k-001-027-E14mDuxCA-mRNA-clustering-0430.csv"
            params: experiment = experiment
            output: "results/mean_outputs_paint_binned/exp_"+experiment+"_mean_dists_all_chrm_cluster{cluster}.png" #, "results/mean_outputs/"+experiment+"_nsums_all_chrm.png"
            resources: mem_mb = 4000
            threads: 20 
            script: "../scripts/calc_mean_dists_missing_0xffff_cluster.jl"

        # plot single chromosome mean pw dists
        rule:
            input:
                expand("resources/filtered_renumbered_dots/"+experiment+"_pos_{pos}_filtered_renumbered.csv", pos=config["experiment_positions"][experiment]),
                "resources/100k-001-027-E14mDuxCA-mRNA-clustering-0430.csv"
            params: experiment = experiment
	    output: "results/mean_outputs/exp_"+experiment+"_mean_dists_{coi}_cluster{cluster}.png"
            resources: mem_mb = 20000
            script: "../scripts/calc_mean_dists_clustered_strands_coi.jl"

        rule:
            input:
                expand("resources/filtered_renumbered_dots/"+experiment+"_pos_{pos}_filtered_renumbered.csv", pos=config["experiment_positions"][experiment]),
                "resources/100k-001-027-E14mDuxCA-mRNA-clustering-0430.csv"
            params: experiment = experiment
            output: "results/prox_freq_outputs/"+experiment+"_prox_counts_{coi}_cluster{cluster}.png", "results/prox_freq_outputs/"+experiment+"_ncoocur_{coi}_cluster{cluster}.png"
            resources: mem_mb = 10000
            script: "../scripts/count_prox_freq_clustered strands_coi.jl"


rule get_E14_rep12_paint_binned_all_chrm_mean_dists_clustered:
    input: 
        expand("results/filtered_renumbered_dots_paint_binned/E14_rep2_replicate1_pos_{pos}_filtered_renumbered.csv", pos=config["experiment_positions"]["E14_rep2_replicate1"]),
        expand("results/filtered_renumbered_dots_paint_binned/E14_rep3_replicate2_pos_{pos}_filtered_renumbered.csv", pos=config["experiment_positions"]["E14_rep3_replicate2"]),
        "resources/100k-001-027-E14mDuxCA-mRNA-clustering-0430.csv"
    params: experiment = experiment
    output: "results/mean_outputs_paint_binned/exp_E14_rep1&2_mean_dists_all_chrm_cluster{cluster}.png" #, "results/mean_outputs/"+experiment+"_nsums_all_chrm.png"
    resources: mem_mb = 4000
    threads: 25 
    script: "../scripts/calc_mean_dists_missing_0xffff_all_reps_all_chrms.jl"

rule get_14_rep12_paint_binned_all_chrm_mean_dists:
    input: 
        expand("results/filtered_renumbered_dots_paint_binned/E14_rep2_replicate1_pos_{pos}_filtered_renumbered.csv", pos=config["experiment_positions"]["E14_rep2_replicate1"]),
        expand("results/filtered_renumbered_dots_paint_binned/E14_rep3_replicate2_pos_{pos}_filtered_renumbered.csv", pos=config["experiment_positions"]["E14_rep3_replicate2"]),
    output: "results/mean_outputs_paint_binned/exp_E14_rep1&2_mean_dists_all_chrm.png"
    resources: mem_mb = 7000
    threads: 20 
    script: "../scripts/calc_mean_dists_missing_0xffff_all_reps.jl"

rule get_14_rep12_binned_all_chrm_mean_dists:
    input: 
        expand("results/filtered_renumbered_dots_binned/E14_rep2_replicate1_pos_{pos}_filtered_renumbered.csv", pos=config["experiment_positions"]["E14_rep2_replicate1"]),
        expand("results/filtered_renumbered_dots_binned/E14_rep3_replicate2_pos_{pos}_filtered_renumbered.csv", pos=config["experiment_positions"]["E14_rep3_replicate2"]),
    output: "results/mean_outputs_binned/exp_E14_rep1&2_mean_dists_all_chrm.png"
    resources: mem_mb = 7000
    threads: 20 
    script: "../scripts/calc_mean_dists_missing_0xffff_all_reps.jl"

rule get_E14_rep12_prox_freq_coi:
    input: 
        expand("resources/filtered_renumbered_dots/E14_rep2_replicate1_pos_{pos}_filtered_renumbered.csv", pos=config["experiment_positions"]["E14_rep2_replicate1"]),
        expand("resources/filtered_renumbered_dots/E14_rep3_replicate2_pos_{pos}_filtered_renumbered.csv", pos=config["experiment_positions"]["E14_rep3_replicate2"]),
    output: "results/prox_freq_outputs/E14_rep1&2_prox_counts_{coi}.png", "results/prox_freq_outputs/E14_rep1&2_ncoocur_{coi}.png"
    resources: mem_mb = 30000
    script: "../scripts/count_prox_freq_strands_coi.jl"

rule get_E14_rep12_prox_freq_all:
    input: 
        expand("results/filtered_renumbered_dots_paint_binned/E14_rep2_replicate1_pos_{pos}_filtered_renumbered.csv", pos=config["experiment_positions"]["E14_rep2_replicate1"]),
        expand("results/filtered_renumbered_dots_paint_binned/E14_rep3_replicate2_pos_{pos}_filtered_renumbered.csv", pos=config["experiment_positions"]["E14_rep3_replicate2"]),
    output: "results/prox_freq_outputs_paint_binned/exp_E14_rep1&2_prox_counts_all.png"
    resources: mem_mb = 10000
    threads: 20
    script: "../scripts/count_prox_freq_all.jl"

rule get_E14_rep12_coocur_all:
    input: 
        expand("results/filtered_renumbered_dots_paint_binned/E14_rep2_replicate1_pos_{pos}_filtered_renumbered.csv", pos=config["experiment_positions"]["E14_rep2_replicate1"]),
        expand("results/filtered_renumbered_dots_paint_binned/E14_rep3_replicate2_pos_{pos}_filtered_renumbered.csv", pos=config["experiment_positions"]["E14_rep3_replicate2"]),
    output: "results/prox_freq_outputs_paint_binned/exp_E14_rep1&2_ncoocur_all.png"
    resources: mem_mb = 10000
    threads: 10
    script: "../scripts/count_cooccurences_all.jl"