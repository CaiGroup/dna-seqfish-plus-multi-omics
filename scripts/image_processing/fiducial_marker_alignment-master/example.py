import pandas as pd
import numpy as np
from copy import deepcopy

from FMAligner import FMAligner

# Datafiles are not in repository because of space
ref_final = pd.read_csv("reference2-points-pos0.csv")
ref_init = pd.read_csv("reference-points-pos0.csv")
ro = pd.read_csv("hybridization-points-pos0.csv")

channels = [1, 2, 3]

all_channel_offsets = []
for channel in channels:
    ch_ref = ref_init.set_index("ch").loc[channel]
    ch_ro = ro.set_index("ch").loc[channel].set_index("hyb")
    ch_ref_final = ref_final.set_index("ch").loc[channel]
    
    #ch_ro = ch_ro.loc[96]
    

    dal = FMAligner(ch_ro, ch_ref, ch_ref_final)
    #dal.load_saved_parameters("ch2_params.txt")
    dal.auto_set_params()
    dal.set_max_lat_offset(20)
    
    dal.set_max_z_offset(8)
    dal2 = deepcopy(dal)
    

    dal.align()
    
    #find if any hybridizations failed to align

    hnf = dal.offsets.loc[np.isnan(dal.offsets['x'])]
    if len(hnf) > 0:
        dal2.set_n_longest_edges(dal.n_longest_edges*30)
        dal2.min_dot_matches = 4
        dal2.min_edge_matches = 2
        #dal2.align()
        for hyb in hnf['hyb']:
            dal.offsets.loc[hyb-1, 'x':] = dal2.align_hyb(hyb)

    
    dal.save_offsets("offsets_ch" + str(channel) + ".csv")
    #dal.save_offsets_SEs("offsets_ch_SEs" + str(channel) + ".csv")
    dal.save_matches("matches_ch" + str(channel) + ".csv")
    dal.save_loocv_errors("loov_errors_ch" + str(channel) + ".csv")
    dal.save_ro_wout_fm("pnts_no_fm_ch" + str(channel) + ".csv")
    
    offsets = dal.offsets
    dal = None
    dal2 = None
    
    offsets.insert(0, column="ch", value=channel)
    all_channel_offsets.append(offsets)
    

agg = pd.concat(all_channel_offsets)
agg.to_csv("all_ch_offsets.csv")
    