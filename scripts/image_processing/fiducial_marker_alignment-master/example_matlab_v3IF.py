import pandas as pd
import os
import sys
import numpy as np
from copy import deepcopy
from datetime import datetime as dt

from FMAligner import FMAligner


quantile = 0

home_dir = sys.argv[1]
init_ref_fname = sys.argv[2]
final_ref_fname = sys.argv[3]
ro_fname = sys.argv[4]
savefname = sys.argv[5]
pos = sys.argv[6]


os.chdir(home_dir)

#ref_final = pd.read_csv(final_ref_fname)
ref_init = pd.read_csv(init_ref_fname)
ro = pd.read_csv(ro_fname)


datestr = dt.strftime(dt.now(), '%Y%m%d')
param_name = '.csv' 
savefname = savefname + datestr + param_name

channels = [3] # for DNA

all_channel_offsets = []
for channel in channels:
    ch_ref = ref_init.set_index("ch").loc[channel]
    ch_ro = ro.set_index("ch").loc[channel].set_index("hyb")
    #ch_ref_final = ref_final.set_index("ch").loc[channel]
    
    ch_ro = ch_ro.loc[ch_ro['int']>np.quantile(ch_ro['int'], quantile)]

    #dal = FMAligner(ch_ro, ch_ref, ch_ref_final)
    dal = FMAligner(ch_ro, ch_ref)
    #dal.load_saved_parameters("qp1_params_ch3_pos10.csv")
    #dal.auto_set_params()
    dal.set_xy_search_error(0.9)
    dal.set_z_search_error(0.5)
    dal.set_max_bright_prop(10)
    dal.set_min_bright_prop(0.4)
    dal.save_params("positions/qp1_params_ch" +str(channel) + "_pos" + str(pos) + ".csv")
    dal.set_max_lat_offset(100) # 100 default
    
    dal.set_max_z_offset(5) # 5 default
    dal.min_dot_matches = 5 # 4-10 for E14, 200-600 for brain, 50 for brain RNA FISH, manually switch.
    dal2 = deepcopy(dal)
    dal.set_n_longest_edges(1000) # 100 for brain, 200-400 for E14

    dal.align()
    
    #find if any hybridizations failed to align

    #hnf = dal.offsets.loc[np.isnan(dal.offsets['x'])]
    #if len(hnf) > 0:
        #dal2.set_n_longest_edges(dal.n_longest_edges*30)
        #dal2.min_dot_matches = 5
        #dal2.min_edge_matches = 3 # 2
        ##dal2.align()
        #for hyb in hnf['hyb']:
            #dal.offsets.loc[hyb-1, 'x':] = dal2.align_hyb(hyb)

    
    dal.save_offsets("positions/offsets_ch" + str(channel) + "_pos" + str(pos) + ".csv")
    #dal.save_offsets_SEs("offsets_ch_SEs" + str(channel) + ".csv")
    dal.save_matches("positions/matches_ch" + str(channel) + "_pos" + str(pos) + ".csv")
    dal.save_loocv_errors("positions/loov_errors_ch" + str(channel) + "_pos" + str(pos) + ".csv")
    dal.save_ro_wout_fm("positions/pnts_no_fm_ch" + str(channel) + "_pos" + str(pos) + ".csv")
    
    offsets = dal.offsets
    dal = None
    dal2 = None
    
    offsets.insert(0, column="ch", value=channel)
    all_channel_offsets.append(offsets)
    

agg = pd.concat(all_channel_offsets)
agg.to_csv("positions/" + savefname)