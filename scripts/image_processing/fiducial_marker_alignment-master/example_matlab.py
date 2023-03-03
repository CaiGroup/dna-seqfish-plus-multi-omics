import pandas as pd
import os
import sys

from FMAligner import FMAligner

home_dir = sys.argv[1]
init_ref_fname = sys.argv[2]
final_ref_fname = sys.argv[3]
ro_fname = sys.argv[4]
savefname = sys.argv[5]


os.chdir(home_dir)

ref_final = pd.read_csv(final_ref_fname)
ref_init = pd.read_csv(init_ref_fname)
ro = pd.read_csv(ro_fname)

channels = [1, 2, 3]

for channel in channels:
    ch_ref = ref_init.set_index("ch").loc[channel]
    ch_ro = ro.set_index("ch").loc[channel].set_index("hyb")
    ch_ref_final = ref_final.set_index("ch").loc[channel]

    dal = FMAligner(ch_ro, ch_ref, ch_ref_final)
    dal.auto_set_params()
    dal.set_max_lat_offset(50)
    dal.set_n_longest_edges(5000)
    dal.set_max_z_offset(8)
    dal.min_dot_matches = 5
    dal.min_edge_matches = 2

    dal.align()
    dal.save_offsets(savefname + "offsets_ch" + str(channel) + ".csv")
    dal.save_matches(savefname + "matches_ch" + str(channel) + ".csv")
    dal.save_loocv_errors(savefname + "loov_errors_ch" + str(channel) + ".csv")
    dal.save_ro_wout_fm(savefname + "pnts_no_fm_ch" + str(channel) + ".csv")