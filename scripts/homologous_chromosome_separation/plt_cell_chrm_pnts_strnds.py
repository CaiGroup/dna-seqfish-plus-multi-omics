#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  2 13:48:09 2021

@author: jonathanwhite
"""

from scipy.spatial import cKDTree as KDTree
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import time

#filename = "LC1-100k-testing-006-2021-05-26-E14-100k-rep2-test-output-finalpoints-paints-decoded.csv"
#filename = "LC1-100k-testing-007-2021-06-29-E14-100k-rep2-test-pos0-output-finalpoints-paints-decoded.csv"
#filename = 'grp_all_chrms_test.csv'
#filename = 'allele_results_all.csv'
#filename = 'results_max_strands2/NMuMG_rep1_pos_9.csv'
#filename = 'sep_res2.csv'
#filename = 'results/NMuMG_rep1_pos_9.csv'
#filename = 'c30chr7_trial.csv'
#filename = 'results/attempt_nbr_fill_E14_rep2_replicate1_pos_1.csv'
#filename = "e14_p1c1chr17_res.csv"
#filename = "e14_res.csv"

filename = 'results/Brain_rep2_pos_3.csv'
filename = '20220201_results/E14_rep3_replicate2_pos_4.csv'
#filename = 'results/E14_rep3_replicate2_pos_8.csv'
#filename = "e14r3r2p8c2c2_res.csv"
filename = 'results/E14_mDuxCa_24hr_rep2_pos_3.csv'
filename = '20220204_results/E14_rep2_replicate1_pos_4.csv'
filename = 'results/E14_rep2_replicate1_pos_4.csv'

filename = 'results/E14_rep3_replicate2_pos_4.csv'



#filename = 'c1chr16_res2.csv'





#ps =  pd.read_csv(filename)
pnts_grpd = pd.read_csv(filename)


#pnts_grpd.loc[:,"x"] *= 103 #nm/pixel
#pnts_grpd.loc[:,"y"] *= 103 #nm/pixel
#pnts_grpd.loc[:,"z"] *= 250 #nm/slice

#pnts_grpd.loc[:,"pos"] = [int(nm.split("-")[1]) for nm in pnts_grpd["name"]]

pnts_grpd.set_index(["cellID", "chrom"], inplace=True)
#chrm = pnts.set_index("chrom").loc['chrX']

atype = {1:'dbscan_allele', 2:'dbscan_ldp_allele', 3:'dbscan_ldp_nbr_allele'}

#def plt_chrom(cell_num, chrm_num, lines=True, a=1):
def plt_chrom(cell_num, chrm_num, a=1):
    #chrm = pnts_grpd.set_index("chrom").loc['chr' + str(chrm_num)]
    chrm = pnts_grpd.loc[cell_num, 'chr' + str(chrm_num)]

    #chrm.sort_values(by='pos', inplace=True)
    chrm.sort_values(by='g', inplace=True)
    if a in atype.keys():
        a = atype[a]

    chrm.reset_index(inplace=True)
    #plt.scatter( chrm['x'],chrm['y'],c=list(chrm['pos']),cmap='Greens')
    plt.scatter( chrm['x'],chrm['y'],c=list(chrm['g']),cmap='Greens')
    #for strand in np.unique(chrm['dbscan_ldp_nbr_allele']):
    if a != 0:
        for strand in np.unique(chrm[a]):
        #for strand in np.unique(chrm['dbscan_ldp_allele']):
        #for strand in np.unique(chrm['dbscan_allele']):
            if strand != -1:
                #plt.plot(chrm.loc[chrm['dbscan_ldp_nbr_allele']==strand,"x"], chrm.loc[chrm['dbscan_ldp_nbr_allele']==strand, "y"],'-')
                plt.plot(chrm.loc[chrm[a]==strand,"x"], chrm.loc[chrm[a]==strand, "y"],'-')

            #plt.plot(chrm.loc[chrm['dbscan_ldp_allele']==strand,"x"], chrm.loc[chrm['dbscan_ldp_allele']==strand, "y"],'-')

            #plt.plot(chrm.loc[chrm['allele']==strand,"x"], chrm.loc[chrm['allele']==strand, "y"],'-')
            #plt.plot(chrm.loc[chrm['dbscan_allele']==strand,"x"], chrm.loc[chrm['dbscan_allele']==strand, "y"],'.')

    plt.show()

#chrm = pnts.set_index("chrom").loc['chr2']