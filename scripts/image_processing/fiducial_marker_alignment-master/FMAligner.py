import pandas as pd
import numpy as np
from multiprocessing import Pool
from scipy.spatial import cKDTree as KDTree
from copy import copy, deepcopy

class FMAligner:
    """
    This Class aligns dots in a seqFISH experiment image containing fiducial markers, to a reference image containing
    only the fiducial markers
    """
    def __init__(self, ro, ref, ref_final=None):
        """
        Initialize RefAligner object. Save pandas DataFrames of reference and readout points
        :param ro: pandas dataframe of reodout points with columns = ("x", "y", "z", "int") and index "hyb" denoting the
                hybridiztion in which the point was found.
        :param ref:  pandas dataframe of reference points with columns = ("x", "y", "z", "int")
        :param ref_final: Optional- pandas dataframe of final referaence points measured after readout hybridizations. Used
                for automatically setting parameters in the auto_set_parameters method.
        """

        self.ro = ro
        self.ref = ref
        self.ref['n_trav_matched'] = 0
        self.ref['n_unmatched'] = 0

        self.ref.sort_values(by='y', inplace=True)
        #self.ref.sort_values(by=['x', 'y', 'z', 'int'], inplace=True)

        if type(ref_final) == pd.DataFrame:
            self.ref_final = ref_final

        # define data members that will be used in alignment
        self.matchDict = {}
        self.edgeMatches = set()
        self.trav_matches = set()
        self.ref_dots = None
        self.reserved_ref_dots = None
        self.ro_hyb = None
        self.ro_lat_tree = None
        self.trav_ind = 0
        self.traversal_queue = []
        self.n_ambiguous = 0
        self.n_trav_matched = 0
        self.offsets = None
        self.matchesDF = pd.DataFrame(columns=('ref_x', 'ref_y', 'ref_z', 'hyb', 'comp_x', 'comp_y', 'comp_z', 'comp_int',
                                               'aligned_x', 'aligned_y', 'aligned_z'))
        self.averaged_ref = pd.DataFrame(columns=("ref_y", "hyb", "x", "y", "z", "int"))

        # define threshold parameters set by their own methods
        self.xyse = None
        self.xyse_sq = None
        self.zse = None
        self.min_edge_matches = 5
        self.min_bright_prop = None
        self.max_bright_prop = None
        self.n_unmatch_give_up = 100
        self.min_dot_matches = 10
        #self.n_longest_edges = 100
        self.set_n_longest_edges(30*len(ref.x))#500)
        self.max_lat_offset = None
        self.max_z_offset = None
        self.outlier_sd_thresh = 3
        self.min_fm_hyb_matches = int((max(self.ro.index)-min(self.ro.index))/2)


    def auto_set_params(self):
        """
        Matches the initial reference dots to the final reference dots by trying a range of possible search parameters.
        Uses summary statistics of the alignment errors and differences in brightness to set search parameters to help
        in the search for the fiducial markers in the readout hybridizations.
        :return:
        """
        if type(self.ref_final) != pd.DataFrame:
            raise(Exception("Both initial and final fiducial marker locations must be provided to automatically set parameters."))

        # set search errors
        lat_errors = np.arange(0.4,2,0.2)
        z_errors = np.arange(0.2, 2, 0.2)
        matches = np.zeros([len(lat_errors), len(z_errors)])
        max_matches = 0
        mm_lt = None
        mm_zt = None
        opt_matches = None

        for i, lat_error in enumerate(lat_errors):
            for j, z_error in enumerate(z_errors):

                ref_hyb_indexed = deepcopy(self.ref).loc[:, 'x':'int']
                ref_hyb_indexed.loc[:, 'hyb'] = 1
                ref_hyb_indexed = ref_hyb_indexed.set_index('hyb')
                ral = FMAligner(ref_hyb_indexed, deepcopy(self.ref_final))

                ral.set_min_dot_matches(4)
                ral.set_min_bright_prop(0.1)
                ral.set_max_bright_prop(10)
                ral.set_min_edge_match(2)
                ral.set_n_unmatch_give_up(100)
                ral.set_min_fm_hyb_matches(None)
                ral.set_n_longest_edges(100)
                ral.set_outlier_sd_thresh(3.0)

                ral.set_xy_search_error(lat_error)
                ral.set_z_search_error(z_error)
                ral.align_hyb(1, print_hyb=False)
                matches[i,j] = len(ral.matchesDF.index)
                if len(ral.matchesDF.index) > max_matches:
                    mm_lt = lat_error
                    mm_zt = z_error
                    opt_matches = deepcopy(ral.matchesDF)

        opt_matches["x_diff"] = opt_matches["aligned_x"] - opt_matches["ref_x"]
        opt_matches["y_diff"] = opt_matches["aligned_y"] - opt_matches["ref_y"]

        opt_matches["lat_diff"] = np.sqrt(opt_matches["x_diff"]**2 + opt_matches["y_diff"]**2)
        opt_matches["z_diff"] = opt_matches["aligned_z"] - opt_matches["ref_z"]
        opt_matches["int_diff"]= opt_matches["comp_int"]-opt_matches["ref_int"]

        x_offset = np.mean(opt_matches["comp_x"]-opt_matches["ref_x"])
        y_offset = np.mean(opt_matches["comp_y"]-opt_matches["ref_y"])
        z_offset = np.mean(opt_matches["comp_z"]-opt_matches["ref_z"])
        lat_offset = np.sqrt(x_offset**2 + y_offset**2)

        self.set_z_search_error(np.quantile(np.abs(opt_matches["z_diff"]), 0.9))
        self.set_xy_search_error(np.quantile(opt_matches["lat_diff"], 0.9))

        dfactors = opt_matches["ref_int"]/opt_matches["comp_int"]

        self.set_min_bright_prop(np.min(dfactors/1.5))
        self.set_max_bright_prop(np.max(1/dfactors))
        #if len(opt_matches["ref_y"]) > self.n_longest_edges*4:
            #self.set_n_longest_edges(len(opt_matches["ref_y"])*4)

        #f len(opt_matches["ref_y"]) > self.min_dot_matches*2:
            #self.min_dot_matches = np.round(len(opt_matches["ref_y"])/2)

            #if len(opt_matches["ref_y"]) > self.min_edge_matches*8:
                #self.min_edge_matches = np.round(len(opt_matches["ref_y"])/8)
        
        print(matches)
        return (opt_matches, mm_lt, mm_zt, x_offset, y_offset, z_offset, lat_offset)

    def load_saved_parameters(self, params_file_name):
        self.set_saved_parameters(pd.read_csv(params_file_name))

    def set_saved_parameters(self, params_df):
        self.set_xy_search_error(params_df["xyse"][0])
        self.set_z_search_error(params_df["zse"][0])
        self.set_max_bright_prop(params_df["max_bright_prop"][0])
        self.set_min_bright_prop(params_df["min_bright_prop"][0])
        #self.min_dot_matches = int(params_df["min_dot_matches"][0])
        #self.min_edge_matxhes = int(params_df["min_edge_matches"][0])
        #self.set_n_longest_edges(int(params_df["n_longest_edges"][0]))

        #if params_df["min_dot_matches"][0] > self.min_dot_matches:
            #self.min_dot_matches = int(params_df["min_dot_matches"][0])
            #if params_df["min_edge_matches"][0] > self.min_edge_matches:
                #self.min_edge_matxhes = int(params_df["min_edge_matches"][0])

    # Parameter setting methods
    def set_n_longest_edges(self, nle):
        self.n_longest_edges = nle
        self.ref_edges = self._find_ref_edges(nle)

    def set_xy_search_error(self, xyse):
        self.xyse = xyse
        self.xyse_sq = xyse**2

    def set_z_search_error(self, zse):
        self.zse = zse

    def set_min_edge_match(self, min_edge_matches):
        self.min_edge_matches = min_edge_matches

    def set_min_bright_prop(self, mbp):
        self.min_bright_prop = mbp

    def set_max_bright_prop(self, mbp):
        self.max_bright_prop = mbp

    def set_n_unmatch_give_up(self, nugu):
        self.n_unmatch_give_up = nugu

    def set_min_dot_matches(self,mdm):
        self.min_dot_matches = mdm
        
    def set_max_lat_offset(self, mlo):
        self.max_lat_offset = mlo
        
    def set_max_z_offset(self, mzo):
        self.max_z_offset = mzo

    def set_outlier_sd_thresh(self, ost):
        self.outlier_sd_thresh = ost

    def set_min_fm_hyb_matches(self, mfhm):
        self.min_fm_hyb_matches = mfhm
    
    def align(self, hybs=None):
        """
        aligns hybridizations in the readout points to the reference
        :param hybs: list of hybridizations to align. If None, aligns all hybridizations.
        :return:
        """
        if not hybs:
            #hybs = list(range(1, 1 + max(self.ro.index)))
            hybs = np.unique(self.ro.index)

        offsets = []
        for hyb in hybs:
            hyb_offsets = self.align_hyb(hyb)
            hyb_offsets = [hyb] + hyb_offsets
            offsets.append(hyb_offsets)
        
        self.unfiltered_matches = deepcopy(self.matchesDF)
        self.matchesDF = self._filter_matches(self.matchesDF)
        #self.unaveraged_offsets = pd.DataFrame(offsets, columns=('hyb', 'x', 'y', 'z', 'x_SE', 'y_SE', 'z_SE', 'n_matches'))
        self.offsets = pd.DataFrame(offsets, columns=('hyb', 'x', 'y', 'z', 'x_SE', 'y_SE', 'z_SE', 'n_matches'))
        self._offsets = self._est_offsets_all_matches(self.matchesDF)
        self.loocv_alignment_errors = self._est_error_loocv(self.matchesDF)


        #return self.ref_aved_offsets #.unaveraged_offsets
        #self.offsets.set_index('hyb')
        #if __name__ == '__main__':
            #pool = Pool(n_processes)
            #pool.map(self.align_hyb, hybs)
        return self.offsets, self.loocv_alignment_errors

    def save_ro_wout_fm(self, name):
        matches = deepcopy(self.matchesDF)
        matches["x"] = matches["comp_x"]
        matches["y"] = matches["comp_y"]
        matches["z"] = matches["comp_z"]
        matches["int"] = matches["comp_int"]

        points_ind = self.ro.set_index(["x", "y", "z", "int"], append=True)
        matches_ind = matches.set_index(["hyb", "x", "y", "z", "int"])

        # ToDO: Figure out why this generates errors that we are ignoring for. Possibly missing some fiducial markers
        points_filtered = points_ind.drop(matches_ind.index, errors="ignore")

        assert self.ro.shape[0] - points_filtered.shape[0] <= matches.shape[0] #ToDo: assert ==, not <=

        tosave = points_filtered.reset_index()

        tosave.to_csv(name)

    def save_offsets(self, filename):
        self.offsets.to_csv(filename, index=True)

    def save_loocv_errors(self, filename):
        self.loocv_alignment_errors.to_csv(filename, index=False)

    def save_matches(self, save_name):
        self.matchesDF.to_csv(save_name, index=False)

    def save_params(self, save_name):
        params = pd.DataFrame()
        params["xyse"] = [self.xyse]
        params["zse"] = [self.zse]
        params["max_bright_prop"] = [self.max_bright_prop]
        params["min_bright_prop"] = [self.min_bright_prop]
        params["min_dot_matches"] = [self.min_dot_matches]
        params["min_edge_matches"] = [self.min_edge_matches]
        params["n_longest_edges"] = [self.n_longest_edges]
        params.to_csv(save_name, index=False)

    def _est_offsets_all_matches(self, _matches):
        """
        From a dataframe of fiducial marker matches, estimate offsets
        :param _matches: DataFrame of matched fiducial markers
        :return: DataFrame of offsets
        """
        matches = _matches.set_index('hyb')
        matches["x_offset"] = matches['comp_x'] - matches['ref_x']
        matches["y_offset"] = matches['comp_y'] - matches['ref_y']
        matches["z_offset"] = matches['comp_z'] - matches['ref_z']
        mfg = matches.groupby('hyb')
        means = mfg.mean()
        offsets = means[['x_offset', 'y_offset', 'z_offset']]
        return offsets

    def _est_error_loocv(self, _matches):
        """
        Estimate the average alignment error of dots between hybridizations using leave-one-out cross validation on the
        matched fiducial markers
        :param _matches: DataFrame of matched fiducial markers
        :return: dataframe of loocv alignments for each fiducial marker to the reference and the alignment errors
        """

        matches = _matches.set_index(["ref_x", "ref_y", "ref_z", "ref_int"], drop=False)
        loocv_alignments = pd.DataFrame(columns=[ "ref_x", "ref_y", "ref_z", "ref_int", "hyb",
                                                "loocv_align_x", "loocv_align_y", "loocv_align_z", "comp_int",
                                               "x_diff", "y_diff", "z_diff"])

        for fm_ref in np.unique(matches.index):
            drop_fm_matches = matches.drop(fm_ref)
            loo_offsets = self._est_offsets_all_matches(drop_fm_matches)

            # align fm with offests calculated from other fiducial markers
            drpd_fm = matches.loc[fm_ref]
            drpd_fm.set_index("hyb", drop=False, inplace=True)
            drpd_fm["loocv_align_x"] = drpd_fm["comp_x"] - loo_offsets.loc[drpd_fm["hyb"]]["x_offset"]
            drpd_fm["loocv_align_y"] = drpd_fm["comp_y"] - loo_offsets.loc[drpd_fm["hyb"]]["y_offset"]
            drpd_fm["loocv_align_z"] = drpd_fm["comp_z"] - loo_offsets.loc[drpd_fm["hyb"]]["z_offset"]
            drpd_fm["x_diff"] = drpd_fm["loocv_align_x"] - drpd_fm["ref_x"]
            drpd_fm["y_diff"] = drpd_fm["loocv_align_y"] - drpd_fm["ref_y"]
            drpd_fm["z_diff"] = drpd_fm["loocv_align_z"] - drpd_fm["ref_z"]

            loocv_alignments = loocv_alignments.append(drpd_fm[["ref_x", "ref_y", "ref_z", "ref_int", "hyb",
                                               "loocv_align_x", "loocv_align_y", "loocv_align_z", "comp_int",
                                                "x_diff", "y_diff", "z_diff"]])

        return loocv_alignments

    def _find_ref_edges(self, n_longest, ref_df = None):
        """
        Returns the n_longest edges (position separation vectors) in the graph of reference dots.
        :param n_longest: the number of longest edges to report
        :param ref_df: optional- pandas DataFrame of reference points. If not supplied, uses self.ref.
        :return:
        """
        n_corner = int(np.ceil(np.sqrt(n_longest)))
        # find dots closest to each corner
        print('Finding Reference Edges')
        if ref_df is None:
            ref_dists = copy(self.ref)
        else:
            ref_dists = copy(ref_df)
        ref_dists['ul_dist'] = 0
        ref_dists['ur_dist'] = 0
        ref_dists['bl_dist'] = 0
        ref_dists['br_dist'] = 0

        #find distances the distance between each dot and each of the four corners of the image
        for i, dot in self.ref.iterrows():
            #ref_dists['ul_dist'].loc[i] = dot.x**2 + dot.y**2
            #ref_dists['ur_dist'].loc[i] = (2048 - dot.x)**2 + dot.y**2
            #ref_dists['bl_dist'].loc[i] = dot.x**2 + (2048 - dot.y)**2
            #ref_dists['br_dist'].loc[i] = (2048 - dot.x)**2 + (2048 - dot.y)**2

            ref_dists.loc[i, 'ul_dist'].loc[i] = dot.x ** 2 + dot.y ** 2
            ref_dists.loc[i, 'ur_dist'].loc[i] = (2048 - dot.x) ** 2 + dot.y ** 2
            ref_dists.loc[i, 'bl_dist'].loc[i] = dot.x ** 2 + (2048 - dot.y) ** 2
            ref_dists.loc[i, 'br_dist'].loc[i] = (2048 - dot.x) ** 2 + (2048 - dot.y) ** 2

        ul_dots = ref_dists.sort_values(by='ul_dist').iloc[:n_corner]
        ur_dots = ref_dists.sort_values(by='ur_dist').iloc[:n_corner]
        bl_dots = ref_dists.sort_values(by='bl_dist').iloc[:n_corner]
        br_dots = ref_dists.sort_values(by='br_dist').iloc[:n_corner]

        edges_array = np.zeros([2*n_corner**2, 12])
        i = 0
        for j, ul_dot in ul_dots.iterrows():
            for k, br_dot in br_dots.iterrows():
                edges_array[i, 0:4] = np.array(ul_dot['x':'int'])
                edges_array[i, 4:8] = np.array(br_dot['x':'int'])

                # now add col_dist, y_dist and absolute distance
                x_dist = br_dot['x'] - ul_dot['x']
                y_dist = br_dot['y'] - ul_dot['y']
                z_dist = br_dot['z'] - ul_dot['z']
                length = np.sqrt(y_dist ** 2 + x_dist ** 2)
                edges_array[i, 8:12] = np.array([y_dist, x_dist, z_dist, length])
                i += 1

        for j, ur_dot in ur_dots.iterrows():
            for k, bl_dot in bl_dots.iterrows():
                edges_array[i, 0:4] = np.array(ur_dot['x':'int'])
                edges_array[i, 4:8] = np.array(bl_dot['x':'int'])

                # now add x_dist, y_dist and absolute distance
                x_dist = bl_dot['x'] - ur_dot['x']
                y_dist = bl_dot['y'] - ur_dot['y']
                z_dist = bl_dot['z'] - ur_dot['z']
                length = np.sqrt(y_dist ** 2 + x_dist ** 2)
                edges_array[i, 8:12] = np.array([x_dist, y_dist, z_dist, length])
                i += 1

        edges_df = pd.DataFrame(data=edges_array,
                              columns=['u_x','u_y',  'u_z', 'u_int', 'v_x', 'v_y', 'v_z', 'v_int',
                                    'x_dist', 'y_dist', 'z_dist', 'length'])
        edges_df.sort_values(by='length', ascending=False, inplace=True)

        edges_df = edges_df.iloc[:n_longest]

        return edges_df

    #def _match_hyb_fiducial_markers(self, hyb):
    def align_hyb(self, hyb, print_hyb = True):
        """
        Finds the alignment of a readout hybridization to the reference
        :param hyb: integer number of the hybridization to align (if not supplying the to_align dataframe)
        :param ref: Optional- pandas DataFrame of the reference dots. If not given, uses self.ref
        :param to_align: Optional - pandas DataFrame of readout dots to align. if not given, aligns the hyb number given.
        :return: list: [x_mean_offset, y_mean_offset, z_mean_offset, x_mean_offset_se, y_mean_offset_se, z_mean_offset_se, <number of dots used to calculate offsets>]
        """
        if print_hyb:
            print('Aligning hybridization', hyb)

        self.ro_hyb = copy(self.ro.loc[hyb])

        self.ro_lat_tree = KDTree(self.ro_hyb.loc[:, 'x':'y'])

        # find a pair of dots in the hybridization readout that match a pair of the reference fiducial markers
        i = 0
        for edge_ind, edge in self.ref_edges.iterrows():
            #print(edge)
            #print('edge number:', edge_ind)
            if self._find_ro_matching_pair(edge):
                break
            i += 1
            if i >= self.n_longest_edges:
                print("Giving up search for matching edges")
                #break
                return [None, None, None, None, None, None, 0]
        offsets_ses = self._est_offsets()
        self._add_hyb_to_match_df(hyb, offsets_ses)
        print('offsets and SEs:', offsets_ses)
        return offsets_ses

    def _find_ro_matching_pair(self, edge):
        """
        Given an edge in the reference graph, look for the corresponding matching edge in the dataframe
        of dots in the readout image. First rule out readout dots based on search parameters, then search candidate
        dots for matching edge vector.
        :param edge: Pandas series with data representing an edge: 'u_x', 'u_y', 'u_z', 'u_int', 'u_x', 'u_y',
        'u_z', 'u_int', 'x_dist', 'y_dist', 'length'
        :return:
        """
        #print('subsetting readout...')
        # ToDo: in C or Cython, subset with a single for loop checking each element for each threshold, and passing if one is failed

        # subset points in read out image to search through by excluding dots with y values that preclude
        # them from being matched to the reference edge
        ###### u is the dot with lower y, v is the dot with higher y
        #change in new version
        # u is the dot with lower x, v is the dot with higher x
        u_r_max = 2048 - edge['x_dist'] + self.xyse
        v_r_min = edge['x_dist'] - self.xyse

        ro_u_candidates = self.ro_hyb.iloc[list(self.ro_hyb['y'] < u_r_max)]
        ro_v_candidates = self.ro_hyb.iloc[list(self.ro_hyb['y'] > v_r_min)]

        if self.max_lat_offset:
            u_allowed = list((edge['u_x'] - ro_u_candidates['x'])**2 + (edge['u_y'] - ro_u_candidates['y'])**2 < self.max_lat_offset**2)
            v_allowed = list((edge['v_x'] - ro_v_candidates['x'])**2 + (edge['v_y'] - ro_v_candidates['y'])**2 < self.max_lat_offset**2)
            ro_u_candidates = ro_u_candidates.iloc[u_allowed]
            ro_v_candidates = ro_v_candidates.iloc[v_allowed]

        if self.max_z_offset:
            u_allowed = list(abs(edge['u_z'] - ro_u_candidates['z']) < self.max_lat_offset)
            v_allowed = list(abs(edge['v_z'] - ro_v_candidates['z']) < self.max_lat_offset)
            ro_u_candidates = ro_u_candidates.iloc[u_allowed]
            ro_v_candidates = ro_v_candidates.iloc[v_allowed]

        ##### exclude dots whose x values preclude them from being matched to the reference edge
        # exclude dots whose y values preclude them from being matched to the reference edge
        if float(edge['u_y']) < float(edge['v_y']):
            # u has a lower x than v, edge x distance is positive
            u_c_max = 2048 - edge['y_dist'] + self.xyse
            v_c_min = edge['y_dist'] - self.xyse

            ro_u_candidates = ro_u_candidates.iloc[list(ro_u_candidates['y'] < u_c_max)]
            ro_v_candidates = ro_v_candidates.iloc[list(ro_v_candidates['y'] > v_c_min)]

        else:
            ##### u has a larger x than v, x distance is negative
            # u has a larger y than v, y distance is negative
            u_c_min = -edge['y_dist'] - self.xyse
            v_c_max = 2048 + edge['y_dist'] + self.xyse

            ro_u_candidates = ro_u_candidates.iloc[list(ro_u_candidates['y'] > u_c_min)]
            ro_v_candidates = ro_v_candidates.iloc[list(ro_v_candidates['y'] < v_c_max)]

        # subset by brightness
        ucands_int_gt_lbnd = np.greater(np.array(ro_u_candidates['int']), float(edge['u_int']) * self.min_bright_prop)
        ucands_int_lt_ubnd = np.less(np.array(ro_u_candidates['int']), float(edge['u_int']) * self.max_bright_prop)
        ro_u_candidates = ro_u_candidates.loc[list(np.logical_and(ucands_int_gt_lbnd, ucands_int_lt_ubnd))]

        v_cands_int_gt_lbnd = np.greater(np.array(ro_v_candidates['int']), float(edge['v_int']) * self.min_bright_prop)
        v_cands_int_lt_ubnd = np.less(np.array(ro_v_candidates['int']), float(edge['v_int']) * self.max_bright_prop)
        ro_v_candidates = ro_v_candidates.loc[list(np.logical_and(v_cands_int_gt_lbnd, v_cands_int_lt_ubnd))]
        
        ro_u_candidates.sort_values(by='x', axis='index', inplace=True)
        ro_v_candidates.sort_values(by='x', axis='index', inplace=True)

        # minimum viable ro v candidate for match to reference edge
        min_v = 0

        # when reference graph is sufficiently traversed, set true to break
        done = False

        for i, udot in ro_u_candidates.iterrows():
            # print(minV, 'of', len(v_cands))
            for j in range(min_v, len(ro_v_candidates)):
                xdist = ro_v_candidates.iloc[j]['x'] - udot['x']
                ydist = ro_v_candidates.iloc[j]['y'] - udot['y']
                zdist = ro_v_candidates.iloc[j]['z'] - udot['z']

                # initial two if elif statements improve runtime
                if xdist < edge['x_dist'] - self.xyse:
                    min_v += 1

                elif xdist > edge['x_dist'] + self.xyse:
                    # print('rdist to vcand out of bounds. Break')
                    break

                elif (xdist - edge['x_dist'])**2 + (ydist - edge['y_dist'])**2 <= self.xyse_sq and \
                        np.abs(zdist - edge['z_dist']) <= self.zse:
                    # found match! run graph traversal
                    done = self._traverse_reference(edge, udot, ro_v_candidates.iloc[j])
                    if done:
                        break
            if done:
                break
        return done

    def _traverse_reference(self, matched_ref_edge, ro_u, ro_v):
        """
        From an initial matched edge dots in the reference image and the reodout image, traverse the
        reference graph in search of the rest of the fiducial markers in the readout image.
        :param matched_ref_edge: Pandas series with data representing an edge: 'u_x', 'u_y', 'u_z', 'u_int', 'u_x', 'u_y',
        'u_z', 'u_int', 'x_dist', 'y_dist', 'length'
        :param ro_u: Pandas series representing u dot matched to reference edge. Has parameters: x, y, z, int
        :param ro_v: Pandas series representing v dot matched to reference edge. Has parameters: x, y, z, int
        :return:
        """
        #print("Starting Traversal")

        # clear match Dict and edge matches
        self.matchDict = {}
        self.edgeMatches = set()
        self.ref_dots = copy(self.ref)
        self.ref_dots.set_index(['x', 'y', 'z', 'int'], drop=False, inplace=True)
        self.reserved_ref_dots = pd.DataFrame(columns=self.ref.columns)
        self.ro_hyb['matches'] = 0

        ref_u_in = matched_ref_edge['u_x':'u_int']
        ref_u_in.rename({'u_x': 'x', 'u_y': 'y', 'u_z': 'z', 'u_int': 'int'}, inplace=True)
        ref_v_in = matched_ref_edge['v_x':'v_int']
        ref_v_in.rename({'v_x': 'x', 'v_y': 'y', 'v_z': 'z', 'v_int': 'int'}, inplace=True)

        # each entry in the traversal queue is a tuple. The first entry is a pd.Series with info on the ref dot,
        # and the second entry is a pd.Series with info in the readout dot
        self.traversal_queue = [(ref_u_in, ro_u), (ref_v_in, ro_v)]
        self.n_ambiguous = 0
        self.n_trav_matched = 0
        self.total_unmatched = 0
        self.trav_ind = 0

        while self.traversal_queue:
            matched_dot = self.traversal_queue.pop(0)
            if tuple(matched_dot[0][:4]) in self.ref_dots.index:
                self.ref_dots.drop(tuple(matched_dot[0][:4]), inplace=True)
                self._find_ro_neighbors(matched_dot)
            
            matched_dot_key = _to_tuples(matched_dot)
            if matched_dot_key in self.matchDict and self.matchDict[matched_dot_key] < self.min_edge_matches:
                self.ref_dots.append(matched_dot[0])

        n_well_matched = self._n_well_matched()
        found_enough = n_well_matched > self.min_dot_matches

        print('n well matched:', n_well_matched, '; n trav matched:', self.n_trav_matched,
              '; n unmatched:', self.total_unmatched, '; n_ambiguous:', self.n_ambiguous)

        return found_enough

    def _find_ro_neighbors(self, matched_dot):
        """
        For a dot matched in the reference and the readout image, search for neighbors in the readout image.
        :param matched_dot: tuple of pandas series representing dots in the reference frame and the readout frame:
            (ref_dot = pd.Series(x, y, z, int), ro_dot = pd.Series(x, y, z, int))
        :return: tuple of pandas series representing dots in the reference frame and the readout frame:
            (ref_dot = pd.Series(x, y, z, int), ro_dot = pd.Series(x, y, z, int))
        """
        self.dot_neighbor_matched = 0
        self.dot_neighbor_unmatched = 0
        matched_dot_key = _to_tuples(matched_dot)

        #print('Searching through', len(self.ref_dots), 'neighbors in reference.')
        i = 0
        keep_going = True
        while i < len(self.ref_dots):
            ref_neighbor = self._trav_next()
            keep_going = self._compare_dots(matched_dot, matched_dot_key, ref_neighbor)
            if not keep_going:
                break
            i += 1
        # if made it to here without making enough matches, see if we can find enough matches in the reserved dots
        if keep_going:
            for i, reserved_ref_dot in self.reserved_ref_dots.iterrows():
                keep_going = self._compare_dots(matched_dot, matched_dot_key, reserved_ref_dot, reserved=True)
                if not keep_going:
                    break

    def _compare_dots(self, matched_dot, matched_dot_key, ref_neighbor, reserved=False):
        """
        Starting with a matched pair of dots between the reference image and the readout image, attempt to find a match
        for another reference fiducial marker dot in the read out image.
        :param matched_dot: tuple of pandas series giving the the coordinates and intensity of the dot in teh reference image
                and its match in the readout image
        :param matched_dot_key: multiindex tuple key for a dot in the matches dictionary (x, y, z, int)
        :param ref_neighbor: Next fiducial marker to search for in the readout image
        :param reserved: set true if looking for neighbors starting from reserved dots
        :return:
        """

        # find position vector separating matched fiducial dot from neighbor in reference image
        xdist = ref_neighbor.x - matched_dot[0].x
        ydist = ref_neighbor.y - matched_dot[0].y
        zdist = ref_neighbor.z - matched_dot[0].z

        # search for neighbor in read out hyb
        search_x = matched_dot[1].x + xdist
        search_y = matched_dot[1].y + ydist
        search_point = [search_x, search_y]

        # search for two closest dots because it is faster than doing a ball search
        nn_dists, nn_inds = self.ro_lat_tree.query(search_point, 2)
        within_error_inds = [nn_inds[i] for i in (0, 1) if nn_dists[i] <= self.xyse and np.abs(self.ro_hyb.iloc[nn_inds[i]]['z'] - zdist - matched_dot[1].z) <= self.zse]
        matches = self.ro_hyb.iloc[within_error_inds]

        if len(matches) == 1:  # we found a match!
            matches = matches.iloc[0]
            new_match = (ref_neighbor['x':'int'], matches['x':'int'])
            if _to_tuples(new_match) not in self.edgeMatches:
                self.n_trav_matched += 1
                self.dot_neighbor_matched += 1
                dropped_neighbor = self._process_match(new_match, matched_dot, reserved)
                if dropped_neighbor:
                    return True
            ref_neighbor_key = tuple(ref_neighbor['x':'int'])
            if not reserved:
                self.ref_dots.loc[ref_neighbor_key, 'n_trav_matched'] += 1

        elif len(matches) > 1:  # ambiguous, add to counter and move on
            self.n_ambiguous += 1
        else:
            self.total_unmatched += 1
            self.dot_neighbor_unmatched += 1
            ref_neighbor_key = tuple(ref_neighbor['x':'int'])
            if not reserved:
                self.ref_dots.loc[ref_neighbor_key, 'n_unmatched'] += 1
                try:
                    no_matches = np.equal(np.array(self.ref_dots.loc[ref_neighbor_key, 'n_trav_matched']), 0)[0]
                except:
                    no_matches = np.equal(self.ref_dots.loc[ref_neighbor_key, 'n_trav_matched'], 0)#
                """
                try:
                    too_many_failed_searches = np.equal(np.array(self.ref_dots.loc[ref_neighbor_key, 'n_unmatched']), self.n_unmatch_give_up)[0]
                except:
                    too_many_failed_searches = np.equal(self.ref_dots.loc[ref_neighbor_key, 'n_unmatched'],
                                                        self.n_unmatch_give_up)

                if no_matches and too_many_failed_searches:
                    drop_ind = self.ref_dots.index.get_loc(ref_neighbor_key)
                    if type(drop_ind) == slice:
                        drop_ind = drop_ind.start
                    self.ref_dots = self.ref_dots.drop(ref_neighbor_key)
                    if drop_ind < self.trav_ind:
                        self.trav_ind -= 1
                    # print('removing', ref_neighbor_key, 'too many failed searches.')
                    return True
                """

            if self.dot_neighbor_matched == 0 and self.dot_neighbor_unmatched >= self.n_unmatch_give_up:
                # print(self.n_umnatch_give_up, 'unmatched: breaking.')
                # return n_matched, n_unmatched
                return False
        if matched_dot_key in self.matchDict and self.matchDict[matched_dot_key] >= self.min_edge_matches:
            # print('Matched enough after looking through', i, 'neighbors. Moving on.')
            # return n_matched, n_unmatched
            return False
        return True

    def _trav_next(self):
        """
        Get the next fiducial dot from the reference image to look for in the readout image
        :return:
        """
        if self.trav_ind >= len(self.ref_dots):
            self.trav_ind = 0
        to_return = self.ref_dots.iloc[self.trav_ind]
        self.trav_ind += 1
        return to_return

    def _process_match(self, new_match, old_match, reserved=False):
        """
        Takes dots matched by a KDTree look up to within self.xyse and within self.zse, does the necessary accounting:
        Adds the new dot the traveral queue if it has not already been added. Adds a vote for each dot match in the edge
        towards being included in the offset calculation
        :param new_match: tuple of pandas series representing dots in the reference frame and the readout frame for the new match:
            (ref_dot = pd.Series(x, y, z, int), ro_dot = pd.Series(x, y, z, int))
        :param old_match: tuple of pandas series representing dots in the reference frame and the readout frame for the \
            previously known match:
            (ref_dot = pd.Series(x, y, z, int), ro_dot = pd.Series(x, y, z, int))
        :param reserved: Set to true if processing a match found by looking through reserve dots
        :return: True if drops neighbor upon neighbor reaching min_edge_matches
        """
        
        new_match_key = _to_tuples(new_match)
        old_match_key = _to_tuples(old_match)
        
        # if match has reached this point, is within traveral error; add to traveral queue is not
        if new_match_key not in self.trav_matches:
            self.traversal_queue.append(new_match)
            self.trav_matches |= {new_match_key}

        if (new_match_key[0], old_match_key[0]) not in self.edgeMatches:
            if new_match_key in self.matchDict:
                self.matchDict[new_match_key] += 1
            else:
                self.matchDict[new_match_key] = 1
            if old_match_key in self.matchDict:
                self.matchDict[old_match_key] += 1
            else:
                self.matchDict[old_match_key] = 1
            self.edgeMatches |= {(old_match_key[0], new_match_key[0]), (new_match_key[0], old_match_key[0])}

            # remove from ref_dots dataframe if reached minimum number of matches to save time searching later
            if self.matchDict[old_match_key] >= self.min_edge_matches and old_match_key[0] in self.ref_dots.index and not reserved:
                self.reserved_ref_dots = self.reserved_ref_dots.append(self.ref_dots.loc[old_match[0]])
                drop_ind = self.ref_dots.index.get_loc(old_match_key[0])
                if type(drop_ind) == slice:
                    drop_ind = drop_ind.start
                self.ref_dots = self.ref_dots.drop(old_match_key[0])
                if drop_ind < self.trav_ind:
                    self.trav_ind -= 1
            if self.matchDict[new_match_key] >= self.min_edge_matches and new_match_key[0] in self.ref_dots.index and not reserved:
                self.reserved_ref_dots = self.reserved_ref_dots.append(self.ref_dots.loc[new_match_key[0]])
                try:
                    drop_ind = self.ref_dots.index.get_loc(new_match_key[0])
                    if type(drop_ind) == slice:
                        drop_ind = drop_ind.start
                    if drop_ind < self.trav_ind:
                        self.trav_ind -= 1
                except AttributeError as e:
                    print(e)

                self.ref_dots = self.ref_dots.drop(new_match_key[0])
                return True

            # reset self.trav_ind if now out of bounds
            if self.trav_ind >= len(self.ref_dots):
                self.trav_ind = 0
        return False

    def _n_well_matched(self):
        """
        Find the number of dots that have been matched by at least the minimum number of edges to be included in
        the offset calculation (self.min_edge_matches).
        :return: integer number edges meeting the threshold to be included in the offset calculation.
        """
        n_well_matched = len([match for match in self.matchDict if self.matchDict[match] >= self.min_edge_matches])
        return n_well_matched

    def _filter_matches(self, matches):
        """
        Removes fiducial markers that were matched in fewer than self.min_fm_hyb_matches hybridizations from the input
        matches dataframe.
        :param matches: Pandas Dataframe of matches for which to remove poor quality fiducial markers.
        :return: subset of input dataframe with poor quality fiduicial markers removed.
        """
        mdf_grp = matches.groupby(["ref_x", "ref_y", "ref_z", "ref_int"])
        sizes = mdf_grp.size()
        matches_filtered = matches.set_index(["ref_x", "ref_y", "ref_z", "ref_int"], drop=False)
        if self.min_fm_hyb_matches:
            ref_to_remove = sizes.index[sizes < self.min_fm_hyb_matches]
            if len(ref_to_remove) == 0:
                return matches_filtered
            return matches_filtered.drop(ref_to_remove)

        return matches_filtered


    def _est_offsets(self):
        """
        Estimates offsets for a hybridization using all fiducial markers matched in that hybridization.
        :return:
        """

        dot_x_offsets, dot_y_offsets, dot_z_offsets, x_mean, xstdv, y_mean, ystdv = self._get_match_coords()

        z_mean = np.mean(dot_z_offsets)
        zstdv = np.std(dot_z_offsets)

        x_mean_se = xstdv / np.sqrt(len(dot_x_offsets))
        y_mean_se = ystdv / np.sqrt(len(dot_y_offsets))
        z_mean_se = zstdv / np.sqrt(len(dot_z_offsets))

        return [x_mean, y_mean, z_mean, x_mean_se, y_mean_se, z_mean_se, len(dot_y_offsets)]


    def _get_match_coords(self):
        """
        Helper method for _est_offsets
        :return:
        """
        xdisp = []
        ydisp = []
        zdisp = []
        nmatches = []
        for match in self.matchDict:
            # only consider dots that were matched more than the minimum allowed number of times
            if self.matchDict[match] >= self.min_edge_matches:
                xdisp.append(match[1][0] - match[0][0])
                ydisp.append(match[1][1] - match[0][1])
                zdisp.append(match[1][2] - match[0][2])
                nmatches.append(self.matchDict[match])

        dot_x_offsets = np.array(xdisp)
        dot_y_offsets = np.array(ydisp)
        dot_z_offsets = np.array(zdisp)
        n_dot_edges = np.array(nmatches)

        # remove any outliers in XY that my have been introduced by mismatches. Don't consider z since z measurements are poor.
        outlier = True
        while outlier:
            # compute summary stats
            x_mean = np.mean(dot_x_offsets)
            y_mean = np.mean(dot_y_offsets)
            xstdv = np.std(dot_x_offsets)
            ystdv = np.std(dot_y_offsets)

            # check if there are lateral outliers
            outliers = np.logical_or(abs(dot_y_offsets - y_mean) > self.outlier_sd_thresh * ystdv, abs(dot_x_offsets - x_mean) > self.outlier_sd_thresh * xstdv)
            outlier = np.any(outliers)

            # remove outliers
            dot_y_offsets = dot_y_offsets[np.logical_not(outliers)]
            dot_x_offsets = dot_x_offsets[np.logical_not(outliers)]
            n_dot_edges = n_dot_edges[np.logical_not(outliers)]

        return dot_x_offsets, dot_y_offsets, dot_z_offsets, x_mean, xstdv, y_mean, ystdv

    def _add_hyb_to_match_df(self, hyb_num, offsets):
        '''
        returns dataframe of fidicuial coordinates in reference image, matched coordinates in hyb images,
        and aligned coordinates off in hyb images. Columns indexed by hyb and then number in hyb
        :param hyb: integer current hyb number
        :return:
        Dataframe as described above
        '''
        self.ref.set_index(['x', 'y', 'z'], inplace=True, drop=False)
        x_offset, y_offset, z_offset, x_mean_se, y_mean_se, z_mean_se, n_dots = offsets
        ref_x = []
        ref_y = []
        ref_z = []
        ref_int = []
        comp_x = []
        comp_y = []
        comp_z = []
        comp_int = []
        for match in self.matchDict:
            if self.matchDict[match] >= self.min_edge_matches:
                ref_x.append(match[0][0])
                ref_y.append(match[0][1])
                ref_z.append(match[0][2])
                ref_int.append(match[0][3])

                comp_x.append(match[1][0])
                comp_y.append(match[1][1])
                comp_z.append(match[1][2])
                comp_int.append(match[1][3])

        matchesDF = pd.DataFrame()  # index=index
        matchesDF['hyb'] = [hyb_num] * len(ref_y)
        matchesDF['ref_x'] = ref_x
        matchesDF['ref_y'] = ref_y
        matchesDF['ref_z'] = ref_z
        matchesDF['ref_int'] = ref_int
        matchesDF['comp_x'] = comp_x
        matchesDF['comp_y'] = comp_y
        matchesDF['comp_z'] = comp_z
        matchesDF['comp_int'] = comp_int
        matchesDF['aligned_x'] = matchesDF['comp_x'] - x_offset
        matchesDF['aligned_y'] = matchesDF['comp_y'] - y_offset
        matchesDF['aligned_z'] = matchesDF['comp_z'] - z_offset

        self.matchesDF = self.matchesDF.append(matchesDF)

        return matchesDF


def _to_tuples(match):
    """
    :param match:tuple of pandas series representing dots in the reference frame and the readout frame for the new match:
        (ref_dot = pd.Series(x, y, z, int), ro_dot = pd.Series(x, y, z, int))
    :return:
    """
    if type(match[0]) == pd.Series:
        ref_tup = tuple(match[0])
    else:
        ref_tup = tuple(match[0].index[0])
    if type(match[1]) == pd.Series:
        ro_tup = tuple(match[1])
    else:
        ro_tup = tuple(match[1].index[0])
    return (ref_tup, ro_tup)