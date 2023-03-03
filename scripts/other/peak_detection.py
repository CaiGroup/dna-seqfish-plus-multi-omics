import pandas as pd
import numpy as np
from scipy.ndimage import uniform_filter1d

def merge_intervals(interval):
    b = []
    for begin,end in sorted(interval):
        if b and b[-1][1] >= begin - 1:
            b[-1][1] = max(b[-1][1], end)
        else:
            b.append([begin, end])
    m_int = np.array(b)
    # may start and end to be int
    start = np.floor(m_int[:, 0]).astype(int)
    end = np.ceil(m_int[:, 1]).astype(int)
    return np.vstack((start, end))


def peak_detection(sdf, marker, chrom, mean_center, prom_thresh, ht_thresh, rel_height, res):
    from scipy.signal import find_peaks, peak_widths
    # smooth the track
    x = sdf[marker].values
    x = uniform_filter1d(x, 4)
    if mean_center == True:
        x = x - x.mean()
    # add begin and end
    x = np.hstack([ [0], x, [0]])
    coord_x = sdf["start"]// res
    coord_x = np.hstack([[coord_x.min() - 1], coord_x, [coord_x.max() + 1]])
    # find peak
    peak, _ = find_peaks(x, prominence = (prom_thresh, None), height = ht_thresh)

    # calculate peak width
    results_half = peak_widths(x, peak, rel_height=rel_height)
    
        # start and end
    start = results_half[1:][1]
    end = results_half[1:][2]
    start = np.maximum(start, np.array([0] * len(start)))
    end = np.minimum(end - 1, np.array([sdf.shape[0] - 2]  * len(start)))
    # merge peaks
    begin_interval = [[start[i], end[i]] for i in range(start.shape[0])]
    
    # with peaks detected
    if len(begin_interval) > 0:
        m_intervals = merge_intervals(begin_interval) 
        if m_intervals.shape == (2,1):
            start = m_intervals[0]
            end = m_intervals[1]
        else:
            start = m_intervals[:, 0]
            end = m_intervals[:, 1]
    
        return (sdf, peak, np.vstack( (np.floor(results_half[1:][1]), np.ceil(results_half[1:][2]))), m_intervals)
    else:
        print (f"{chrom} {marker} has no peak detected")

def stitch_peak(peak_dict, chrom_list, marker, res):
    m_peak_dfs = []
    for chrom in chrom_list:
        if chrom in peak_dict.keys():
            sls = peak_dict[chrom]
            sdf, peak, begin_intervals, m_intervals  = sls[0], sls[1], sls[2], sls[3].astype(int)
            m_peak_start = sdf["start"].values[m_intervals]

            mp_df = pd.DataFrame(m_peak_start).T.rename(columns = {0 : "start", 1: "end"})

            # make it
            mp_df["chrom"] = chrom
            mp_df[f"{res//1000}kb start"] = mp_df["start"] // res
            mp_df[f"{res//1000}kb end"] = mp_df["end"] // res
            m_peak_dfs.append(mp_df)


        else:
            print (f"{marker} {chrom} not exist")
    if len(m_peak_dfs) > 0:
        m_peak_df = pd.concat(m_peak_dfs)
        m_peak_df["len"] = m_peak_df["end"] - m_peak_df["start"]
    return m_peak_df

def melt_interval_df(tdf, res, thresh):
    tdf["len"] = tdf[f"{res//1000}kb end"] - tdf[f"{res//1000}kb start"]
    tdf["chrom_id"] = tdf["chrom"].str[3:].replace("X", 20).astype(int)

    # define short/ long peaks according to chosen threshold
    tdf["short peak"] = (tdf["len"] <= thresh * 200000) * 1

    tdf = tdf.sort_values(by = ["chrom_id", "start"])
    tdf = tdf.reset_index(drop = True).reset_index().rename(columns = {"index" : "peak_id"})
    bins = []
    peak_ids = []
    chroms = []
    shorts = []
    for index, row in tdf.iterrows():
        chrom = row["chrom"]
        start = row[f"{res//1000}kb start"]
        end = row[f"{res//1000}kb end"]
        short = row["short peak"]
        peak_id = row["peak_id"]
        blist = list(range(start, end + 1))
        # print (start, end, blist)
        peak_ids += [peak_id] * len(blist)
        chroms += [chrom] * len(blist)
        bins += blist
        shorts += [short] * len(blist)
    final_peak = pd.DataFrame({f"{res//1000}kb bins" : bins, "peak_id" : peak_ids, "chrom" : chroms, "short" : shorts})
    final_peak[f"{res//1000}kb name"] = final_peak["chrom"] + '-g' + final_peak[f"{res//1000}kb bins"].astype(str)
    return final_peak

# Import the library
import argparse
# Create the parser
parser = argparse.ArgumentParser()
parser.add_argument('--dot_file', type=str, required=True)
parser.add_argument('--meta_file', type=str, required=True)
parser.add_argument('--output_file', type=str, required=True)
parser.add_argument('--marker', type=str, required=True)
parser.add_argument('--peak_params', type=str, required=True)
parser.add_argument('--res', type=int, required=True)
parser.add_argument('--mean_center', type=str, required=False)
parser.add_argument('--thresh', type=int, required=False)

def main():
    args = parser.parse_args()

    dot_file = args.dot_file
    meta_file = args.meta_file
    output_file = args.output_file
    marker = args.marker
    prom_thresh, ht_thresh, rel_height = [float(x) for x in args.peak_params.split(",")]
    res = args.res
    thresh = args.thresh
    mean_center = args.mean_center
    print (f"Detecting peak for marker {marker}")
    chrom_list = [f"chr{x}" for x in range(1, 20)]
    chrom_list.append("chrX")

    # read in dot information
    df = pd.read_csv(dot_file)[[f"{res // 1000}kb name", marker]]
    meta = pd.read_csv(meta_file)
    df = df.merge(meta[["200kb name", "start", "end", "chrom"]], how = "left")
    peak_dict = {}
    for chrom in  chrom_list:
        peak_coord = peak_detection(df[df["chrom"] == chrom], marker, chrom, mean_center, prom_thresh, ht_thresh, rel_height, res)
        if peak_coord:
            peak_dict[chrom] = peak_coord
    tdf = stitch_peak(peak_dict, chrom_list, marker, res)
    melt_interval_df(tdf, res, thresh).to_csv(output_file, index = None)
    print ("All done XD")
    
if __name__ == "__main__":
    main()