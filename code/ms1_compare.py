import pandas as pd
import numpy as np
from .cluster import compare_clusters_by_coor
from .cluster import merge_cluster_pairs

def normalize_rt_ms1(raw_data, rt_tol=10, ms1_tol=0.01, clusters={}):
    """
    get normalized point coordinate by RT and MS1
    :param raw_data: a pandas data frame
    :param rt_tol: RT tolerance, 10s
    :param ms1_tol: MS1 tolerance, 0.01 Da
    :return: cluster id with point coordinate, a dataFrame
    """
    raw_data.loc[:, 'mzmed'] = raw_data.mzmed / float(ms1_tol)
    raw_data.loc[:, 'predict.rt'] = raw_data.loc[:, 'predict.rt'] / float(rt_tol)
    peak_name2point_coor = {}
    cluster2point_coor = {}
    for i in raw_data.index:
        # name = raw_data.index[i]
        coor = np.array([raw_data.loc[i, 'predict.rt'], raw_data.loc[i, 'mzmed']])
        peak_name2point_coor[i] = coor
    for c in clusters:  # set cluster central coordinate for the first time
        cluster_central_coor = peak_name2point_coor[clusters[c].key_point]
        clusters[c].set_central_coor(cluster_central_coor)
        cluster2point_coor[c] = cluster_central_coor
    # return a dataFrame
    return pd.DataFrame(cluster2point_coor, index=['rt', 'mz']).T


def rt_ms1_compare(raw_data, cluster_collection, rt_tol, ms1_tol):
    """
    compare each peak with all other peaks through rt and mz,
    merge the peaks that are very similar(within rt_tol and ms1_tol)
    :param raw_data: a panda data frame
    :param cluster_collection: total Clusters, a list
    :param rt_tol: 10s
    :param ms1_tol: 15ppm
    :return:
    """
    # cluster_collection = cluster_collection
    # cluster2point_coor: the total clusters' normalized coordinate(rt, mz), a dataFrame
    cluster2point_coor = normalize_rt_ms1(raw_data.copy(), rt_tol=rt_tol, ms1_tol=ms1_tol,
                                          clusters=cluster_collection)
    # all the peaks that were marked as key point and be grouped(without isotopes)
    kg_peaks = raw_data[(raw_data.grouped == 1) & (raw_data.key_point == 1)]  # 2356
    ng_peaks = raw_data[raw_data.grouped == 0]  # all no group peaks  # 13229

    ##  grouped key points dereplication
    # ax = plt.subplot(111)
    # plt.scatter(x=key_points_coor['rt'], y=key_points_coor['mz'])
    # plt.savefig('key_points2.png', dpi=200)
    # compare grouped key points inner
    kg_cluster_coor = cluster2point_coor.loc[pd.unique(kg_peaks['cluster_id'])]
    cluster_pairs_1 = compare_clusters_by_coor(grouped_cluster1_coor=kg_cluster_coor,
                                             cluster2_coor=kg_cluster_coor,
                                             cluster_collection=cluster_collection)
    merge_cluster_pairs(cluster_pairs=cluster_pairs_1['cluster_pairs'],
                        clusters=cluster_collection, raw_data=raw_data,
                        cluster2point_coor=cluster2point_coor)

    ## no group points dereplication

    # compare no grouped points inner
    ng_cluster_coor = cluster2point_coor.loc[pd.unique(ng_peaks['cluster_id'])]
    cluster_pairs_2 = compare_clusters_by_coor(grouped_cluster1_coor=ng_cluster_coor,
                                               cluster2_coor=ng_cluster_coor,
                                               cluster_collection=cluster_collection)
    # merge all clusters
    merge_cluster_pairs(cluster_pairs=cluster_pairs_2['cluster_pairs'],
                        clusters=cluster_collection, raw_data=raw_data,
                        cluster2point_coor=cluster2point_coor)

    # compare between grouped key points and no grouped points
    # update coordinate in kg clusters and ng clusters
    kg_peaks = raw_data[(raw_data.grouped == 1) & (raw_data.key_point == 1)] # 2348
    ng_peaks = raw_data[raw_data.grouped == 0] # 13229
    kg_cluster_coor = cluster2point_coor.loc[pd.unique(kg_peaks['cluster_id'])]
    ng_cluster_coor = cluster2point_coor.loc[pd.unique(ng_peaks['cluster_id'])]  # 13026
    cluster_pairs_3 = compare_clusters_by_coor(grouped_cluster1_coor=kg_cluster_coor,
                                               cluster2_coor=ng_cluster_coor,
                                               cluster_collection=cluster_collection)
    # merge all clusters
    merge_cluster_pairs(cluster_pairs=cluster_pairs_3['cluster_pairs'],
                        clusters=cluster_collection, raw_data=raw_data,
                        cluster2point_coor=cluster2point_coor)

    ng_inx = ng_peaks.index.tolist()
    raw_data.loc[ng_inx, 'grouped'] = 1
    return {'raw_data': raw_data, 'cluster_collection': cluster_collection}