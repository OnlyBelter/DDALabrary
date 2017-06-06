import pandas as pd
import numpy as np
# import matplotlib
# import matplotlib.pyplot as plt
from scipy.spatial import distance
import os
import re
from collections import defaultdict

# matplotlib.style.use('ggplot')
root_dir = r'D:\vm_share_folder\DDA_lib\01RT_QC\result'
file_name = 'corrected_pos_1_0.csv'

def read_file(root_d, file_n):
    """
    read file after RT correction
    :param root_d: root dir
    :param file_n: file name
    :return: file data
    """
    file_path = os.path.join(root_d, file_n)
    raw_data = pd.read_csv(file_path)
    part_data = raw_data.loc[:, ['X', 'name', 'mzmed', 'isotopes', 'adduct', 'predict.rt']]
    return(part_data)


class Cluster():
    def __init__(self, cluster_id, sample_name):
        # (id, sample_name) is the primary key for a cluster
        self.id = cluster_id
        self.members = set()  # peak name
        self.isotope = set()
        self.key_point = ''  # has min mz
        self.central_coor = np.array([0, 0])  # normalized rt and mz
        self.ms2_num = []  # ms2 peak number
        self.ms2_info = []  # may has more than one ms2 info, [{'mz': 'intensity'}, {}, ...]
        self.sample_name = sample_name  # sample name, a string

    def add_member(self, new_member):
        # new_member is a set
        self.members = self.members.union(new_member)

    def add_isotope(self, isotope):
        # isotope is a set
        self.isotope = self.isotope.union(isotope)

    def set_key_point(self, key_point):
        self.key_point = key_point

    def set_central_coor(self, coor):
        # coor is a numpy array with two dimensions
        self.central_coor = coor

    def set_sample_name(self, sample_name):
        # sample name is a string
        self.sample_name = sample_name

    def add_ms2_num(self, num):
        # num is a list
        self.ms2_num.append(num)

    def add_ms2_info(self, ms2_info):
        # ms2_info is a dict
        self.ms2_info.append(ms2_info)

    def get_uniso_mem(self): # un isotope member
        return self.members - self.isotope


class ClusterCollection():

    def create_new_cluster_id(self):
        pass


def de_isotope(raw_data, sample_name):
    raw_data['key_point'] = 0  # 0 means this point isn't a key point
    raw_data['grouped'] = 0  # 0 means this point isn't grouped
    raw_data['is_iso'] = 0
    raw_data['cluster_id'] = 'no_group'
    cluster_collection = {}  # cluster id to class Cluster
    name2cluster_id = {}
    cache = {}
    for i in range(len(raw_data)):
        if_null = raw_data.loc[i, :].isnull()
        name = raw_data.loc[i, 'name']
        if not if_null['isotopes']:
            iso_field = raw_data.loc[i, 'isotopes']
            isotope = re.search(r'\[(\d+)\]\[M(\+\d)?\]', iso_field) # has isotope
            if isotope:
                inx = int(isotope.group(1))
                raw_data.loc[i, 'grouped'] = 1
                # if inx not in cluster_collection:
                #     new_cluster = Cluster(inx)
                #     cluster_collection[inx] = new_cluster
                # cluster_collection[inx].add_member(set([name]))
                if isotope.group(2):  # check if has group(2)
                        raw_data.loc[i, 'is_iso'] = 1
                else:  # only has group(1), [M]
                    raw_data.loc[i, 'key_point'] = 1
            else:
                print(raw_data.loc[i, :])
                break
        else:
            inx = 'ng' + str(i)  # not grouped(single) peak cluster index
            raw_data.loc[i, 'key_point'] = 1

        raw_data.loc[i, 'cluster_id'] = inx
        name2cluster_id[name] = inx
        if inx not in cluster_collection:
            cluster_collection[inx] = Cluster(inx, sample_name)
        cluster_collection[inx].add_member(set([name]))
        if raw_data.loc[i, 'is_iso'] == 1:
            cluster_collection[inx].add_isotope(set([name]))
        if raw_data.loc[i, 'key_point'] == 1:
            cluster_collection[inx].set_key_point(name)

    cache['cluster_collection'] = cluster_collection
    cache['name2cluster_id'] = name2cluster_id
    return {'raw_data': raw_data, 'cache': cache}
            

def normalize_rt_ms1(raw_data, rt_tol=10, ms1_tol=0.01, clusters={}):
    """
    get normalized point coordinate by RT and MS1
    :param raw_data: a pandas data frame
    :param rt_tol: RT tolerance, 10s
    :param ms1_tol: MS1 tolerance, 0.01 Da
    :return:
    """
    raw_data.loc[:, 'mzmed'] = raw_data.mzmed / float(ms1_tol)
    raw_data.loc[:, 'predict.rt'] = raw_data.loc[:, 'predict.rt'] / float(rt_tol)
    name2point_coor = {}
    for i in range(len(raw_data)):
        name = raw_data.name[i]
        coor = np.array([raw_data.loc[i, 'predict.rt'], raw_data.loc[i, 'mzmed']])
        name2point_coor[name] = coor
    for c in clusters:  # set cluster central coordinate for the first time
        cluster_central_coor = name2point_coor[clusters[c].key_point]
        clusters[c].set_central_coor(cluster_central_coor)

    # return a dataFrame
    return pd.DataFrame(name2point_coor, index=['rt', 'mz']).T


def get_distance(p_a, p_b):
    """
    Find the Euclidean distances between p_a and p_b
    https://docs.scipy.org/doc/scipy-0.19.0/reference/generated/scipy.spatial.distance.cdist.html
    :param p_a: a numpy array, one point
    :param p_b: a numpy array, can be multiple points
    :return: distance, a scale
    """
    return distance.cdist(p_a, p_b, 'euclidean')


def merge_two_cluster(c1, c2, clusters, all_cluster_id, raw_data):
    """
    merge two different cluster
    :param c1: cluster 1, cluster Class
    :param c2: cluster 2, cluster Class
    :param clusters: all clusters Class
    :param all_cluster_id: a list contains all cluster id
    :param raw_data: raw data, a data frame
    :return: new clusters, new all_cluster_id
    """
    keep_c = c1  # keep this cluster after merge
    del_c = c2  # delete this cluster
    keep_c.add_member(del_c.members)
    keep_c.add_isotope(del_c.isotope)
    keep_c.set_central_coor((keep_c.central_coor + del_c.central_coor)/2.0)
    #TODO: deal with ms2
    peak_name = del_c.key_point
    the_col = raw_data[raw_data.name==peak_name]
    if the_col['grouped'].tolist() == [0]:  # change 'grouped' status in raw_data
        inx = the_col.index.tolist()
        raw_data.loc[inx, 'grouped'] = 1
    del clusters[del_c.id]
    all_cluster_id.remove(del_c.id)


def get_cluster_central_coor(raw_data, clusters):
    """
    get grouped clusters' central coordinate
    :param raw_data: a dataFrame
    :param name2coor: each peak's normalized coordinate(rt, mz), a dataFrame
    :param clusters: all cluster class, cluster_collection;
    :param name2isotope: can check if this peak is an isotope
    :return: grouped cluster to coordinate, a data frame, (cluster<index>, rt, mz)
    """
    g_cluster_id = []  # grouped cluster id
    g_cluster_coor = {'rt': [], 'mz': []}  # grouped cluster coordinate
    # grouped key points' name list
    key_gr_name = list(raw_data[(raw_data.grouped==1) & (raw_data.key_point==1)].loc[:, 'name'])
    for i in clusters:
        c = clusters[i]  # c is a cluster class
        if c.key_point in key_gr_name:
            g_cluster_id.append(i)
            g_cluster_coor['rt'].append(clusters[i].central_coor[0])
            g_cluster_coor['mz'].append(clusters[i].central_coor[1])
    clusters_coor = pd.DataFrame({'cluster_id': g_cluster_id,
                         'rt': g_cluster_coor['rt'], 'mz': g_cluster_coor['mz']})
    clusters_coor = clusters_coor.set_index(['cluster_id'])
    return clusters_coor[['rt', 'mz']]


def compare_clusters_by_coor(grouped_cluster1_coor, cluster2_coor, cluster_collection, name2cluster_id):
    """
    compare two clusters(each may contain lot of points) by coordinate(normalized rt/mz pair)
    and find pair clusters list that distance less than 1
    :param grouped_cluster1_coor: a data frame contains grouped clusters' coordinate
    :param cluster2_coor: another data frame, can same as cluster1_coor
    :return: a list has element of cluster pair, [(cluster_1, cluster_2), (cluster_m, cluster_n), ...]
    """
    cluster_pairs = []
    points_coor_array_1 = grouped_cluster1_coor.values
    points_coor_array_2 = cluster2_coor.values
    points_dis = get_distance(points_coor_array_1, points_coor_array_2)
    if grouped_cluster1_coor is cluster2_coor:  # inner clusters, a square matrix
        inx_l = np.tril_indices(points_dis.shape[0])  # the indices for the lower-triangle of points_dis
        points_dis[inx_l] = 10  # give lower-triangle indices a value bigger than 1
    less_than_one_inx = np.where(points_dis < 1)  # all the points' location that distance is less that 1
    print(less_than_one_inx)
    
    for i in range(len(less_than_one_inx[0])):
        _1 = less_than_one_inx[0][i]
        _2 = less_than_one_inx[1][i]
        #TODO: here has some problem, we may can try cluster id as index for both cluster2
        peak1 = grouped_cluster1_coor.iloc[_1].name   
        peak2 = cluster2_coor.iloc[_2].name
        c1 = cluster_collection[name2cluster_id[peak1]]
        c2 = cluster_collection[name2cluster_id[peak2]]
        # need to compare MS2
        # merge two points that the distance less than 1
        print(peak1, peak2)
        cluster_pairs.append((c1, c2))
        name2cluster_id[peak2] = c1.id
    return {'cluster_pairs': cluster_pairs, 'points_dis': points_dis}

def rt_ms1_compare(raw_data, cache, rt_tol, ms1_tol):
    cluster_collection = cache['cluster_collection']
    name2cluster_id = cache['name2cluster_id']
    all_cluster_id = list(cluster_collection.keys())
    # name2point_coor: the total points' normalized coordinate(rt, mz), a dataFrame
    name2point_coor = normalize_rt_ms1(raw_data.copy(), rt_tol=rt_tol, ms1_tol=ms1_tol, clusters=cluster_collection)
    # all the peaks that were marked as key point and be grouped
    key_points_g = raw_data[(raw_data.grouped==1) & (raw_data.key_point==1)]
    key_points_g_coor = name2point_coor.loc[key_points_g['name']]
    ##  grouped key points dereplication
    # ax = plt.subplot(111)
    # plt.scatter(x=key_points_coor['rt'], y=key_points_coor['mz'])
    # plt.savefig('key_points2.png', dpi=200)
    cluster_pairs_1 = compare_clusters_by_coor(grouped_cluster1_coor=key_points_g_coor,
                                             cluster2_coor=key_points_g_coor,
                                             cluster_collection=cluster_collection,
                                             name2cluster_id=name2cluster_id)
    # merge two points that the distance less than 1
    for (c1, c2) in cluster_pairs_1['cluster_pairs']:
        # We can change cluster_collection and all_cluster_id in function merge_two_cluster directly
        # always merge c2 to c1
        merge_two_cluster(c1, c2, clusters=cluster_collection, all_cluster_id=all_cluster_id, raw_data=raw_data)
    # also need to check the distance among key_points
    # ...
    ## no group point dereplication
    # grouped clusters central coordinate
    no_group_peaks = raw_data[raw_data.grouped==0]  # all no group peaks  
    no_group_peaks_coor = name2point_coor.loc[no_group_peaks['name']]
    max_cluster_num = 2000.0
    segment_num = int(np.ceil(len(no_group_peaks_coor)/max_cluster_num))
    inx_list = np.arange(segment_num, dtype=int)
    inx_list = np.append(inx_list, [-1])
    clusters = cluster_collection
    gc_cen_coor = get_cluster_central_coor(raw_data, cluster_collection)
    all_cluster_pairs = []
    all_points_dis = []
    i = 0
    for i in range(segment_num):
        i_for = inx_list[i] * int(max_cluster_num)
        i_back = inx_list[i+1] * int(max_cluster_num)
        if inx_list[i+1] == -1:
            i_back = inx_list[i+1]
        print(i_for, i_back)
        current_no_group_peaks_coor = no_group_peaks_coor[i_for:i_back]
        print(len(current_no_group_peaks_coor))
        grouped_cluster1_coor=gc_cen_coor
        cluster2_coor=current_no_group_peaks_coor
        cluster_pairs_2 = compare_clusters_by_coor(grouped_cluster1_coor=gc_cen_coor,
                                                 cluster2_coor=current_no_group_peaks_coor,
                                                 cluster_collection=cluster_collection,
                                                 name2cluster_id=name2cluster_id)
        all_cluster_pairs += cluster_pairs_2['cluster_pairs']
        all_points_dis.append(cluster_pairs_2['points_dis'])
    print(all_cluster_pairs)
    # merge two points that the distance less than 1
    if all_cluster_pairs:
        for (c1, c2) in all_cluster_pairs:
            # We can change cluster_collection and all_cluster_id in function merge_two_cluster directly
            # always merge c2 to c1
            merge_two_cluster(c1, c2, clusters=cluster_collection, all_cluster_id=all_cluster_id, raw_data=raw_data)



print('starting...')
raw_data = read_file(root_dir, file_name)
de_isotope_result = de_isotope(raw_data, 'urine1')
raw_data = de_isotope_result['raw_data']
cache = de_isotope_result['cache']
rt_tol = 15
ms1_tol = 0.01
rt_ms1_compare(raw_data, cache, rt_tol=rt_tol, ms1_tol=ms1_tol)

raw_data.to_csv(os.path.join(root_dir, 'result.csv'))








