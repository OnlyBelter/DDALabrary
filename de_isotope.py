import pandas as pd
import numpy as np
# import matplotlib
# import matplotlib.pyplot as plt
from scipy.spatial import distance
import os
import re

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
    """
    de-isotope and create cluster_collection which contains all peaks
    :param raw_data:
    :param sample_name: each sample should has an unique name
    :return:
    """
    raw_data['key_point'] = 0  # 0 means this point isn't a key point
    raw_data['grouped'] = 0  # 0 means this point isn't grouped
    raw_data['is_iso'] = 0
    raw_data['cluster_id'] = 'no_group'
    cluster_collection = {}  # cluster id to class Cluster
    # name2cluster_id = {}  # whole peak name with cluster id contains isotope
    # cache = {}
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
        # name2cluster_id[name] = inx
        if inx not in cluster_collection:
            cluster_collection[inx] = Cluster(inx, sample_name)
        cluster_collection[inx].add_member(set([name]))
        if raw_data.loc[i, 'is_iso'] == 1:
            cluster_collection[inx].add_isotope(set([name]))
        if raw_data.loc[i, 'key_point'] == 1:
            cluster_collection[inx].set_key_point(name)
    raw_data = raw_data.set_index(['name'])
    # cache['cluster_collection'] = cluster_collection
    # cache['name2cluster_id'] = name2cluster_id
    return {'raw_data': raw_data, 'cluster_collection': cluster_collection}


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


def get_distance(p_a, p_b):
    """
    Find the Euclidean distances between p_a and p_b
    https://docs.scipy.org/doc/scipy-0.19.0/reference/generated/scipy.spatial.distance.cdist.html
    :param p_a: a numpy array, one point
    :param p_b: a numpy array, can be multiple points
    :return: distance, a scale
    """
    return distance.cdist(p_a, p_b, 'euclidean')


def merge_cluster_pairs(cluster_pairs, clusters, raw_data, cluster2point_coor):
    """
    We can change clusters and raw_data in this function directly
    merge two different clusters that the distance less than 1
    :param cluster_pairs: a list contains cluster pairs, eg [(c1, c2), (c4, c5), (c6, c7), ...]
                          c1: cluster 1, a cluster Class
                          c2: cluster 2, a cluster Class
                          always merge c2 to c1
    :param all_cluster_id: a list contains all cluster id
    :param raw_data: raw data, a data frame
    :return: update clusters and raw_data
    eg: cluster ng0, ng1, ng2, ng3, ng4 is one cluster
    The distance between  ng0  and  ng1  is  0.0674124947206
    The distance between  ng0  and  ng2  is  0.133845686437
    The distance between  ng0  and  ng3  is  0.200344858931
    The distance between  ng0  and  ng4  is  0.266925584437
    The distance between  ng1  and  ng2  is  0.0666883381443
    The distance between  ng1  and  ng3  is  0.133344817214
    The distance between  ng1  and  ng4  is  0.200007691142
    The distance between  ng2  and  ng3  is  0.0666666854167
    The distance between  ng2  and  ng4  is  0.133333344268
    The distance between  ng3  and  ng4  is  0.0666666667867

    """
    for (c1, c2) in cluster_pairs:
        # (c1, c2) = cluster_pairs[0]
        all_cluster_id = cluster2point_coor.index
        if (c1.id in all_cluster_id) and (c2.id in all_cluster_id):
            keep_c = c1  # keep this cluster after merge
            del_c = c2  # delete this cluster
            del_c_members = del_c.members
            keep_c.add_member(del_c_members)
            keep_c.add_isotope(del_c.isotope)
            keep_c.set_central_coor((keep_c.central_coor + del_c.central_coor)/2.0)
            # update coordinate in cluster2point_coor
            cluster2point_coor.loc[keep_c.id] = keep_c.central_coor
            #TODO: deal with ms2
            raw_data.loc[list(del_c_members), 'cluster_id'] = keep_c.id
            raw_data.loc[list(del_c_members), 'key_point'] = 0
            print('merge ', del_c.id, 'to ', keep_c.id)
            if del_c.id in clusters: del clusters[del_c.id]
            if del_c.id in cluster2point_coor.index:
                cluster2point_coor.drop([del_c.id], inplace=True)


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


def compare_clusters_by_coor(grouped_cluster1_coor, cluster2_coor, cluster_collection):
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
    # print(less_than_one_inx)
    
    for i in range(len(less_than_one_inx[0])):
        _1 = less_than_one_inx[0][i]
        _2 = less_than_one_inx[1][i]
        # cluster id as index
        cluster_id1 = grouped_cluster1_coor.iloc[_1].name   
        cluster_id2 = cluster2_coor.iloc[_2].name
        c1 = cluster_collection[cluster_id1]
        c2 = cluster_collection[cluster_id2]
        # need to compare MS2
        # merge two points that the distance less than 1
        print('The distance between ', cluster_id1, ' and ', cluster_id2, ' is ', points_dis[_1, _2])
        cluster_pairs.append((c1, c2))
        # name2cluster_id[peak2] = c1.id
    return {'cluster_pairs': cluster_pairs}

def rt_ms1_compare(raw_data, cluster_collection, rt_tol, ms1_tol):
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


print('starting...')
raw_data = read_file(root_dir, file_name)
de_isotope_result = de_isotope(raw_data, 'urine1')
raw_data = de_isotope_result['raw_data']
cluster_collection = de_isotope_result['cluster_collection']
rt_tol = 15
ms1_tol = 0.01
rt_ms1_compare_result = rt_ms1_compare(raw_data, cluster_collection, rt_tol=rt_tol, ms1_tol=ms1_tol)
raw_data2 = rt_ms1_compare_result['raw_data']
raw_data2.to_csv(os.path.join(root_dir, 'result3.csv'))








