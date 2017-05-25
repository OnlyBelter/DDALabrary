import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.spatial import distance
import os
import re
from collections import defaultdict

matplotlib.style.use('ggplot')
root_dir = r'D:\vm_share_folder\DDA_lib\01RT_QC\result'
file_name = 'corrected_pos_1.csv'

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
    def __init__(self, cluster_id):
        self.id = cluster_id
        self.members = set()  # peak name
        self.isotope = set()
        self.key_point = ''  # has min mz
        self.central_coor = []
        self.ms2_num = 0
        self.sample_name = ''  # sample name

    def add_member(self, new_member):
        # new_member is a set
        self.members = self.members.union(new_member)

    def add_isotope(self, isotope):
        # isotope is a set
        self.isotope = self.isotope.union(isotope)

    def set_key_point(self, key_point):
        self.key_point = key_point

    def set_central_coor(self, coor):
        # coor is a list has two elements
        self.set_central_coor = coor

    def set_sample_name(self, sample_name):
        # sample name is a string
        self.sample_name = sample_name

    def set_ms2_num(self, num):
        self.ms2_num = num

    def get_members(self):
        return self.members

    def get_uniso_mem(self): # un isotope member
        return self.members - self.isotope

    def get_central_coor(self):
        return self.central_coor


class ClusterCollection():

    def create_new_cluster_id(self):
        pass


def de_isotope(raw_data):
    raw_data['key_point'] = 0  # 0 means this point isn't a key point
    raw_data['grouped'] = 0  # 0 means this point isn't grouped
    raw_data['is_iso'] = 0
    raw_data['cluster_id'] = 'no_group'
    cluster_collection = {}  # cluster id to class Cluster
    name2cluster_id = {}
    cache = {}
    for i in range(len(raw_data)):
        if_null = raw_data.loc[i, :].isnull()
        if not if_null['isotopes']:
            iso_field = raw_data.loc[i, 'isotopes']
            isotope = re.search(r'\[(\d+)\]\[M(\+\d)?\]', iso_field) # has isotope
            if isotope:
                inx = int(isotope.group(1))
                name = raw_data.loc[i, 'name']
                if inx not in cluster_collection:
                    new_cluster = Cluster(inx)
                    cluster_collection[inx] = new_cluster
                cluster_collection[inx].add_member(set([name]))
                raw_data.loc[i, 'grouped'] = 1
                raw_data.loc[i, 'cluster_id'] = inx
                name2cluster_id[name] = inx
                if isotope.group(2):  # check if has group(2)
                        raw_data.loc[i, 'is_iso'] = 1
                        cluster_collection[inx].add_isotope(set([name]))
                else:  # only has group(1), [M]
                    raw_data.loc[i, 'key_point'] = 1
                    cluster_collection[inx].set_key_point(name)
    cache['cluster_collection'] = cluster_collection
    cache['name2cluster_id'] = name2cluster_id
    return {'raw_data': raw_data, 'cache': cache}
            

def normalize_rt_ms1(raw_data, rt_tol=10, ms1_tol=0.01):
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


def merge_two_points():
    pass


def get_cluster_central(name2coor, name2cluster_id, cluster_collection, name2isotope, cluster2coor={}):
    """
    calculate a cluster's central coordinate without isotope
    :param name2coor: each peak's normalized coordinate(rt, mz), a dataFrame
    :param name2cluster_id:
    :param cluster_collection:
    :param name2isotope: can check if this peak is an isotope
    :return: cluster to coordinate, each cluster has one central coordinate
    """
    if not cluster2coor:  # for the first time to get cluster2coor
        for (i, names) in cluster_collection:
            pass
        



def rt_ms1_compare(raw_data, cache, rt_tol, ms1_tol):
    cluster_collection = cache['cluster_collection']
    name2cluster_id = cache['name2cluster_id']
    # name2point_coor: the total points' normalized coordinate(rt, mz), a dataFrame
    name2point_coor = normalize_rt_ms1(raw_data.copy(), rt_tol=rt_tol, ms1_tol=ms1_tol)
    key_points = raw_data[raw_data.key_point==1]  # all the peaks that were marked as key point
    no_group_peaks = raw_data[raw_data.grouped==0]  # all no group peaks    
    key_points_coor = name2point_coor.loc[key_points['name']]
    no_group_peaks_coor = name2point_coor.loc[key_points['name']]
    # ax = plt.subplot(111)
    
    # key points dereplication
    plt.scatter(x=key_points_coor['rt'], y=key_points_coor['mz'])
    plt.savefig('key_points2.png', dpi=200)
    key_points_coor_array = key_points_coor.values
    key_points_inner_dis = get_distance(key_points_coor_array, key_points_coor_array)
    inx_l = np.tril_indices(key_points_inner_dis.shape[0])  # the indices for the lower-triangle of key_points_inner_dis
    key_points_inner_dis[inx_l] = 10  # give lower-triangle indices a value bigger than 1
    less_than_one_inx = np.where(key_points_inner_dis < 1)  # all the points' location that distance is less that 1
    plt.scatter(x=less_than_one_inx[0], y= less_than_one_inx[1])
    name2isotope = raw_data.set_index('name')['is_iso'].to_dict()
    for i in range(len(less_than_one_inx[0])):
        _1 = less_than_one_inx[0][i]
        _2 = less_than_one_inx[1][i]
        # need to compare MS2
        # merge two points that the distance less than 1
        print(key_points_coor.iloc[_1].name, key_points_coor.iloc[_2].name)
    # also need to check the distance among key_points
    # ...
    
    
    
    
raw_data = read_file(root_dir, file_name)
de_isotope_result = de_isotope(raw_data)
raw_data = de_isotope_result['raw_data']
cache = de_isotope_result['cache']
rt_tol = 15
ms1_tol = 0.01

raw_data.to_csv(os.path.join(root_dir, 'result.csv'))








