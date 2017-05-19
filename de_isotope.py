import pandas as pd
import numpy as np
from scipy.spatial import distance
import os
import re
from collections import defaultdict


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


def de_isotope(raw_data):
    raw_data['key_point'] = 0  # 0 means this point isn't a key point
    raw_data['grouped'] = 0  # 0 means this point isn't grouped
    raw_data['is_iso'] = 0
    raw_data['cluster_id'] = 'no_group'
    cluster_collection = defaultdict(list)
    name2cluster_id = {}
    cache = {}
    for i in range(len(raw_data)):
        if_null = raw_data.loc[i, :].isnull()
        if not if_null['isotopes']:
            iso_field = raw_data.loc[i, 'isotopes']
            isotope = re.search(r'\[(\d+)\]\[M(\+\d)?\]', iso_field)
            if isotope:
                inx = isotope.group(1)
                name = raw_data.loc[i, 'name']
                cluster_collection[inx].append(name)
                raw_data.loc[i, 'grouped'] = 1
                raw_data.loc[i, 'cluster_id'] = inx
                name2cluster_id[name] = inx
                if isotope.group(2):  # check if has group(2)
                        raw_data.loc[i, 'is_iso'] = 1
                else:  # only has group(1), [M]
                    raw_data.loc[i, 'key_point'] = 1
    cache['cluster_collection'] = cluster_collection
    cache['name2cluster_id'] = name2cluster_id
    return(raw_data, cache)
            

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
    return name2point_coor


def get_distance(p_a, p_b):
    """
    Find the Euclidean distances between p_a and p_b
    https://docs.scipy.org/doc/scipy-0.19.0/reference/generated/scipy.spatial.distance.cdist.html
    :param p_a: a numpy array, one point
    :param p_b: a numpy array, can be multiple points
    :return: distance, a scale
    """
    return distance.cdist(p_a, p_b, 'euclidean')


def rt_ms1_compare(raw_data, cache, rt_tol, ms1_tol):
    cluster_collection = cache['cluster_collection']
    name2cluster_id = cache['name2cluster_id']
    name2point_coor = normalize_rt_ms1(raw_data, rt_tol=rt_tol, ms1_tol=ms1_tol)
    key_points = raw_data[raw_data.key_point==1]  # all the peaks that were marked as key point
    no_group_peaks = raw_data[raw_data.grouped==0]  # all no group peaks    
    # also need to check the distance among key_points
    # ...
    
    
    
    
raw_data = read_file(root_dir, file_name)
de_isotope_result = de_isotope(raw_data)
raw_data.to_csv(os.path.join(root_dir, 'result.csv'))








