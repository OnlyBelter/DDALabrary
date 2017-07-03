import os
import re
try:
    from cluster import MS2Info
    from public_helper import get_mz_pair
except:
    from .cluster import MS2Info
    from .public_helper import get_mz_pair
import json
import numpy as np
import pandas as pd
from numpy.linalg import norm



root_d = r'D:\vm_share_folder\DDA_lib\02MS2_match'
ms2_spec_file = 'urine1_spectra.msp'
sample_name = 'urine1'

def read_msp(root_d, ms2_f, sample_name):
    """
    read msp file which contains ms2 information
    :param root_d:
    :param ms2_f:
    :param sample_name: unique with each other
    :return: {peak_name: MS2Info Class, ..., peak_name: MS2Info Class}
    """
    all_ms2 = {}
    with open(os.path.join(root_d, ms2_f), 'r') as f_handle:
        inx = ''
        for each_line in f_handle:
            each_line = each_line.strip()
            search_inx = re.search(r'NAME:\s+(M\d+T\d+)', each_line)
            if search_inx:
                inx = search_inx.group(1)
                all_ms2[inx] = MS2Info(inx, sample_name)
            if re.search(r'^\d+', each_line):
                mz, intensity = each_line.split(' ')
                all_ms2[inx].add_mz(mz)
                all_ms2[inx].add_int(intensity)
    return all_ms2


def print_ms2_info(root_d, ms2_dic, file_n):
    """
    print ms2 info as json file
    :param root_d:
    :param ms2_dic: {'peak_name': MS2Info class, ..., }
    :param file_n: file name
    :return:
    """
    ms2_info_dic = {}
    ms2_info = ms2_dic
    for i in ms2_info:
        ms2_info_dic[i] = {'mz': [], 'int': [], 'max_int': 0, 'segment_num': 0}
        ms2_info_dic[i]['mz'] = ms2_info[i].get_mz()
        ms2_info_dic[i]['int'] = ms2_info[i].get_int()
        ms2_info_dic[i]['max_int'] = ms2_info[i].get_max_int()
        ms2_info_dic[i]['segment_num'] = ms2_info[i].get_segment_num()
        ms2_info_dic[i]['sample_name'] = ms2_info[i].sample_name
    with open(os.path.join(root_d, file_n), 'a') as f_handle:
        json.dump(ms2_info_dic, f_handle, indent=2)


def ms2_pre_process(root_d, ms2_f, output_f, sample_name, tol_mz=0.01):

    print('Start to read msp file...')
    ms2_info_all = read_msp(root_d, ms2_f, sample_name)
    # print_ms2_info(root_d, ms2_info_all.copy(), file_n=sample_name + '_total_ms2.json')

    print('Start to remove noise...')
    # remove noise, intensity lower than 1% of the highest intensity in this peak removed
    for ms2 in ms2_info_all:
        ms2_info = ms2_info_all[ms2]
        int_array = np.array(ms2_info.get_int())  # ms2 intensity
        int_percentage = int_array/float(ms2_info.get_max_int())
        need_remove_inx = np.where(int_percentage < 0.01)
        ms2_info.remove_by_inx(need_remove_inx)
    # print_ms2_info(root_d, ms2_info_all.copy(), file_n=sample_name + '_after_remove_noise.json')

    print('Start to remove ring effect...')
    # remove ring effect, merge fragments within the same m/z bin across spectra
    # and the removed one's intensity is less than 20% of the highest one's in this bin
    for ms2 in ms2_info_all:
        ms2_info = ms2_info_all[ms2]
        int_array = np.array(ms2_info.get_int()) # intensity
        mz_array = np.array(ms2_info.get_mz())
        mz_pairs = get_mz_pair(mz_array, tol_mz=0.01) # m/z pair in the same m/z bin
        need_remove_inx = []
        if mz_pairs:
            for pair in mz_pairs:
                small_int_inx = 0
                int_pair = int_array[pair['inx']]
                if int_pair[1] < int_pair[0]:
                    small_int_inx = 1
                int_percentage = int_pair/float(np.max(int_pair))
                # print(int_percentage)
                if int_percentage[small_int_inx] <= .2:
                    # print(int_percentage.shape)
                    # mz_array[pair['inx'][small_int_inx]] = 0
                    need_remove_inx.append(pair['inx'][small_int_inx])
        if need_remove_inx:
            # print(ms2)
            # print(mz_array)
            # print(need_remove_inx)
            ms2_info.remove_by_inx(need_remove_inx)
    # print_ms2_info(root_d, ms2_info_all.copy(), file_n=sample_name + '_after_remove_ring_effect2.json')


def dop_product(lib_peak, exp_peak, match_method):
    """
    this function can calculate two peaks' dot product in two ways - forward and reverse
    :param lib_peak: a pandas data frame, contains mz-intensity list
    :param exp_peak:
    :param match_method: forward or reverse
    :return: dot product, a scalar
    """
    # print(lib_peak.mz**3)
    lib_peak_vector = np.multiply(np.array(lib_peak.mz)**3, np.array(lib_peak.int)**0.6)
    exp_peak_vector = np.multiply(np.array(exp_peak.mz)**3, np.array(exp_peak.int)**0.6)
    inx = np.empty((0))
    if match_method == 'forward':
        inx = np.append(inx, np.where(exp_peak_vector != 0))
    elif match_method == 'reverse':
        inx = np.append(inx, np.where(lib_peak_vector != 0))
    else:
        print('Please input correct match method, only \'forward\' or \'reverse\' can be accepted')
    # print('inx is', inx.astype(int))
    v1 = lib_peak_vector[inx.astype(int)]
    v2 = exp_peak_vector[inx.astype(int)]
    # print(v1)
    # print(v2)
    return np.dot(v1, v2)/(norm(v1)*norm(v2))
    # print(type(lib_peak_vector))
    # print(type(exp_peak_vector))


def ms2_compare(lib_peak_name, exp_peak_name, ms2_info_all, compare_method):
    """
    this function is used to compare two different peaks, if they have similar ms1 and rt
    :param lib_peak_name:
    :param exp_peak_name:
    :param ms2_info_all: MS2 information, a json file
    :param compare_method: forward or reverse
    :return:
    """
    with open(ms2_info_all, 'r') as f_handle:
        ms2_info = json.load(f_handle)
    lib_peak = ms2_info[lib_peak_name]
    exp_peak = ms2_info[exp_peak_name]
    # get 'mz' to 'int' key-value pair
    lib_peak_dic = {lib_peak['mz'][i]: lib_peak['int'][i] for i in range(len(lib_peak['mz']))}
    exp_peak_dic = {exp_peak['mz'][i]: exp_peak['int'][i] for i in range(len(exp_peak['mz']))}
    mz_union_list = np.unique(np.array(list(lib_peak_dic.keys()) + list(exp_peak_dic.keys())))
    peak_union_df = pd.DataFrame(data=[mz_union_list], index=['mz']).T
    lib_peak_union = peak_union_df.copy()
    exp_peak_union = peak_union_df.copy()
    lib_peak_union['int'] = [lib_peak_dic.get(mz, 0) for mz in mz_union_list]
    exp_peak_union['int'] = [exp_peak_dic.get(mz, 0) for mz in mz_union_list]
    # print(mz_union_list)
    # print(lib_peak_union)
    # print(exp_peak_union)
    # get mz pairs that need to merge
    mz_pairs = get_mz_pair(mz_union_list, tol_mz=0.01)
    # print(mz_pairs)
    need_remove_inx = []
    for pair in mz_pairs:
        # remain bigger intensity and remove smaller mz in a pair
        inx = pair['inx'] # a list contains two elements, [0, 1]
        p1_pair_int = lib_peak_union.loc[inx].int
        p2_pair_int = exp_peak_union.loc[inx].int
        lib_peak_union.loc[inx, 'int'] = max(p1_pair_int)
        exp_peak_union.loc[inx, 'int'] = max(p2_pair_int)
        need_remove_inx.append(inx[0])
    lib_peak_union = lib_peak_union.drop(need_remove_inx)
    exp_peak_union = exp_peak_union.drop(need_remove_inx)
    # print(lib_peak_union)
    # print(exp_peak_union)
    # print(need_remove_inx)
    # need to calculate dot product between these two peaks
    dp = dop_product(lib_peak=lib_peak_union, exp_peak=exp_peak_union, match_method='forward')
    # dp = dop_product(lib_peak=lib_peak_union, exp_peak=exp_peak_union, match_method='reverse')
    print('dp is', dp)



# ms2_pre_process(root_d, ms2_spec_file, output_f='after_processed_ms2.json', sample_name=sample_name)
ms2_info_file = os.path.join(root_d, 'urine1_after_remove_ring_effect2.json')
ms2_compare('M137T542', 'M137T555', ms2_info_file, 'forward')
