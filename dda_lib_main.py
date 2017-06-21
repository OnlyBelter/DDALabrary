import os
import pandas as pd
import argparse
import json
import numpy as np
# import matplotlib
# import matplotlib.pyplot as plt

from code.de_isotope import de_isotope
from code.ms1_compare import rt_ms1_compare


# root_dir = r'D:\vm_share_folder\DDA_lib\01RT_QC\result'
# file_name = 'corrected_pos_1_0.csv'

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


class ClusterCollection():

    def create_new_cluster_id(self):
        pass


def main(params):
    print('start to read MetAnalyzer result...')
    root_dir = params['file_path']
    raw_data = read_file(root_dir, params['file_name'])
    print('start to de-isotope...')
    de_isotope_result = de_isotope(raw_data, 'urine1')
    raw_data = de_isotope_result['raw_data']
    cluster_collection = de_isotope_result['cluster_collection']
    rt_tol = 15
    ms1_tol = 0.01
    rt_ms1_compare_result = rt_ms1_compare(raw_data, cluster_collection, rt_tol=rt_tol, ms1_tol=ms1_tol)
    raw_data2 = rt_ms1_compare_result['raw_data']
    raw_data2.to_csv(os.path.join(root_dir, 'ms1_compare_result.csv'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file_path', type=str,
                        help='the result file comes from MetAnalyser and need to correct RT')
    parser.add_argument('file_name', default='result.csv', type=str, help='file name')
    args = parser.parse_args([r'D:\vm_share_folder\DDA_lib\01RT_QC\result', 'corrected_pos_1.csv'])
    params = vars(args)  # convert to ordinary dict
    print('parsed parameters: ')
    print(json.dumps(params, indent=2))
    main(params)








