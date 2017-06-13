from scipy.spatial import distance
import numpy as np


def get_distance(p_a, p_b):
    """
    Find the Euclidean distances between p_a and p_b
    https://docs.scipy.org/doc/scipy-0.19.0/reference/generated/scipy.spatial.distance.cdist.html
    :param p_a: a numpy array, one point or multiple points
    :param p_b: a numpy array, one point or multiple points
    :return: distance, a scale
    """
    print(p_a.shape, p_b.shape)
    p_b_lines_num = p_b.shape[0]
    part_num = 1
    if p_b_lines_num > 1500:
        part_num = divmod(p_b_lines_num, 1500)[0]
    p_b_list = np.array_split(p_b, part_num)
    print(len(p_b_list))
    result = []
    for i in p_b_list:
        print(i.shape)
        dis = distance.cdist(p_a, i, 'euclidean')
        result.append(dis)
    return np.hstack(result)
