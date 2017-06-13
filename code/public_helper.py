from scipy.spatial import distance
import numpy as np


def get_distance(p_a, p_b):
    """
    Python 3 needs smaller memory
    Find the Euclidean distances between p_a and p_b
    https://docs.scipy.org/doc/scipy-0.19.0/reference/generated/scipy.spatial.distance.cdist.html
    :param p_a: a numpy array, one point or multiple points
    :param p_b: a numpy array, one point or multiple points
    :return: distance, a scale
    """
    print(p_a.shape, p_b.shape)
    p_a_lines_num = p_a.shape[0]
    part_num = 1
    if p_a_lines_num > 2048:
        part_num = divmod(p_a_lines_num, 2048)[0]
    p_a_list = np.array_split(p_a, part_num)
    print(len(p_a_list))
    result = []
    for i in p_a_list:
        print(i.shape)
        dis = distance.cdist(i, p_b, 'euclidean')
        result.append(dis)
    return np.vstack(result)
