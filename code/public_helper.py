from scipy.spatial import distance
import numpy as np


def get_distance(p_a, p_b):
    """
    Python 3 needs smaller memory
    Find the Euclidean distances between p_a and p_b
    https://docs.scipy.org/doc/scipy-0.19.0/reference/generated/scipy.spatial.distance.cdist.html
    :param p_a: a numpy array, one point or multiple points
    :param p_b: a numpy array, one point or multiple points
    :return: distance among each point, a matrix
    """
    # print(p_a.shape, p_b.shape)
    p_a_lines_num = p_a.shape[0]
    part_num = 1
    if p_a_lines_num > 2048:
        part_num = divmod(p_a_lines_num, 2048)[0]
    p_a_list = np.array_split(p_a, part_num)
    # print(len(p_a_list))
    result = []
    for i in p_a_list:
        # print(i.shape)
        dis = np.float32(distance.cdist(i, p_b, 'euclidean'))
        result.append(dis)
    return np.vstack(result)


def get_mz_pair(mz_list, tol_mz):
    """
    this function can calculate m/z distance in MS2 for each peak
    :param mz_list: a numpy array, m/z list
    :return: same m/z bin ms pairs
    """
    new_ms_list = np.zeros([len(mz_list), 2])
    new_ms_list[:, 0] = mz_list
    ms_distance = get_distance(new_ms_list, new_ms_list)
    inx_l = np.tril_indices(ms_distance.shape[0])  # the indices for the lower-triangle of points_dis
    ms_distance[inx_l] = tol_mz + 1  # give lower-triangle indices a value bigger than tol_mz
    less_than_one_inx = np.where(ms_distance <= tol_mz)  # all the points' location that distance is less that tol_mz
    mz_pair = []
    if less_than_one_inx:
        _ltoi = less_than_one_inx
        half_pair_a = mz_list[_ltoi[0]]
        half_pair_b = mz_list[_ltoi[1]]
        for i in np.arange(len(half_pair_a)):
            mz_pair.append({'inx': [_ltoi[0][i], _ltoi[1][i]], 'mz': [round(half_pair_a[i], 4), round(half_pair_b[i], 4)]})
    return mz_pair

