from scipy.spatial import distance


def get_distance(p_a, p_b):
    """
    Find the Euclidean distances between p_a and p_b
    https://docs.scipy.org/doc/scipy-0.19.0/reference/generated/scipy.spatial.distance.cdist.html
    :param p_a: a numpy array, one point
    :param p_b: a numpy array, can be multiple points
    :return: distance, a scale
    """
    return distance.cdist(p_a, p_b, 'euclidean')
