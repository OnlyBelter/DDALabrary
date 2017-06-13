import re
from .cluster import Cluster


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
