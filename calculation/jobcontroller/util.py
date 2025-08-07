import copy
import itertools

import numpy as np
import pandas as pd


def join_separated_group_info(separated_group_info):
    """
    Converts
    [{"group":0, "labels":[1, 2], "adj_matrix":np.array([[0, 1],[1, 0]])},
     {"group":0, "labels":[3, 4], "adj_matrix":np.array([[0, 1],[1, 0]])},
     ...]
    to
    [{"group":0, "labels":1_2, "adj_matrix":0_1_1_0},
     {"group":0, "labels":3_4, "adj_matrix":0_1_1_0},
     ...]
    """
    separated_group_info = copy.deepcopy(separated_group_info)

    for line in separated_group_info:
        line["labels"] = "_".join(
            [str(value) for value in line["labels"]]
        )
        line["adj_matrix"] = "_".join(
            [str(value) for value in line["adj_matrix"].flatten()]
        )

    return separated_group_info


def join_group_info(group_info):
    """
    Converts
    [{"group":0, "adj_matrix":np.array([[0, 1],[1, 0]])},
     {"group":0, "adj_matrix":np.array([[0, 1],[1, 0]])},
     ...]
    to
    [{"group":0, "adj_matrix":0_1_1_0},
     {"group":0, "adj_matrix":0_1_1_0},
     ...]
    """
    group_info = copy.deepcopy(group_info)

    for line in group_info:
        line["adj_matrix"] = "_".join(
            [str(value) for value in line["adj_matrix"].flatten()]
        )

    return group_info


def split_group_info(group_info):
    """
    The opposite transformation of join_group_info.
    """
    group_info = copy.deepcopy(group_info)

    for line in group_info:
        splitted_adj_matrix = line["adj_matrix"].split("_")
        l = int(np.sqrt(len(splitted_adj_matrix)))
        line["adj_matrix"] = np.array(
            [int(value) for value in splitted_adj_matrix]
        ).reshape(l, l)

    return group_info


def split_group_info_df(group_info_df):
    group_info = group_info_df.to_dict(orient="records")
    group_info = split_group_info(group_info)
    return pd.DataFrame(group_info)


def split_separated_group_info(separated_group_info):
    """
    The opposite transformation of join_separated_group_info.
    """
    separated_group_info = copy.deepcopy(separated_group_info)

    for line in separated_group_info:
        line["labels"] = [
            int(value) for value in line["labels"].split("_")
        ]
        l = len(line["labels"])
        line["adj_matrix"] = np.array(
            [int(value) for value in line["adj_matrix"].split("_")]
        ).reshape(l, l)

    return separated_group_info


def split_separated_group_info_df(separated_group_info_df):
    separated_group_info = separated_group_info_df.to_dict(orient="records")
    separated_group_info = split_separated_group_info(separated_group_info)
    return pd.DataFrame(separated_group_info)


def dict_product(dictionary):
    """
    Converts
    dictionary = {"a":[0, 1], "b":[0, 1, 2]}
    to
    iter([
        {"a":0, "b":0},
        {"a":0, "b":1},
        {"a":0, "b":2},
        {"a":1, "b":0},
        {"a":1, "b":1},
        {"a":1, "b":2}
    ])
    """
    keys = dictionary.keys()
    values_list = dictionary.values()

    for values in itertools.product(*values_list):
        product = dict(zip(keys, values))
        yield product


def dihedral_angles_product(lines):
    """
    Converts
    lines = [
        (0, 1, {'dihedral_angles': [0]}),
        (0, 2, {'dihedral_angles': [0]}),
        (0, 3, {'dihedral_angles': [0, 180]}),
        (0, 4, {'dihedral_angles': [0]}),
        (0, 5, {'dihedral_angles': [0, 180]})
    ]
    to
    iter([
        [
            (0, 1, {'dihedral_angle': 0}),
            (0, 2, {'dihedral_angle': 0}),
            (0, 3, {'dihedral_angle': 0}),
            (0, 4, {'dihedral_angle': 0}),
            (0, 5, {'dihedral_angle': 0})
        ],
        [
            (0, 1, {'dihedral_angle': 0}),
            (0, 2, {'dihedral_angle': 0}),
            (0, 3, {'dihedral_angle': 0}),
            (0, 4, {'dihedral_angle': 0}),
            (0, 5, {'dihedral_angle': 180})
        ],
        [
            (0, 1, {'dihedral_angle': 0}),
            (0, 2, {'dihedral_angle': 0}),
            (0, 3, {'dihedral_angle': 180}),
            (0, 4, {'dihedral_angle': 0}),
            (0, 5, {'dihedral_angle': 0})
        ],
        [
            (0, 1, {'dihedral_angle': 0}),
            (0, 2, {'dihedral_angle': 0}),
            (0, 3, {'dihedral_angle': 180}),
            (0, 4, {'dihedral_angle': 0}),
            (0, 5, {'dihedral_angle': 180})
        ]
    ])
    """
    dihedral_angles_list = [line[2]["dihedral_angles"] for line in lines]

    for dihedral_angles in itertools.product(*dihedral_angles_list):
        product = [
            (line[0], line[1], {"dihedral_angle":dihedral_angle})
            for line, dihedral_angle in zip(lines, dihedral_angles)
        ]
        yield product

def get_key_from_value(d, val):
    keys = [k for k, v in d.items() if v == val]

    if keys:
        return keys[0]
    else:
        return None
