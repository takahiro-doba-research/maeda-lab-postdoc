import copy
import itertools

import numpy as np
import pandas as pd

from .data import covalent_radii


def join_group_data(group_data):
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
    group_data = copy.deepcopy(group_data)

    for group_dict in group_data:
        group_dict["adj_matrix"] = "_".join(
            [str(value) for value in group_dict["adj_matrix"].flatten()]
        )

    return group_data

def R_between(atom0, atom1):
    R = covalent_radii(atom0) + covalent_radii(atom1)
    return R
