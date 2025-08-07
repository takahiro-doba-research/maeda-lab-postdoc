import os
import re
import time

import networkx as nx
import numpy as np
import pandas as pd
import networkx.algorithms.isomorphism as iso

from .coordinate import Coordinate
from .output import Output
from .util import join_group_data


def read_eq_list(path):

    eq_list = EQList()
    eq_list.name = os.path.splitext(os.path.basename(path))[0]

    with open(path, "r") as f:
        lines = f.readlines()

    eq_indices = [i for i, line in enumerate(lines) if "#" in line]

    if eq_indices:

        for i in range(len(eq_indices)):
            eq = Output()

            if i != len(eq_indices) - 1:
                eq_block = lines[eq_indices[i]:eq_indices[i+1]-1]
            else:
                eq_block = lines[eq_indices[i]:]

            # Set Output.name
            eq.name = int(re.findall(r"EQ(.+),", eq_block[0])[0].strip())
            
            # Set Output.coordinate
            E_index = [
                i for i, line in enumerate(eq_block) if "Energy" in line
            ]
            coordinate = eq_block[1:E_index[0]]
            coordinate = [line.strip().split() for line in coordinate]
            coordinate = [
                [i+1, line[0],
                 float(line[1]), float(line[2]), float(line[3]), []]
                for i, line in enumerate(coordinate)
            ]
            df = pd.DataFrame(
                coordinate,
                columns=["label", "atom", "x", "y", "z", "note"]
            ).set_index("label")
            eq.coordinate = Coordinate(df)

            # Set Output.E
            eq.E = float(
                re.findall(r"\((.+):", eq_block[E_index[0]])[0].strip()
            )

            # Set Output.spin
            spin = [line for line in eq_block if "Spin(**2)" in line]
            eq.spin = float(re.findall(r"=(.+)", spin[0])[0].strip())

            # Set Output.E_zpv
            E_zpv = [line for line in eq_block if "ZPVE" in line]
            eq.E_zpv = float(re.findall(r"=(.+)", E_zpv[0])[0].strip())

            # Set Output.normal_mode
            normal_mode_index = [
                i for i, line in enumerate(eq_block)
                if "Normal mode eigenvalues" in line
            ]
            normal_mode = eq_block[normal_mode_index[0]+1:]
            normal_mode = [
                float(value)
                for line in normal_mode
                for value in line.strip().split()
            ]
            eq.normal_mode = np.array(normal_mode)

            # Add to eq_list
            eq_list[eq.name] = eq

    return eq_list


class EQList:

    def __init__(self, name=None, eqs=None):
        self.name = name
        self.eqs = eqs or {}

    def __len__(self):
        return len(self.eqs)

    def __getitem__(self, key):
        return self.eqs[key]

    def __setitem__(self, key, value):
        self.eqs[key] = value

    def keys(self):
        return self.eqs.keys()

    def values(self):
        return self.eqs.values()

    def items(self):
        return self.eqs.items()

    def show(self, names):

        if isinstance(names, int):
            names = [names]

        for name in names:
            eq = self.eqs[name]
            eq.coordinate.show(f"EQ{name}")
            time.sleep(2)

    def shorter_distance(self, index0, index1, distance):
        eqs_shorter_distance = [
            name for name, eq in self.items()
            if eq.coordinate.get_distance(index0, index1) < distance
        ]
        return eqs_shorter_distance

    def longer_distance(self, index0, index1, distance):
        eqs_longer_distance = [
            name for name, eq in self.items()
            if distance < eq.coordinate.get_distance(index0, index1)
        ]
        return eqs_longer_distance

    def between_distance(self, index0, index1, distance0, distance1):
        eqs_between_distance = [
            name for name, eq in self.items()
            if distance0 < eq.coordinate.get_distance(index0, index1) < distance1
        ]
        return eqs_between_distance

    def between_dihedral_angle(
            self, index0, index1, index2, index3,
            dihedral_angle0, dihedral_angle1):
        eqs_between_dihedral_angle = [
            name for name, eq in self.items()
            if dihedral_angle0 < eq.coordinate.get_dihedral_angle(
                index0, index1, index2, index3) < dihedral_angle1
        ]
        return eqs_between_dihedral_angle

    def adjust_energy(self, name=None):

        if name == None:
            eq_std = min(self.values(), key=lambda eq:eq.E)
            E_std = eq_std.E

        else:
            E_std = self[name].E

        for eq in self.values():
            eq.E = (eq.E - E_std) * 2625.5

    def save_eqs(self, path):
        os.makedirs(f"{path}/EQ", exist_ok=True)
        header = ["# STO-3G\n", "\n", "title\n", "\n", "0 1\n"]
        footer = ["\n"]

        for name, eq in self.items():
            eq.coordinate.save(
                path=f"{path}/EQ/{name}.com",
                header=header,
                footer=footer,
                ignore_notes=True
            )

    def set_group(self):
        adj_matrices = {}
        group_count = 0

        for name, eq in self.items():
            matrix = eq.coordinate.get_adj_matrix()

            for group, adj_matrix in adj_matrices.items():

                if np.all(matrix == adj_matrix):
                    eq.group = group
                    break

            else:
                adj_matrices[group_count] = matrix
                eq.group = group_count
                group_count += 1

    def set_group_isomorphic(self):
        graphs = {}
        group_count = 0
        nm = iso.categorical_node_match("atom", None)

        for name, eq in self.items():
            graph = eq.coordinate.to_graph()

            for group, g in graphs.items():

                if nx.is_isomorphic(graph, g, node_match=nm):
                    eq.group = group
                    break

            else:
                graphs[group_count] = graph
                eq.group = group_count
                group_count += 1

    def group(self, path):
        folder = f"{path}/{self.name}_group"
        os.makedirs(folder, exist_ok=True)
        header = ["# STO-3G\n", "\n", "title\n", "\n", "0 1\n"]
        footer = ["\n"]
        eq_data = []
        group_data = []
        group_count = 0

        for name, eq in self.items():
            labels = eq.coordinate.label
            adj_matrix = eq.coordinate.get_adj_matrix()

            for group_dict in group_data:

                if np.all(adj_matrix == group_dict["adj_matrix"]):
                    eq_data.append({"EQ":name, "group":group_dict["group"]})
                    eq.coordinate.save(
                        f"{folder}/{group_dict['group']}/{name}.com",
                        header=header,
                        footer=footer,
                        ignore_notes=True
                    )
                    break

            else:
                eq_data.append({"EQ":name, "group":group_count})
                group_data.append(
                    {"group":group_count, "adj_matrix":adj_matrix}
                )
                os.makedirs(f"{folder}/{group_count}", exist_ok=True)
                eq.coordinate.save(
                    f"{folder}/{group_count}/{name}.com",
                    header=header,
                    footer=footer,
                    ignore_notes=True
                )
                group_count += 1

        eq_df = pd.DataFrame(eq_data)
        eq_df.to_csv(
            f"{folder}/EQ.csv",
            columns=["EQ", "group"],
            index=False
        )
        group_data = join_group_data(group_data)
        group_df = pd.DataFrame(group_data)
        group_df.to_csv(
            f"{folder}/group.csv",
            columns=["group", "adj_matrix"],
            index=False
        )
