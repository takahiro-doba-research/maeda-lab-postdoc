import copy
import subprocess

import networkx as nx
import numpy as np
import pandas as pd
import networkx.algorithms.isomorphism as iso
from scipy.spatial import distance

import calculationfile.affine_transformation as at
import calculationfile.data as data


class Coordinate:

    def __init__(self, df=None):
        self.df = df

        if not isinstance(self.df, pd.core.frame.DataFrame):
            self.df = pd.DataFrame(
                columns=["label", "atom", "x", "y", "z", "note"]
            )

    def __len__(self):
        return len(self.df)

    @property
    def label(self):
        return list(self.df.index)

    @property
    def atom(self):
        return list(self.df["atom"])

    @property
    def note(self):
        return list(self.df["note"].copy(deep=True))

    @label.setter
    def label(self, labels):
        self.df.index = labels

    @atom.setter
    def atom(self, atoms):
        self.df["atom"] = atoms

    @note.setter
    def note(self, notes):
        self.df["note"] = notes

    def aslabel(self, coordinate):
        self.df.index = list(coordinate.df.index)

    def asatom(self, coordinate):
        self.df["atom"] = list(coordinate.df["atom"])

    def asnote(self, coordinate):
        self.df["note"] = list(coordinate.df["note"].copy(deep=True))

    def get_atoms(self, labels):
        return list(self.df.loc[labels, "atom"])

    def get_notes(self, labels):
        return list(self.df.loc[labels, "note"].copy(deep=True))

    def copy(self):
        return Coordinate(self.df.copy(deep=True))

    def extract(self, indices):
        df = self.df.loc[indices, :]
        return Coordinate(df.copy(deep=True))
    
    def drop(self, indices):
        self.df.drop(index=indices, inplace=True)

    def drop_notes(self):
        self.df.drop(columns="note", inplace=True)

    def reset_label(self):
        self.df.set_axis(range(1, len(self.df)+1), axis=0, inplace=True)

    def get_indices(self, atom):
        return list(self.df.index[self.df["atom"]==atom])

    def get_to_x_axis(self, indices, x, angle):
        xyzs = self.df.loc[indices, ["x", "y", "z"]].to_numpy()
        thetax = angle * 2.0 * np.pi / 360.0
        matrix = at.get_to_x_axis(xyzs, x, thetax)
        return matrix

    def get_to_x_axis_minus(self, indices, x, angle):
        xyzs = self.df.loc[indices, ["x", "y", "z"]].to_numpy()
        thetax = angle * 2.0 * np.pi / 360.0
        matrix = at.get_to_x_axis_minus(xyzs, x, thetax)
        return matrix

    def transform(self, matrix):
        xyzs = self.df[["x", "y", "z"]].to_numpy()
        xyzs = at.transform(xyzs, matrix)
        self.df[["x", "y", "z"]] = xyzs

    def transform_to_x_axis(self, ref_indices, x, angle):
        matrix = self.get_to_x_axis(ref_indices, x, angle)
        self.transform(matrix)

    def transform_to_x_axis_minus(self, ref_indices, x, angle):
        matrix = self.get_to_x_axis_minus(ref_indices, x, angle)
        self.transform(matrix)

    def transform_partial(self, indices, matrix):
        xyzs = self.df.loc[indices, ["x", "y", "z"]].to_numpy()
        xyzs = at.transform(xyzs, matrix)
        self.df.loc[indices, ["x", "y", "z"]] = xyzs

    def transform_to_x_axis_partial(
            self, affine_indices, ref_indices, x, angle):
        matrix = self.get_to_x_axis(ref_indices, x, angle)
        self.transform_partial(affine_indices, matrix)

    def transform_to_x_axis_minus_partial(
            self, affine_indices, ref_indices, x, angle):
        matrix = self.get_to_x_axis_minus(ref_indices, x, angle)
        self.transform_partial(affine_indices, matrix)

    def get_distance_matrix(self):
        xyzs = self.df[["x", "y", "z"]].to_numpy()
        return distance.cdist(xyzs, xyzs)

    def get_covalent_radii_matrix(self):
        covalent_radii = [
            data.covalent_radii(atom) for atom in self.df["atom"]
        ]
        xx, yy = np.meshgrid(covalent_radii, covalent_radii)
        matrix = xx + yy
        np.fill_diagonal(matrix, 0.0)
        return matrix

    def get_adj_matrix(self, threshold=1.25):
        r_matrix = self.get_distance_matrix()
        R_matrix = self.get_covalent_radii_matrix()
        adj_matrix = (r_matrix < R_matrix * threshold)
        return adj_matrix.astype(int)

    def is_overlapped(self, threshold=0.8):
        r_matrix = self.get_distance_matrix()
        R_matrix = self.get_covalent_radii_matrix()
        return np.any(r_matrix < R_matrix * threshold)

    def to_graph(self):
        matrix = self.get_adj_matrix()
        G = nx.from_numpy_matrix(matrix, create_using=nx.Graph())
        for _, _, weight in G.edges(data=True):
            weight.clear()

        # Relabel nodes with self.df.index.
        mapping = dict(zip(G, self.df.index))
        G = nx.relabel_nodes(G, mapping)

        # Set "atom" attributes to nodes.
        node_attrs = self.df["atom"].to_dict()
        nx.set_node_attributes(G, node_attrs, "atom")

        return G

    def separate(self):
        G = self.to_graph()
        indices_list = [
            list(indices) for indices in nx.connected_components(G)
        ]
        indices_list = sorted(indices_list, key=len, reverse=True)
        coordinates = [self.extract(indices) for indices in indices_list]
        return coordinates

    def get_distance(self, index0, index1):
        xyzs = self.df.loc[[index0, index1], ["x", "y", "z"]].to_numpy()
        return distance.euclidean(xyzs[0], xyzs[1])

    def get_covalent_bond_length(self, index0, index1):
        covalent_radius0 = data.covalent_radii(self.df.at[index0, "atom"])
        covalent_radius1 = data.covalent_radii(self.df.at[index1, "atom"])
        return covalent_radius0 + covalent_radius1

    def has_a_bond(self, index0, index1, threshold=1.25):
        r = self.get_distance(index0, index1)
        R = self.get_covalent_bond_length(index0, index1)
        return r < R * threshold

    def adj(self, index, threshold=1.25):
        xyz = self.df.loc[[index], ["x", "y", "z"]].to_numpy()
        xyzs = self.df[["x", "y", "z"]].to_numpy()
        rs = distance.cdist(xyz, xyzs)[0]

        covalent_radius = data.covalent_radii(self.df.at[index, "atom"])
        covalent_radii = np.array(
            [data.covalent_radii(atom) for atom in self.df["atom"]]
        )
        Rs = covalent_radius + covalent_radii

        indices = list(self.df.index[rs < Rs * threshold])
        indices.remove(index)

        return indices

    def get_angle(self, index0, index1, index2):
        xyzs = self.df.loc[
            [index0, index1, index2], ["x", "y", "z"]
        ].to_numpy()
        v0 = xyzs[0] - xyzs[1]
        v1 = xyzs[2] - xyzs[1]
        inner_product = np.inner(v0, v1)
        l0 = np.linalg.norm(v0)
        l1 = np.linalg.norm(v1)
        theta = np.arccos(inner_product/(l0*l1))
        return np.rad2deg(theta)

    def get_dihedral_angle(self, index0, index1, index2, index3):
        xyzs = self.df.loc[
            [index0, index1, index2, index3], ["x", "y", "z"]
        ].to_numpy()
        theta = at.get_dihedral_angle(xyzs[[1, 2, 0, 3]])
        return np.rad2deg(theta)

    def save(self, path, header=None, footer=None, ignore_notes=False):

        if header == None:
            header = []
        if footer == None:
            footer = []

        if ignore_notes:
            lines = [
                "{:2s}   {: 16.12f}   {: 16.12f}   {: 16.12f}\n".format(
                    *self.df.iloc[i, 0:4]
                )
                for i in range(len(self.df))
            ]
        else:
            lines = [
                "{:2s}   {: 16.12f}   {: 16.12f}   {: 16.12f}{}\n".format(
                    *self.df.iloc[i, 0:4],
                    "".join([f" {note}" for note in self.df.iat[i, 4]])
                )
                for i in range(len(self.df))
            ]

        with open(path, "w") as f:
            f.writelines(header)
            f.writelines(lines)
            f.writelines(footer)

    def show(self, name=None):
        trash_path = f"/Users/takahirodoba/trash/{name}.com"
        header = ["# STO-3G\n", "\n", "title\n", "\n", "0 1\n"]
        footer = ["\n"]
        self.save(trash_path, header, footer, ignore_notes=True)
        subprocess.run(
            f"open -a /Applications/gv/gview.app {trash_path}", shell=True
        )


def concat(coordinates):
    """
    coordinates : list of Coordinate object
    """
    dfs = [coordinate.df for coordinate in coordinates]
    df = pd.concat(dfs, axis=0, copy=True)
    return Coordinate(df)


def same_structure(coordinate0, coordinate1):
    """
    Loose comparision of two coordinates.
    """
    graph0 = coordinate0.to_graph()
    graph1 = coordinate1.to_graph()
    nm = iso.categorical_node_match("atom", None)
    return nx.is_isomorphic(graph0, graph1, node_match=nm)


def same_connectivity(coordinate0, coordinate1):
    """
    Strict comparision of two coordinates.
    """
    label_bool = all(coordinate0.df.index == coordinate1.df.index)
    matrix0 = coordinate0.get_adj_matrix()
    matrix1 = coordinate1.get_adj_matrix()
    matrix_bool = all(matrix0 == matrix1)
    return label_bool and matrix_bool
