import os
import re
import time

import numpy as np
import pandas as pd

from .coordinate import Coordinate
from .output import Output
from .util import join_group_data


def read_pt_list(path):

    pt_list = PTList()
    pt_list.name = os.path.splitext(os.path.basename(path))[0]

    with open(path, "r") as f:
        lines = f.readlines()

    pt_indices = [i for i, line in enumerate(lines) if "#" in line]

    if pt_indices:

        for i in range(len(pt_indices)):
            pt = Output()

            if i != len(pt_indices) - 1:
                pt_block = lines[pt_indices[i]:pt_indices[i+1]-1]
            else:
                pt_block = lines[pt_indices[i]:]

            # Set Output.name
            pt.name = int(re.findall(r"TS(.+),", pt_block[0])[0].strip())
            
            # Set Output.coordinate
            E_index = [
                i for i, line in enumerate(pt_block) if "Energy" in line
            ]
            coordinate = pt_block[1:E_index[0]]
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
            pt.coordinate = Coordinate(df)

            # Set Output.E
            pt.E = float(
                re.findall(r"\((.+):", pt_block[E_index[0]])[0].strip()
            )

            # Set Output.spin
            spin = [line for line in pt_block if "Spin(**2)" in line]
            pt.spin = float(re.findall(r"=(.+)", spin[0])[0].strip())

            # Set Output.E_zpv
            E_zpv = [line for line in pt_block if "ZPVE" in line]
            pt.E_zpv = float(re.findall(r"=(.+)", E_zpv[0])[0].strip())

            # Set Output.normal_mode
            normal_mode_index = [
                i for i, line in enumerate(pt_block)
                if "Normal mode eigenvalues" in line
            ]
            normal_mode = pt_block[normal_mode_index[0]+1:-1]
            normal_mode = [
                float(value)
                for line in normal_mode
                for value in line.strip().split()
            ]
            pt.normal_mode = np.array(normal_mode)

            # Set Output.connection
            connection = re.findall(r":(.+)-(.+)", pt_block[-1])[0]
            connection = [
                node.strip() if "??" in node else int(node.strip())
                for node in connection
            ]
            pt.connection = tuple(connection)

            # Add to pt_list
            pt_list[pt.name] = pt

    return pt_list


class PTList:

    def __init__(self, name=None, pts=None):
        self.name = name
        self.pts = pts or {}

    def __len__(self):
        return len(self.pts)

    def __getitem__(self, key):
        return self.pts[key]

    def __setitem__(self, key, value):
        self.pts[key] = value

    def keys(self):
        return self.pts.keys()

    def values(self):
        return self.pts.values()

    def items(self):
        return self.pts.items()

    def show(self, names):

        if isinstance(names, int):
            names = [names]

        for name in names:
            pt = self.pts[name]
            pt.coordinate.show(f"PT{name}")
            time.sleep(2)

    def shorter_distance(self, index0, index1, distance):
        pts_shorter_distance = [
            name for name, pt in self.items()
            if pt.coordinate.get_distance(index0, index1) < distance
        ]
        return pts_shorter_distance

    def longer_distance(self, index0, index1, distance):
        pts_longer_distance = [
            name for name, pt in self.items()
            if distance < pt.coordinate.get_distance(index0, index1)
        ]
        return pts_longer_distance

    def between_distance(self, index0, index1, distance0, distance1):
        pts_between_distance = [
            name for name, pt in self.items()
            if distance0 < pt.coordinate.get_distance(index0, index1) < distance1
        ]
        return pts_between_distance

    def between_dihedral_angle(
            self, index0, index1, index2, index3,
            dihedral_angle0, dihedral_angle1):
        pts_between_dihedral_angle = [
            name for name, pt in self.items()
            if dihedral_angle0 < pt.coordinate.get_dihedral_angle(
                index0, index1, index2, index3) < dihedral_angle1
        ]
        return pts_between_dihedral_angle

    def adjust_energy(self, name=None):

        if name == None:
            eq_std = min(self.values(), key=lambda pt:pt.E)
            E_std = eq_std.E

        else:
            E_std = self[name].E

        for pt in self.values():
            pt.E = (pt.E - E_std) * 2625.5

    def save_pts(self, path):
        os.makedirs(f"{path}/PT", exist_ok=True)
        header = ["# STO-3G\n", "\n", "title\n", "\n", "0 1\n"]
        footer = ["\n"]

        for name, pt in self.items():
            pt.coordinate.save(
                path=f"{path}/PT/{name}.com",
                header=header,
                footer=footer,
                ignore_notes=True
            )
