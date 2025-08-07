import os
import re

import numpy as np
import pandas as pd

from .coordinate import Coordinate
from .output import Output
from .optimization_path import OptPath


def read_pt(path):
    file = PT()
    file.name = os.path.splitext(os.path.basename(path))[0]

    with open(path, "r") as f:
        lines = f.readlines()

    node_indices = [i for i, line in enumerate(lines) if "#" in line]

    if node_indices:
        pt_path = OptPath()
        pt_path.name = file.name

        for i in range(len(node_indices)):
            node = Output()

            if i != len(node_indices) - 1:
                node_block = lines[node_indices[i]:node_indices[i+1]-1]
            else:
                node_block = lines[node_indices[i]:]

            # Set Output.name
            node.name = int(re.findall(r"NODE (\d+)", node_block[0])[0].strip())
            
            # Set Output.coordinate
            E_index = [
                i for i, line in enumerate(node_block) if "ENERGY" in line
            ]
            coordinate = node_block[1:E_index[0]-1]
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
            node.coordinate = Coordinate(df)

            # Set Output.E
            node.E = float(
                re.findall(r"\((.+):", node_block[E_index[0]])[0].strip()
            )

            # Set Output.spin
            spin = [line for line in node_block if "Spin(**2)" in line]
            node.spin = float(re.findall(r"\)(.+)", spin[0])[0].strip())

            # Set Output.normal_mode
            normal_mode_index = [
                i for i, line in enumerate(node_block)
                if "NORMAL MODE EIGENVALUE" in line
            ]
            normal_mode = node_block[normal_mode_index[0]+1:]
            normal_mode = [
                float(value)
                for line in normal_mode if line
                for value in line.strip().split()
            ]
            node.normal_mode = np.array(normal_mode)

            # Append the complete node
            pt_path.outputs[i] = node
        
        # Set EQList.nodes
        file.pt_path = pt_path

    return file


class PT:

    def __init__(self, name=None, pt_path=None):
        self.name = name
        self.pt_path = pt_path
