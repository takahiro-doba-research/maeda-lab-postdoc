import os

import networkx as nx
import pandas as pd

from .coordinate_connection import CoordinateConnection


def read_gaussian_sp_input_connection(path):
    file = GaussianSpInputConnection()
    file.name = os.path.splitext(os.path.basename(path))[0]
    
    with open(path, "r") as f:
        lines = f.readlines()

    sharp_index = [i for i, line in enumerate(lines) if "#" in line][0]
    blank_indices = [i for i, line in enumerate(lines) if line == "\n"]

    start_index = sharp_index + 5
    stop_index = [index for index in blank_indices if index > start_index][0]

    # Set GaussianSpInput.coordinate
    coordinate = lines[start_index:stop_index]
    coordinate = [line.strip().split() for line in coordinate]
    coordinate = [
        [i+1, line[0], float(line[1]), float(line[2]), float(line[3]),
         [int(value) for value in line[4:]]]
        for i, line in enumerate(coordinate)
    ]
    df = pd.DataFrame(
        coordinate,
        columns=["label", "atom", "x", "y", "z", "note"]
    ).set_index("label")
    file.coordinate = CoordinateConnection(df)

    file.header = lines[0:start_index]
    file.footer = lines[stop_index:]
    
    return file


class GaussianSpInputConnection:

    def __init__(self, name=None, header=None, coordinate=None, footer=None):
        self.name = name
        self.header = header or []
        self.coordinate = coordinate
        self.footer = footer or []

        if not isinstance(self.coordinate, CoordinateConnection):
            self.coordinate = CoordinateConnection()
            
        

        

    