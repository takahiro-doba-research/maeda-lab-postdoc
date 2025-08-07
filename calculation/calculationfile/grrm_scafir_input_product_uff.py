import os
import re

import pandas as pd

from .coordinate_uff import CoordinateUff


def read_grrm_scafir_input_product_uff(path):
    file = GrrmScafirInputProductUff()
    file.name = os.path.splitext(os.path.basename(path))[0]
    
    with open(path, "r") as f:
        lines = f.readlines()

    sharp_index = [i for i, line in enumerate(lines) if "#" in line][0]
    options_index = [
        i for i, line in enumerate(lines) if "product" in line.lower()
    ][0]

    start_index = sharp_index + 3
    stop_index = options_index

    # Set GaussianSpInput.coordinate
    coordinate = lines[start_index:stop_index]
    coordinate = [line.strip().split() for line in coordinate]
    coordinate = [
        [i+1,
         line[0], float(line[1]), float(line[2]), float(line[3]),
         line[4:]]
        for i, line in enumerate(coordinate)
    ]
    df = pd.DataFrame(
        coordinate,
        columns=["label", "atom", "x", "y", "z", "note"]
    ).set_index("label")
    file.coordinate = CoordinateUff(df)

    file.header = lines[0:start_index]
    file.footer = lines[stop_index:]
    
    return file


class GrrmScafirInputProductUff:

    def __init__(self, name=None, header=None, coordinate=None, footer=None):
        self.name = name
        self.header = header or []
        self.coordinate = coordinate
        self.footer = footer or []

        if not isinstance(self.coordinate, CoordinateUff):
            self.coordinate = CoordinateUff()
