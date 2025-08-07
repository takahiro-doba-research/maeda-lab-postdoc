import os
import re

import numpy as np
import pandas as pd

from .coordinate import Coordinate
from .output import Output
from .optimization_path import OptPath


def read_grrm_lup_output(
        path, lups0=True, lups1=True, ircs=True,
        status=True, time=True):
    
    file = GrrmLupOutput()
    file.name = os.path.splitext(os.path.basename(path))[0]

    with open(path, "r") as f:
        lines = f.readlines()

    if lups0:
        start_indices = [
            i for i, line in enumerate(lines)
            if "LUP-path optimization" in line
        ]
        stop_indices = [
            i for i, line in enumerate(lines)
            if "The first criteria met" in line
        ]
        lup0_block = lines[start_indices[0]:stop_indices[0]]
        itr_indices = [
            i for i, line in enumerate(lup0_block)
            if "LUP-path optimization" in line
        ]
        profile_indices = [
            i for i, line in enumerate(lup0_block)
            if "---Profile of LUP path" in line
        ]
        blank_indices = [
            i for i, line in enumerate(lup0_block)
            if line == "\n"
        ]
        stop_blank_indices = [
            [
                blank_index for blank_index in blank_indices
                if blank_index > profile_index
            ][0]
            for profile_index in profile_indices
        ]

        for i in range(len(itr_indices)):
            opt_path = OptPath()
            itr = int(
                re.findall(
                    r"ITR.(.+)of", lup0_block[itr_indices[i]]
                )[0].strip()
            )

            # Set OptPath.name
            opt_path.name = itr

            # Set OptPath.outputs
            path_block = lup0_block[itr_indices[i]:stop_blank_indices[i]]
            node_indices = [
                i for i, line in enumerate(path_block)
                if "#" in line
            ]
            energy_indices = [
                i for i, line in enumerate(path_block)
                if "ENERGY" in line
            ]
            opt_profile_indices = [
                i for i, line in enumerate(path_block)
                if "---Profile of LUP path" in line
            ]

            for j in range(len(node_indices)):
                output = Output()
                node = int(
                    re.findall(
                        r"# NODE(.+)", path_block[node_indices[j]]
                    )[0].strip()
                )

                # Set Output.name
                output.name = node

                # Set Output.coordinate
                coordinate = path_block[node_indices[j]+1:energy_indices[j]]
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
                output.coordinate = Coordinate(df)

                # Set Output.E
                output.E = float(
                    re.findall(
                        r"=(.+)", path_block[energy_indices[j]]
                    )[0].strip()
                )

                # Append Output object
                opt_path.outputs[node] = output

            # Set OptPath.profile
            profile_lines = path_block[opt_profile_indices[0]+2:]
            profile_data = [line.split() for line in profile_lines]
            profile_data = [
                [int(line[0]), float(line[1]), float(line[2])]
                for line in profile_data
            ]
            profile_df = pd.DataFrame(
                profile_data,
                columns=["Node", "Length (ang)", "Energy (hartree)"]
            )
            opt_path.profile = profile_df

            # Append OptPath object
            file.lups0[itr] = opt_path

    if lups1:
        start_indices = [
            i for i, line in enumerate(lines)
            if "The first criteria met" in line
        ]
        stop_indices = [
            i for i, line in enumerate(lines)
            if "---Approximate TS geometry" in line
        ]
        lup1_block = lines[start_indices[0]+3:stop_indices[0]]
        itr_indices = [
            i for i, line in enumerate(lup1_block)
            if "LUP-path optimization" in line
        ]
        profile_indices = [
            i for i, line in enumerate(lup1_block)
            if "---Profile of LUP path" in line
        ]
        blank_indices = [
            i for i, line in enumerate(lup1_block)
            if line == "\n"
        ]
        stop_blank_indices = [
            [
                blank_index for blank_index in blank_indices
                if blank_index > profile_index
            ][0]
            for profile_index in profile_indices
        ]

        for i in range(len(itr_indices)):
            opt_path = OptPath()
            itr = int(
                re.findall(
                    r"ITR.(.+)of", lup1_block[itr_indices[i]]
                )[0].strip()
            )

            # Set OptPath.name
            opt_path.name = itr

            # Set OptPath.outputs
            path_block = lup1_block[itr_indices[i]:stop_blank_indices[i]]
            node_indices = [
                i for i, line in enumerate(path_block)
                if "#" in line
            ]
            energy_indices = [
                i for i, line in enumerate(path_block)
                if "ENERGY" in line
            ]
            opt_profile_indices = [
                i for i, line in enumerate(path_block)
                if "---Profile of LUP path" in line
            ]

            for j in range(len(node_indices)):
                output = Output()
                node = int(
                    re.findall(
                        r"# NODE(.+)", path_block[node_indices[j]]
                    )[0].strip()
                )

                # Set Output.name
                output.name = node

                # Set Output.coordinate
                coordinate = path_block[node_indices[j]+1:energy_indices[j]]
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
                output.coordinate = Coordinate(df)

                # Set Output.E
                output.E = float(
                    re.findall(
                        r"=(.+)", path_block[energy_indices[j]]
                    )[0].strip()
                )

                # Append Output object
                opt_path.outputs[node] = output

            # Set OptPath.profile
            profile_lines = path_block[opt_profile_indices[0]+2:]
            profile_data = [line.split() for line in profile_lines]
            profile_data = [
                [int(line[0]), float(line[1]), float(line[2])]
                for line in profile_data
            ]
            profile_df = pd.DataFrame(
                profile_data,
                columns=["Node", "Length (ang)", "Energy (hartree)"]
            )
            opt_path.profile = profile_df

            # Append OptPath object
            file.lups1[itr] = opt_path

        # Set OptPath.eqs for the OptPath
        app_eq_indices = [
            i for i, line in enumerate(lines)
            if "---Approximate EQ geometry" in line
        ]
        app_ts_indices = [
            i for i, line in enumerate(lines)
            if "---Approximate TS geometry" in line
        ]
        energy_indices = [
            i for i, line in enumerate(lines)
            if "ENERGY" in line
        ]
        app_eq_stop_index = [
            [
                energy_index for energy_index in energy_indices
                if energy_index > app_eq_index
            ][0]
            for app_eq_index in app_eq_indices
        ]
        app_ts_stop_index = [
            [
                energy_index for energy_index in energy_indices
                if energy_index > app_ts_index
            ][0]
            for app_ts_index in app_ts_indices
        ]

        for i in range(len(app_eq_indices)):
            eq = Output()

            # Set Output.name
            eq.name = int(
                re.findall(r"NODE(.+)\)", lines[app_eq_indices[i]])[0].strip()
            )

            # Set Output.coordinate
            coordinate = lines[app_eq_indices[i]+1:app_eq_stop_index[i]]
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
                re.findall(r"\((.+):", lines[app_eq_stop_index[i]])[0].strip()
            )

            # Append to OptPath.eqs
            last_itr = max(file.lups1.keys())
            file.lups1[last_itr].eqs[eq.name] = eq

        for i in range(len(app_ts_indices)):
            ts = Output()

            # Set Output.name
            ts.name = int(
                re.findall(r"NODE(.+)\)", lines[app_ts_indices[i]])[0].strip()
            )

            # Set Output.coordinate
            coordinate = lines[app_ts_indices[i]+1:app_ts_stop_index[i]]
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
            ts.coordinate = Coordinate(df)

            # Set Output.E
            ts.E = float(
                re.findall(r"\((.+):", lines[app_ts_stop_index[i]])[0].strip()
            )

            # Append to OptPath.tss
            last_itr = max(file.lups1.keys())
            file.lups1[last_itr].tss[ts.name] = ts
            
    if status:
        normal_index = [line for line in lines if "Normal termination" in line]

        if normal_index:
            file.status = "Normal termination"

        else:
            file.status = "Error termination"

    if time:
        time = [line for line in lines if "TOTAL ELAPSED TIME" in line]

        if time:
            file.time = float(re.findall(r":(.+)SEC.", time[0])[0].strip())

    return file


class GrrmLupOutput:

    def __init__(
            self, name=None, lups0=None, lups1=None, ircs=None,
            status=None, time=None):
            
        self.name = name
        self.lups0 = lups0 or {}  # Dictionary of OptPath objects
        self.lups1 = lups1 or {}  # Dictionary of OptPath objects
        self.ircs = ircs or {}  # Dictionary of OptPath objects
        self.status = status
        self.time = time
