import glob
import os
import re
from multiprocessing import Pool

import numpy as np
import pandas as pd

import calculationfile as cf

from .gaussian_sp_grrm_min_job import gaussian_sp_grrm_min_job
from .util import join_separated_group_info


class SeparationJob:

    def __init__(self, wd):
        self.wd = wd

    def separate(self):
        os.makedirs(f"{self.wd}/separated_group")
        paths = sorted(
            glob.glob(f"{self.wd}/EQ/*.com"),
            key=lambda x:int(os.path.splitext(os.path.basename(x))[0])
        )
        separated_eq_data = []
        separated_group_data = []
        separated_group_info = []
        separated_eq_count = 0
        separated_group_count = 0

        for path in paths:
            file = cf.read_gaussian_sp_input(path)
            coordinates = file.coordinate.separate()

            for coordinate in coordinates:
                labels = coordinate.label
                adj_matrix = coordinate.get_adj_matrix()

                for line in separated_group_info:

                    if (labels == line["labels"]
                            and np.all(adj_matrix == line["adj_matrix"])):

                        output_folder = "{}/separated_group/{}".format(
                            self.wd, line["separated_group"]
                        )
                        coordinate.save(
                            f"{output_folder}/{separated_eq_count}.com",
                            header=file.header,
                            footer=file.footer,
                            ignore_notes=True
                        )

                        separated_eq_data.append({
                            "EQ":file.name,
                            "separated_EQ":separated_eq_count
                        })
                        separated_group_data.append({
                            "separated_EQ":separated_eq_count,
                            "separated_group":line["separated_group"]
                        })
                        
                        separated_eq_count += 1
                        break

                else:
                    output_folder = "{}/separated_group/{}".format(
                        self.wd, separated_group_count
                    )
                    os.makedirs(output_folder)
                    coordinate.save(
                        f"{output_folder}/{separated_eq_count}.com",
                        header=file.header,
                        footer=file.footer,
                        ignore_notes=True
                    )

                    separated_eq_data.append({
                        "EQ":file.name,
                        "separated_EQ":separated_eq_count
                    })
                    separated_group_data.append({
                        "separated_EQ":separated_eq_count,
                        "separated_group":separated_group_count
                    })
                    separated_group_info.append({
                        "separated_group":separated_group_count,
                        "labels":labels,
                        "adj_matrix":adj_matrix
                    })

                    separated_eq_count += 1
                    separated_group_count += 1

        separated_eq_df = pd.DataFrame(separated_eq_data)
        separated_eq_df.to_csv(
            f"{self.wd}/separated_EQ.csv",
            columns=["EQ", "separated_EQ"],
            index=False
        )
        separated_group_df = pd.DataFrame(separated_group_data)
        separated_group_df.to_csv(
            f"{self.wd}/separated_group.csv",
            columns=["separated_EQ", "separated_group"],
            index=False
        )
        separated_group_info = join_separated_group_info(separated_group_info)
        separated_group_info_df = pd.DataFrame(separated_group_info)
        separated_group_info_df.to_csv(
            f"{self.wd}/separated_group_info.csv",
            columns=["separated_group", "labels", "adj_matrix"],
            index=False
        )

    def optimize(self, processes,
            sp_method_path, sp_option, min_method_path, min_option):
        sp_method = cf.read_gaussian_sp_input(f"{self.wd}/{sp_method_path}")
        min_method = cf.read_grrm_min_input(f"{self.wd}/{min_method_path}")
        wds = sorted(
            glob.glob(f"{self.wd}/separated_group/*"),
            key=lambda x:int(os.path.basename(x))
        )
        iterable = [
            (wd, sp_method, sp_option, min_method, min_option) for wd in wds
        ]
        pool = Pool(processes=processes)
        pool.starmap(gaussian_sp_grrm_min_job, iterable, 1)

    def analyze(self):
        """
        Analyzes the optimization results and exports
            separated_EQ_energy.csv
            separated_group_energy.csv
            separated_EQ
        separated_group.csv is required.
        """
        # Read separated_group.csv.
        separated_group_df = pd.read_csv(
            f"{self.wd}/separated_group.csv",
            index_col="separated_group"
        )

        # Data for separated_EQ_energy.csv.
        separated_eq_energy_data = []

        # Data for separated_group_energy.csv.
        separated_group_energy_data = []

        # Folder for separated_EQ.
        separated_eq_folder = f"{self.wd}/separated_EQ"
        os.makedirs(separated_eq_folder, exist_ok=True)

        for separated_group in set(separated_group_df.index):
            separated_group_line = {
                "separated_group":separated_group
            }
            separated_eqs = list(
                separated_group_df.loc[[separated_group], "separated_EQ"]
            )

            for separated_eq in separated_eqs:
                separated_eq_line = {
                    "separated_group":separated_group,
                    "separated_EQ":separated_eq
                }

                # Aanalyze SP output.
                sp_path = "{}/separated_group/{}/{}_SP.log".format(
                    self.wd, separated_group, separated_eq
                )

                if os.path.exists(sp_path):
                    sp_file = cf.read_gaussian_sp_output(
                        sp_path, pid=False, status=True, E=True, time=True
                    )

                    if sp_file.status == "Normal termination":
                        separated_eq_line["SP_status"] = sp_file.status
                        separated_eq_line["SP_E"] = sp_file.E
                        separated_eq_line["SP_time"] = sp_file.time

                    elif sp_file.status == "Error termination":
                        separated_eq_line["SP_status"] = sp_file.status
                        separated_eq_line["SP_time"] = sp_file.time

                    else:
                        pass

                # Analyze MIN output.
                min_path = "{}/separated_group/{}/{}_MIN.log".format(
                    self.wd, separated_group, separated_eq
                )

                if os.path.exists(min_path):
                    min_file = cf.read_grrm_min_output(
                        min_path,
                        status=True, itrs=False, opt=True,
                        hessian=False, freq=True, thermo=False, time=True
                    )

                    if min_file.status == "Minimum point was found":
                        separated_eq_line["MIN_status"] = min_file.status
                        separated_eq_line["MIN_freq0"] = min_file.freq[0]
                        separated_eq_line["MIN_E"] = min_file.E
                        separated_eq_line["MIN_time"] = min_file.time

                        # Add data to separated_group_energy.csv
                        separated_group_line["separated_EQ"] = str(separated_eq)
                        separated_group_line["MIN_E"] = min_file.E

                        # Save the optimized coordinate to separated_EQ folder
                        header = ["# STO-3G\n", "\n", "title\n", "\n", "0 1\n"]
                        footer = ["\n"]
                        min_file.coordinate.save(
                            f"{separated_eq_folder}/{separated_eq}.com",
                            header=header,
                            footer=footer,
                            ignore_notes=True
                        )

                    elif min_file.status == "The structure is dissociating":
                        separated_eq_line["MIN_status"] = min_file.status
                        separated_eq_line["MIN_time"] = min_file.time

                    elif min_file.status == "Saddle point was found":
                        separated_eq_line["MIN_status"] = min_file.status
                        separated_eq_line["MIN_freq0"] = min_file.freq[0]
                        separated_eq_line["MIN_time"] = min_file.time

                    elif min_file.status == "Error termination":
                        separated_eq_line["MIN_status"] = min_file.status

                    else:
                        pass

                separated_eq_energy_data.append(separated_eq_line)

            separated_group_energy_data.append(separated_group_line)

        separated_eq_energy_df = pd.DataFrame(separated_eq_energy_data)
        separated_eq_energy_df = separated_eq_energy_df.sort_values(
            ["separated_group", "separated_EQ"]
        )
        separated_eq_energy_df.to_csv(
            f"{self.wd}/separated_EQ_energy.csv",
            columns=[
                "separated_group", "separated_EQ",
                "SP_status", "SP_E", "SP_time",
                "MIN_status", "MIN_freq0", "MIN_E", "MIN_time"
            ],
            index=False
        )
        separated_group_energy_df = pd.DataFrame(separated_group_energy_data)
        separated_group_energy_df = separated_group_energy_df.sort_values(
            "separated_group"
        )
        separated_group_energy_df.to_csv(
            f"{self.wd}/separated_group_energy.csv",
            columns=["separated_group", "separated_EQ", "MIN_E"],
            index=False
        )