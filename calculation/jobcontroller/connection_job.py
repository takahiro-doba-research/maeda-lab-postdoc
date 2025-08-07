import glob
import os
from multiprocessing import Pool

import pandas as pd
from IPython.display import display

import calculationfile as cf

from .gaussian_sp_grrm_min_job import gaussian_sp_grrm_min_job
from .util import split_separated_group_info_df, dict_product


class ConnectionJob:

    def __init__(self, wd):
        self.wd = wd

    def add_note(self, note_path):
        # Read csv files.
        separated_group_df = pd.read_csv(
            f"{self.wd}/separated_group.csv",
            index_col="separated_EQ"
        )
        separated_group_info_df = pd.read_csv(
            f"{self.wd}/separated_group_info.csv"
        )
        separated_group_info_df = split_separated_group_info_df(
            separated_group_info_df
        ).set_index("separated_group")

        # Read the note file.
        note_file = cf.read_gaussian_sp_input_connection(
            f"{self.wd}/{note_path}"
        )

        paths = sorted(
            glob.glob(f"{self.wd}/separated_EQ/*"),
            key=lambda x:int(os.path.splitext(os.path.basename(x))[0])
        )
        
        for path in paths:
            file = cf.read_gaussian_sp_input(path)
            group = separated_group_df.at[int(file.name), "separated_group"]
            labels = separated_group_info_df.at[group, "labels"]
            file.coordinate.note = note_file.coordinate.get_notes(labels)
            file.coordinate.save(
                path,
                header=file.header,
                footer=file.footer,
                ignore_notes=False
            )

    def structure_check(self, note_path, sub_folders):
        """
        Check whether all of the notes are given correctly.
        """
        # Read the structure in the note file.
        note_file = cf.read_gaussian_sp_input_connection(
            f"{self.wd}/{note_path}"
        )
        note_file.coordinate.read_connectivity()

        # Read substituents.
        sub_files_2d = {}

        for sub_folder in sub_folders:
            sub_paths = sorted(
                glob.glob(f"{self.wd}/{sub_folder}/*"),
                key=lambda x:int(os.path.splitext(os.path.basename(x))[0])
            )
            sub_files = []

            for sub_path in sub_paths:
                sub_file = cf.read_gaussian_sp_input_connection(sub_path)
                sub_file.coordinate.read_connectivity()
                sub_files.append(sub_file)

            sub_files_2d[sub_folder] = sub_files

        for substituent in sub_files_2d.keys():
            my_folder = f"{self.wd}/structure_check/{substituent}"
            os.makedirs(my_folder, exist_ok=True)
            other_substituents = [
                sub for sub in sub_files_2d.keys() if sub != substituent
            ]
            other_files = [sub_files_2d[sub][0] for sub in other_substituents]
            other_coordinates = [file.coordinate for file in other_files]

            for file in sub_files_2d[substituent]:
                cooridnates = [
                    note_file.coordinate, file.coordinate, *other_coordinates
                ]
                coordinate = cf.connect(cooridnates)
                coordinate.save(
                    f"{my_folder}/{file.name}.com",
                    header=note_file.header,
                    footer=note_file.footer,
                    ignore_notes=True
                )

    def make_conformer(self, sub_folders):
        """
        Parameters
        --------
        sub_folders : list
            List of the names of the folders of substituents.
            e.g. ["arene", "alkene", "PG", "backbone", "pyridone"]
        """
        # Read csv files.
        # separated_EQ -> separated_group
        # separated_group -> separated_group_info
        separated_group_df = pd.read_csv(
            f"{self.wd}/separated_group.csv",
            index_col="separated_EQ"
        )
        separated_group_info_df = pd.read_csv(
            f"{self.wd}/separated_group_info.csv"
        )
        separated_group_info_df = split_separated_group_info_df(
            separated_group_info_df
        ).set_index("separated_group")

        # Determine my_name of each substituent and return the mapping.
        # e.g. 
        # mapping = {"arene":1, "alkene":2, "PG":3, "backbone":4, "pyridone":5}
        mapping = {}

        for sub_folder in sub_folders:
            first_path = glob.glob(f"{self.wd}/{sub_folder}/*")[0]
            first_file = cf.read_gaussian_sp_input_connection(first_path)
            first_file.coordinate.read_connectivity()
            mapping[sub_folder] = first_file.coordinate.my_name

        # Read substituents.
        sub_files_2d = {}

        for sub_folder in sub_folders:
            sub_paths = sorted(
                glob.glob(f"{self.wd}/{sub_folder}/*"),
                key=lambda x:int(os.path.splitext(os.path.basename(x))[0])
            )
            sub_files = []

            for sub_path in sub_paths:
                sub_file = cf.read_gaussian_sp_input_connection(sub_path)
                sub_file.coordinate.read_connectivity()
                sub_files.append(sub_file)

            sub_files_2d[sub_folder] = sub_files

        derivative_data = []
        conformer_data = []
        derivative_count = 0
        conformer_count = 0

        separated_eq_paths = sorted(
            glob.glob(f"{self.wd}/separated_EQ/*"),
            key=lambda x:int(os.path.splitext(os.path.basename(x))[0])
        )

        for separated_eq_path in separated_eq_paths:
            separated_eq_file = cf.read_gaussian_sp_input_connection(
                separated_eq_path
            )
            separated_group = separated_group_df.at[
                int(separated_eq_file.name), "separated_group"
            ]
            labels = separated_group_info_df.at[
                separated_group, "labels"
            ]
            separated_eq_file.coordinate.label = labels
            separated_eq_file.coordinate.read_connectivity()

            used_sub_folders = [
                folder for folder, name in mapping.items()
                if name in separated_eq_file.coordinate.other_names
            ]
            used_sub_files_2d = {
                folder:files for folder, files in sub_files_2d.items()
                if folder in used_sub_folders
            }
            """
            e.g.
            used_sub_files_2d = {
                "arene":[<GaussianSpInputConnection object at xxxxxxxxxxx>],
                "alkene":[<GaussianSpInputConnection object at xxxxxxxxxxx>],
                "PG":[<GaussianSpInputConnection object at xxxxxxxxxxx>],
                "backbone":[
                    <GaussianSpInputConnection object at xxxxxxxxxxx>,
                    <GaussianSpInputConnection object at xxxxxxxxxxx>,
                    <GaussianSpInputConnection object at xxxxxxxxxxx>,
                    ...
                ],
                "pyridone":[
                    <GaussianSpInputConnection object at xxxxxxxxxxx>,
                    <GaussianSpInputConnection object at xxxxxxxxxxx>,
                    <GaussianSpInputConnection object at xxxxxxxxxxx>,
                    ...
                ],
            }
            """

            for used_sub_files in dict_product(used_sub_files_2d):
                os.makedirs(f"{self.wd}/derivative/{derivative_count}")
                derivative_data.append({
                    "separated_EQ":separated_eq_file.name,
                    **{
                        folder:file.name
                        for folder, file in used_sub_files.items()
                    },
                    "derivative":derivative_count
                })

                dihedral_angles_2d = {
                    folder:dihedral_angles[2]["dihedral_angles"]
                    for folder, file in used_sub_files.items()
                    for dihedral_angles in file.coordinate.dihedral_angles
                }
                """
                e.g.
                dihedral_angles_2d = {
                    "arene":[0],
                    "alkene":[0],
                    "PG":[0],
                    "backbone":[0, 120, 240],
                    "pyridone":[0, 120, 240],
                }
                """

                for dihedral_angles in dict_product(dihedral_angles_2d):
                    files = [separated_eq_file, *used_sub_files.values()]
                    coordinates = [file.coordinate for file in files]
                    dihedral_angles_edge = [
                        (
                            mapping[folder],
                            separated_eq_file.coordinate.my_name,
                            {"dihedral_angle":dihedral_angle}
                        )
                        for folder, dihedral_angle in dihedral_angles.items()
                    ]
                    connected_coordinate = cf.connect(
                        coordinates,
                        dihedral_angles=dihedral_angles_edge
                    )
                    output_path = "{}/derivative/{}/{}.com".format(
                        self.wd, derivative_count, conformer_count
                    )
                    connected_coordinate.save(
                        output_path,
                        header=separated_eq_file.header,
                        footer=separated_eq_file.footer,
                        ignore_notes=True
                    )

                    conformer_data.append({
                        "derivative":derivative_count,
                        **dihedral_angles,
                        "conformer":conformer_count
                    })

                    conformer_count += 1

                derivative_count += 1
        
        derivative_df = pd.DataFrame(derivative_data)
        derivative_df.to_csv(
            f"{self.wd}/derivative.csv",
            float_format="%.0f",
            columns=["separated_EQ", *sub_folders, "derivative"],
            index=False
        )
        conformer_df = pd.DataFrame(conformer_data)
        conformer_df.to_csv(
            f"{self.wd}/conformer.csv",
            float_format="%.0f",
            columns=["derivative", *sub_folders, "conformer"],
            index=False
        )

    def optimize(self, processes,
            sp_method_path, sp_option, min_method_path, min_option):
        sp_method = cf.read_gaussian_sp_input(f"{self.wd}/{sp_method_path}")
        min_method = cf.read_grrm_min_input(f"{self.wd}/{min_method_path}")
        wds = sorted(
            glob.glob(f"{self.wd}/derivative/*"),
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
            conformer_energy.csv
            derivetive_energy.csv
            conformer
        conformer.csv is required.
        """
        # Read conformer.csv.
        conformer_df = pd.read_csv(
            f"{self.wd}/conformer.csv",
            index_col="derivative"
        )

        # Data for conformer_energy.csv.
        conformer_energy_data = []

        # Data for derivative_energy.csv.
        derivative_energy_data = []

        # Folder for conformer.
        conformer_folder = f"{self.wd}/conformer"
        os.makedirs(conformer_folder, exist_ok=True)

        for derivative in set(conformer_df.index):
            derivative_line = {
                "derivative":derivative
            }
            conformers = list(
                conformer_df.loc[[derivative], "conformer"]
            )

            for conformer in conformers:
                conformer_line = {
                    "derivative":derivative,
                    "conformer":conformer
                }

                # Aanalyze SP output.
                sp_path = "{}/derivative/{}/{}_SP.log".format(
                    self.wd, derivative, conformer
                )

                if os.path.exists(sp_path):
                    sp_file = cf.read_gaussian_sp_output(
                        sp_path, pid=False, status=True, E=True, time=True
                    )

                    if sp_file.status == "Normal termination":
                        conformer_line["SP_status"] = sp_file.status
                        conformer_line["SP_E"] = sp_file.E
                        conformer_line["SP_time"] = sp_file.time

                    elif sp_file.status == "Error termination":
                        conformer_line["SP_status"] = sp_file.status
                        conformer_line["SP_time"] = sp_file.time

                    else:
                        pass

                # Analyze MIN output.
                min_path = "{}/derivative/{}/{}_MIN.log".format(
                    self.wd, derivative, conformer
                )

                if os.path.exists(min_path):
                    min_file = cf.read_grrm_min_output(
                        min_path,
                        status=True, itrs=False, opt=True,
                        hessian=False, freq=True, thermo=False, time=True
                    )

                    if (min_file.status == "Minimum point was found"
                            or min_file.status == "Saddle point was found"):
                        conformer_line["MIN_status"] = min_file.status
                        conformer_line["MIN_freq0"] = min_file.freq[0]
                        conformer_line["MIN_E"] = min_file.E
                        conformer_line["MIN_time"] = min_file.time

                        # Add data to derivative_energy.csv
                        derivative_line["conformer"] = str(conformer)
                        derivative_line["MIN_E"] = min_file.E

                        # Save the optimized coordinate to conformer folder
                        header = ["# STO-3G\n", "\n", "title\n", "\n", "0 1\n"]
                        footer = ["\n"]
                        min_file.coordinate.save(
                            f"{conformer_folder}/{conformer}.com",
                            header=header,
                            footer=footer,
                            ignore_notes=True
                        )

                    elif min_file.status == "The structure is dissociating":
                        conformer_line["MIN_status"] = min_file.status
                        conformer_line["MIN_time"] = min_file.time

                    elif min_file.status == "Error termination":
                        conformer_line["MIN_status"] = min_file.status

                    else:
                        pass

                conformer_energy_data.append(conformer_line)

            derivative_energy_data.append(derivative_line)

        conformer_energy_df = pd.DataFrame(conformer_energy_data)
        conformer_energy_df = conformer_energy_df.sort_values(
            ["derivative", "conformer"]
        )
        conformer_energy_df.to_csv(
            f"{self.wd}/conformer_energy.csv",
            columns=[
                "derivative", "conformer",
                "SP_status", "SP_E", "SP_time",
                "MIN_status", "MIN_freq0", "MIN_E", "MIN_time"
            ],
            index=False
        )
        derivative_energy_df = pd.DataFrame(derivative_energy_data)
        derivative_energy_df = derivative_energy_df.sort_values(
            "derivative"
        )
        derivative_energy_df.to_csv(
            f"{self.wd}/derivative_energy.csv",
            columns=["derivative", "conformer", "MIN_E"],
            index=False
        )
