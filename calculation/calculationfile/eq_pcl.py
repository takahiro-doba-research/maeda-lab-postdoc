import os

import pandas as pd


def read_eq_pcl(path):
    file = EQPcl()
    file.name = os.path.splitext(os.path.basename(path))[0]

    with open(path, "r") as f:
        lines = f.readlines()

    data = [[i, int(line)] for i, line in enumerate(lines)]
    file.df = pd.DataFrame(
        data, columns=["EQ", "group"]
    ).set_index("EQ")
    
    return file


class EQPcl:

    def __init__(self, name=None, df=None):
        self.name = name
        self.df = df

    def __getitem__(self, item):

        if isinstance(item, int):
            return self.df.at[item, "group"]

        else:
            return list(self.df.loc[item, "group"])

    def to_dict(self):
        return self.df["group"].to_dict()

    def to_group_dict(self):
        group_dict = {}

        for eq, group in self.df.itertuples():
            
            if group in group_dict.keys():
                group_dict[group].append(eq)

            else:
                group_dict[group] = [eq]

        return group_dict

    def number_of_groups(self):
        return len(set(self.df["group"]))

    def eqs_in_group(self, group):
        eqs = [eq for eq, g in self.df.itertuples() if g == group]
        return eqs

