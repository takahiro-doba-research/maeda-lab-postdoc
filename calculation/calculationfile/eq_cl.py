import os
import re

import pandas as pd


def read_eq_cl(path):
    file = EQCl()
    file.name = os.path.splitext(os.path.basename(path))[0]

    with open(path, "r") as f:
        lines = f.readlines()

    file.neq = int(lines[0].strip())

    data = [
        int(line.strip()) for line in lines[1:]
        if re.fullmatch(r"\d+\n", line)
    ]
    data = [[i, line] for i, line in enumerate(data)]
    file.df = pd.DataFrame(
        data, columns=["EQ", "group"]
    ).set_index("EQ")
    
    return file


class EQCl:

    def __init__(self, name=None, neq=0, df=None):
        self.name = name
        self.neq = neq
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

