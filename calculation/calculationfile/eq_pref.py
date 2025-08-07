import os
import re

import pandas as pd


def read_eq_pref(path):
    file = EQPref()
    file.name = os.path.splitext(os.path.basename(path))[0]

    with open(path, "r") as f:
        lines = f.readlines()

    file.neq = int(re.findall(r"# of EQs: (.+)\n", lines[0])[0])

    data = [
        [i, int(line.split()[0]), float(line.split()[1])]
        for i, line in enumerate(lines[1:])
    ]
    file.df = pd.DataFrame(
        data, columns=["EQ", "tag", "preference"]
    ).set_index("EQ")
    
    return file


class EQPref:

    def __init__(self, name=None, neq=0, df=None):
        self.name = name
        self.neq = neq
        self.df = df

    def __getitem__(self, item):

        if isinstance(item, int):
            return self.df.at[item, "preference"]

        else:
            return list(self.df.loc[item, "preference"])

    def max(self):
        return self.df["preference"].max()

    def max_eq(self):
        return self.df["preference"].idxmax()

    def max_eqs(self):
        max_pref = self.df["preference"].max()
        max_eqs = self.df.index[self.df["preference"]==max_pref]
        return list(max_eqs)

    def tag(self, eqs):

        if isinstance(eqs, int):
            return self.df.at[eqs, "tag"]

        else:
            return list(self.df.loc[eqs, "tag"])

    def to_dict(self):
        return self.df["preference"].to_dict()

    def save(self, path):
        lines = [
            f"{tag:3d} {value:18.12f}\n"
            for tag, value in zip(self.df["tag"], self.df["preference"])
        ]

        with open(path, "w") as f:
            f.write(f"# of EQs: {self.neq}\n")
            f.writelines(lines)

    
