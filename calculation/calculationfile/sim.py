import os
import re

import pandas as pd
from IPython.display import display


def read_sim(path):
    file = Sim()
    file.name = os.path.splitext(os.path.basename(path))[0]

    with open(path, "r") as f:
        lines = f.readlines()

    """
    Python 3.8 -
    data = [
        [int(value[0][0]), float(value[0][1])]
        for line in lines
        if (value := re.findall(
            r"Traffic volume for EQ -\s*(\d+) : (.+)", line
        ))
    ]
    """
    data = []

    for line in lines:
        value = re.findall(r"Traffic volume for EQ -\s*(\d+) : (.+)", line)

        if value:
            data.append([int(value[0][0]), float(value[0][1])])

    file.df = pd.DataFrame(
        data, columns=["EQ", "traffic"]
    ).set_index("EQ")
    
    return file


class Sim:

    def __init__(self, name=None, df=None):
        self.name = name
        self.df = df

    def __getitem__(self, item):

        if isinstance(item, int):
            return self.df.at[item, "traffic"]

        else:
            return list(self.df.loc[item, "traffic"])

    def max(self):
        return self.df["traffic"].max()

    def max_eq(self):
        return self.df["traffic"].idxmax()

    def max_eqs(self):
        max_traffic = self.df["traffic"].max()
        max_eqs = self.df.index[self.df["traffic"]==max_traffic]
        return list(max_eqs)

    def to_dict(self):
        return self.df["traffic"].to_dict()


def max_sim(files):
    merged_df = files[0].df.copy()

    for file in files[1:]:
        merged_df = pd.merge(merged_df, file.df, on="EQ", how="outer")

    s = merged_df.max(axis=1)
    df = pd.DataFrame(s, columns=["traffic"])

    return Sim(df=df)
