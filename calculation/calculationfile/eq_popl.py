import os

import pandas as pd


def read_eq_popl(path):
    file = EQPopl()
    file.name = os.path.splitext(os.path.basename(path))[0]

    with open(path, "r") as f:
        lines = f.readlines()

    data = [
        [i, float(line.split()[-1])]
        for i, line in enumerate(lines)
    ]
    file.df = pd.DataFrame(
        data, columns=["EQ", "population"]
    ).set_index("EQ")
    
    return file


class EQPopl:

    def __init__(self, name=None, df=None):
        self.name = name
        self.df = df

    def __getitem__(self, item):

        if isinstance(item, int):
            return self.df.at[item, "population"]

        else:
            return list(self.df.loc[item, "population"])

    def max(self):
        return self.df["population"].max()

    def max_eq(self):
        return self.df["population"].idxmax()

    def max_eqs(self):
        max_popl = self.df["population"].max()
        max_eqs = self.df.index[self.df["population"]==max_popl]
        return list(max_eqs)

    def to_dict(self):
        return self.df["population"].to_dict()
