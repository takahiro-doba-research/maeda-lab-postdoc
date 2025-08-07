import re

from .coordinate import Coordinate


class CoordinateUff(Coordinate):

    def __init__(self, df=None):
        super().__init__(df)

    def simplify_atom(self):
        self.df["atom"] = [
            re.findall(r"(.+?)-", atom)[0]
            for atom in self.df["atom"]
        ]
