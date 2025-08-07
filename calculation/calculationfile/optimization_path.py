import subprocess

from .data import atomic_number


class OptPath:

    def __init__(
            self, name=None, outputs=None, eqs=None, tss=None, profile=None):

        self.name = name
        self.outputs = outputs or {}
        self.eqs = eqs or {}
        self.tss = tss or {}
        self.profile = profile

    def to_gv(self, path):
        header = [" #p\n\n"]
        footer = [
""" GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad
 Normal termination of Gaussian 16"""
        ]

        noutput = len(self.outputs)
        lines = []

        for key, output in self.outputs.items():
            coordinates = [
                " {: 6d}        {: 3d}           0     {: 11.6f} {: 11.6f} {: 11.6f}\n".format(
                    index, atomic_number(atom), x, y, z
                )
                for index, atom, x, y, z
                in zip(
                    output.coordinate.df.index,
                    output.coordinate.df["atom"],
                    output.coordinate.df["x"],
                    output.coordinate.df["y"],
                    output.coordinate.df["z"]
                )
            ]
            lines += [
f""" GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad
                          Input orientation:                          
 ---------------------------------------------------------------------
 Center     Atomic      Atomic             Coordinates (Angstroms)
 Number     Number       Type             X           Y           Z
 ---------------------------------------------------------------------
 ---------------------------------------------------------------------

 GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad
 Step number   1 out of a maximum of   2 on scan point  {key+1: 4d} out of  {noutput: 4d}

 GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad
                          Input orientation:                          
 ---------------------------------------------------------------------
 Center     Atomic      Atomic             Coordinates (Angstroms)
 Number     Number       Type             X           Y           Z
 ---------------------------------------------------------------------\n""",
*coordinates,
f""" ---------------------------------------------------------------------
 SCF Done:  E(UB3LYP) = {output.E: 15.12f}     A.U.

 GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad
 Step number   2 out of a maximum of   2 on scan point  {key+1: 4d} out of  {noutput: 4d}\n\n"""
            ]

        with open(path, "w") as f:
            f.writelines(header)
            f.writelines(lines)
            f.writelines(footer)

    def show(self):
        trash_path = f"/Users/takahirodoba/trash/{self.name}.log"
        self.to_gv(trash_path)
        subprocess.run(
            f"open -a /Applications/gv/gview.app {trash_path}", shell=True
        )

