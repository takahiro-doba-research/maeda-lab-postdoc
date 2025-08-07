import os
import re

import numpy as np
import pandas as pd

from .coordinate import Coordinate
from .output import Output


def read_grrm_min_output(
        path, status=True, itrs=True, opt=True,
        hessian=True, freq=True, thermo=True, time=True):
    
    file = GrrmMinOutput()
    file.name = os.path.splitext(os.path.basename(path))[0]

    with open(path, "r") as f:
        lines = f.readlines()

    if status:
        minimum_index = [
            i for i, line in enumerate(lines)
            if "Minimum point was found" in line
        ]
        dissociating_index = [
            i for i, line in enumerate(lines)
            if "The structure is dissociating" in line
        ]
        saddle_index = [
            line for line in lines
            if "Saddle point was found" in line
        ]

        if minimum_index:
            file.status = "Minimum point was found"
        elif dissociating_index:
            file.status = "The structure is dissociating"
        elif saddle_index:
            file.status = "Saddle point was found"
        else:
            file.status = "Error termination"

    if itrs:
        itr_indices = [i for i, line in enumerate(lines) if "#" in line]
        blank_indices = [i for i, line in enumerate(lines) if line == "\n"]

        if itr_indices:
            itrs = []

            for i in range(len(itr_indices)):
                itr = Output()
                stop_index = [
                    index for index in blank_indices
                    if index > itr_indices[i]
                ]
                itr_block = lines[itr_indices[i]:stop_index[0]]

                # Set Output.name
                itr.name = int(
                    re.findall(r"# ITR.(.+)", itr_block[0])[0].strip()
                )

                # Set Output.coordinate
                item_index = [
                    i for i, line in enumerate(itr_block) if "Item" in line
                ]
                coordinate = itr_block[1:item_index[0]]
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
                itr.coordinate = Coordinate(df)

                # Set Output.E
                E = [line for line in itr_block if "ENERGY" in line]
                itr.E = float(
                    re.findall(r"ENERGY(.+)\(", E[0])[0].strip()
                )

                # Set Output.spin
                spin = [line for line in itr_block if "Spin(**2)" in line]
                itr.spin = float(
                    re.findall(r"\)(.+)", spin[0])[0].strip()
                )

                # Set Output.normal_mode
                normal_mode_index = [
                    i for i, line in enumerate(itr_block)
                    if "NORMAL MODE EIGENVALUE" in line
                ]
                normal_mode = itr_block[normal_mode_index[0]+1:]
                normal_mode = [
                    float(value)
                    for line in normal_mode
                    for value in line.strip().split()
                ]
                itr.normal_mode = np.array(normal_mode)

                # Append the complete itr
                itrs.append(itr)

            # Set GRRMMINFile.itrs
            file.itrs = itrs

    if opt:
        start_index = [
            i for i, line in enumerate(lines)
            if "Optimized structure" in line
        ]
        stop_index = [
            i for i, line in enumerate(lines)
            if "was found" in line
        ]

        if stop_index:
            opt = Output()
            opt_block = lines[start_index[0]+1:stop_index[0]-1]

            # Set Output.df
            E_index = [
                i for i, line in enumerate(opt_block) if "ENERGY" in line
            ]
            coordinate = opt_block[0:E_index[0]]
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
            opt.coordinate = Coordinate(df)

            # Set Output.E
            opt.E = float(
                re.findall(r"=(.+)", opt_block[E_index[0]])[0].strip()
            )

            # Set Output.spin
            spin = [line for line in opt_block if "Spin(**2)" in line]
            opt.spin = float(re.findall(r"=(.+)", spin[0])[0].strip())

            # Set Output.E_zpv
            E_zpv = [line for line in opt_block if "ZPVE" in line]
            opt.E_zpv = float(re.findall(r"=(.+)", E_zpv[0])[0].strip())

            # Set Output.gradient
            gradient_index = [
                i for i, line in enumerate(opt_block)
                if "GRADIENT VECTOR" in line
            ]
            hessian_index = [
                i for i, line in enumerate(opt_block)
                if "HESSIAN MATRIX" in line
            ]
            gradient = [
                float(line.strip())
                for line in opt_block[gradient_index[0]+1:hessian_index[0]]
            ]
            opt.gradient = np.array(gradient)

            # Set Output.normal_mode
            normal_mode_index = [
                i for i, line in enumerate(opt_block)
                if "NORMAL MODE EIGENVALUE" in line
            ]
            normal_mode = opt_block[normal_mode_index[0]+1:]
            normal_mode = [
                float(value)
                for line in normal_mode
                for value in line.strip().split()
            ]
            opt.normal_mode = np.array(normal_mode)

            # Set Output.hessian
            if hessian:
                hessian = opt_block[hessian_index[0]+1:normal_mode_index[0]]
                hessian = [
                    [float(value) for value in line.strip().split()]
                    for line in hessian
                ]
                head_indices = [
                    i for i, line in enumerate(hessian) if len(line) == 1
                ]

                if len(head_indices) == 1:
                    shape0 = len(hessian)
                    shape1 = len(hessian[-1])
                    hessian = [
                        line + [0] * (shape1 - len(line))
                        for line in hessian
                    ]
                    opt.hessian = np.array(hessian)

                else:
                    shape0 = head_indices[1] - head_indices[0]
                    hessian_blocks = []

                    for i in range(len(head_indices)):

                        if i != len(head_indices) - 1:
                            hessian_block = hessian[
                                head_indices[i]:head_indices[i+1]
                            ]
                        else:
                            hessian_block = hessian[
                                head_indices[i]:normal_mode_index[0]
                            ]

                        shape1 = len(hessian_block[-1])
                        zeros = np.zeros((shape0-len(hessian_block), shape1))
                        hessian_block = [
                            line + [0] * (shape1 - len(line))
                            for line in hessian_block
                        ]
                        hessian_block = np.array(hessian_block)
                        hessian_block = np.concatenate([zeros, hessian_block])
                        hessian_blocks.append(hessian_block)

                    opt.hessian = np.concatenate(hessian_blocks, 1)

            # Set thermochemistry
            if thermo:
                thermo_indices = [
                    i for i, line in enumerate(lines)
                    if "Thermochemistry" in line
                ]
                before = [
                    float(
                        re.findall(
                            r"=(.+)", lines[thermo_indices[0]+i+2]
                        )[0].strip().split()[0]
                    )
                    for i in range(13)
                ]
                after = [
                    float(
                        re.findall(
                            r"=(.+)", lines[thermo_indices[1]+i+2]
                        )[0].strip().split()[0]
                    )
                    for i in range(13)
                ]
                opt.E_zpv = before[0]
                opt.E_E_zpv = before[1]
                opt.E_tr = before[2]
                opt.E_rot = before[3]
                opt.E_vib = before[4]
                opt.H_corr = before[5]
                opt.H = before[6]
                opt.S_el = before[7]
                opt.S_tr = before[8]
                opt.S_rot = before[9]
                opt.S_vib = before[10]
                opt.G_corr = before[11]
                opt.G = before[12]

                opt.E_zpv_r = after[0]
                opt.E_E_zpv_r = after[1]
                opt.E_tr_r = after[2]
                opt.E_rot_r = after[3]
                opt.E_vib_r = after[4]
                opt.H_corr_r = after[5]
                opt.H_r = after[6]
                opt.S_el_r = after[7]
                opt.S_tr_r = after[8]
                opt.S_rot_r = after[9]
                opt.S_vib_r = after[10]
                opt.G_corr_r = after[11]
                opt.G_r = after[12]

            # Set frequencies
            if freq:
                freq = [line for line in lines if "Freq." in line]
                freq = [
                    float(value)
                    for line in freq
                    for value in re.findall(r":(.+)", line)[0].strip().split()
                ]
                opt.freq = np.array(freq)

            # Set GRRMMINFile.opt
            file.opt = opt

    if time:
        time = [line for line in lines if "TOTAL ELAPSED TIME" in line]

        if time:
            file.time = float(re.findall(r":(.+)SEC.", time[0])[0].strip())

    return file


class GrrmMinOutput:

    def __init__(
            self, name=None, status=None, itrs=None, opt=None, time=None):
            
        self.name = name
        self.status = status
        self.itrs = itrs or []
        self.opt = opt
        self.time = time

        if not isinstance(self.opt, Output):
            self.opt = Output()

    @property
    def coordinate(self):
        return self.opt.coordinate

    @property
    def E(self):
        return self.opt.E

    @property
    def freq(self):
        return self.opt.freq
