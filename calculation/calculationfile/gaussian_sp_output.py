import os
import re

from .output import Output


def read_gaussian_sp_output(
        path, pid=True, status=True, E=True, time=True):

    file = GaussianSpOutput()
    file.name = os.path.splitext(os.path.basename(path))[0]
    
    with open(path, "r") as f:
        lines = f.readlines()

    if pid:
        pid = [line for line in lines if "PID=" in line]

        if pid:
            file.pid = int(re.findall(r"PID=(.+).", pid[0])[0].strip())

    if status:
        normal_index = [
            i for i, line in enumerate(lines) if "Normal termination" in line
        ]
        error_index = [
            i for i, line in enumerate(lines) if "Error termination" in line
        ]

        if normal_index:
            file.status = "Normal termination"
        elif error_index:
            file.status = "Error termination"
        else:
            file.status = "Running"

    if E:
        E = [line for line in lines if "SCF Done:" in line]

        if E:
            file.opt.E = float(
                re.findall(r"=(.+)A.U.", E[0])[0].strip()
            )

    if time:
        time = [line for line in lines if "Job cpu time:" in line]

        if time:
            time = re.findall(
                "Job cpu time:(.+)days(.+)hours(.+)minutes(.+)seconds.",
                time[0]
            )
            time = [float(string.strip()) for string in time[0]]
            time = time[0] * 86400 + time[1] * 3600 + time[2] * 60 + time[3]
            file.time = time

    return file


class GaussianSpOutput:

    def __init__(
            self, name=None, pid=None, status=None, opt=None, time=None):

        self.name = name
        self.pid=pid
        self.status = status
        self.opt = opt
        self.time = time

        if not isinstance(self.opt, Output):
            self.opt = Output()

    @property
    def E(self):
        if self.status == "Normal termination":
            return self.opt.E
        else:
            return None

