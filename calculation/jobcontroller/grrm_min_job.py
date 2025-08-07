import glob
import os
import subprocess
import time

import calculationfile as cf


class GrrmMinJob:

    def __init__(self, wd, input_, output):
        """
        self.wd : str
            Working directory. e.g. "."
        input_ : str
            The name of the input file with extension. e.g. "td_c1_100.com"
        """
        self.wd = wd
        self.input_ = input_
        self.output = output

    def submit(self, option=""):
        subprocess.run(
            [f"GRRMs {os.path.splitext(self.input_)[0]} {option}"],
            capture_output=True,
            shell=True,
            cwd=self.wd,
            check=True,
            text=True
        )
        print(f"{self.wd}/{self.input_}: Submitted.")

    def wait(self, walltime=172800):

        for _ in range(walltime):
            
            if f"{self.wd}/{self.output}" in glob.glob(f"{self.wd}/*"):
                file = cf.read_grrm_min_output(
                    f"{self.wd}/{self.output}",
                    status=True,
                    itrs=False,
                    opt=False,
                    hessian=False,
                    freq=False,
                    thermo=False,
                    time=False
                )

                if file.status == "Minimum point was found":
                    print(f"{self.wd}/{self.output}: Minimum point was found.")

                elif file.status == "The structure is dissociating":
                    print(f"{self.wd}/{self.output}: The structure is dissociating.")

                elif file.status == "Error termination":
                    print(f"{self.wd}/{self.output}: Error termination.")

                else:
                    print(f"{self.wd}/{self.output}: Unexpected error.")

                break

            else:
                pass

            time.sleep(1)
