import glob
import os
import subprocess
import time

import calculationfile as cf
    

class GaussianSpJob:

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
            [f'bsub {option} "g16 {self.input_}"'],
            capture_output=True,
            shell=True,
            cwd=self.wd,
            check=True,
            text=True
        )
        print(f"{self.wd}/{self.input_}: Submitted.")

    def wait(self, walltime=86400):

        for _ in range(walltime):
            
            if f"{self.wd}/{self.output}" in glob.glob(f"{self.wd}/*"):
                file = cf.read_gaussian_sp_output(
                    f"{self.wd}/{self.output}",
                    pid=False,
                    status=True,
                    E=False,
                    time=False
                )
                
                if file.status == "Normal termination":
                    print(f"{self.wd}/{self.output}: Normal termination.")
                    break

                elif file.status == "Error termination":
                    print(f"{self.wd}/{self.output}: Error termination.")
                    break

                elif file.status == "Running":
                    pass

                else:
                    print(f"{self.wd}/{self.output}: Unexpected error.")
                    break

            else:
                pass

            time.sleep(1)

