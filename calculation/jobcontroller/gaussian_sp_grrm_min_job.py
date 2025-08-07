import glob
import os
import re

import calculationfile as cf

from .gaussian_sp_job import GaussianSpJob
from .grrm_min_job import GrrmMinJob


def gaussian_sp_grrm_min_job(wd, sp_file, sp_option, min_file, min_option):
    paths = sorted(
        glob.glob(f"{wd}/*.com"),
        key=lambda x:os.path.splitext(os.path.basename(x))[0]
    )

    if len(paths) == 1:
        file = cf.read_gaussian_sp_input(paths[0])
        file.coordinate.save(
            f"{wd}/{file.name}_MIN.com",
            header=min_file.header[:-1]+file.header[-1:],
            footer=min_file.footer,
            ignore_notes=True
        )
        grrm_min_job = GrrmMinJob(
            wd, f"{file.name}_MIN.com", f"{file.name}_MIN.log"
        )
        grrm_min_job.submit(min_option)
        grrm_min_job.wait()

    else:

        for path in paths:
            file = cf.read_gaussian_sp_input(path)
            file.coordinate.save(
                f"{wd}/{file.name}_SP.com",
                header=sp_file.header[:-1]+file.header[-1:],
                footer=sp_file.footer,
                ignore_notes=True
            )
            gaussian_sp_job = GaussianSpJob(
                wd, f"{file.name}_SP.com", f"{file.name}_SP.log"
            )
            gaussian_sp_job.submit(sp_option)
            gaussian_sp_job.wait()

        output_paths = sorted(glob.glob(f"{wd}/*_SP.log"))
        """
        Python 3.8 -
        output_files = [
            file for output_path in output_paths
            if (file := cf.read_gaussian_sp_output(
                output_path, pid=False, status=True, E=True, time=False
            )).status == "Normal termination"
        ]
        """
        output_files = []

        for output_path in output_paths:
            file = cf.read_gaussian_sp_output(
                output_path, pid=False, status=True, E=True, time=False
            )

            if file.status == "Normal termination":
                output_files.append(file)

        if output_files:
            lowest_energy_file = min(output_files, key=lambda x:x.E)
            name = re.findall(r"(.+)_SP", lowest_energy_file.name)[0]
            file = cf.read_gaussian_sp_input(f"{wd}/{name}.com")
            file.coordinate.save(
                f"{wd}/{name}_MIN.com",
                header=min_file.header[:-1]+file.header[-1:],
                footer=min_file.footer,
                ignore_notes=True
            )
            grrm_min_job = GrrmMinJob(wd, f"{name}_MIN.com", f"{name}_MIN.log")
            grrm_min_job.submit(min_option)
            grrm_min_job.wait()
