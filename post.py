from supplements.time_energy_gen import read_log_or_out
from supplements.bin_tc_error import tb_P_inv
from squid.utils.units import convert_energy
from squid import files

import matplotlib.pyplot as plt
from pathlib import Path
import os.path as osp
import numpy as np
import click
import yaml
import os
import re
import pandas as pd

BASE_PATH = Path(__file__).parent
LAMMPS_DIR = BASE_PATH / "lammps"
STRUCT_DIR = BASE_PATH / "structs"
CONFIG_PATH = BASE_PATH / "config"

thermo_dix = {
    "step": "Step",
    "etotal": "TotEng",
    "ke": "KinEng",
    "pe": "PotEng",
    "temp": "Temp",
    "press": "Press",
}

def_idx = {
    "equil": [0, 1],
    "npt": [0, 1],
    "onlyz": [0, 1],
    "uncoupled": [0, 1],
    "zfirst": [0, 1]
}


def read_one_dir(
    job_name, thermo_format, indices=[0, 1], unit_conversion=None, plot_log=True
):
    job_dir = LAMMPS_DIR / job_name
    data_list = []
    fptr = str(job_dir / (job_name + ".log"))
    thermo_format = [thermo_dix[k] for k in thermo_format]
    with open(fptr) as fobj:
        column_indices = indices[:]
        step_and_vals = [
            [float(row[i]) for i in column_indices]
            for row in read_log_or_out(fobj, thermo_format)
        ]
    step, yvals = zip(*step_and_vals)

    if unit_conversion is not None and len(indices) == 2:
        yvals = [unit_conversion(yval) for yval in yvals]

    data_list = yvals[:]
    return step, data_list


def read_ti(job_base, thermo_format, nframes):
    pmf = []
    errors = []

    for i in range(nframes):
        dir_name = job_base + f"_{i}"
        data_list = read_one_dir(dir_name, thermo_format, indices=[2, 5])
        f1x, f2x = data_list[2], data_list[5]
        m_f1x, m_f2x = np.mean(f1x[-10000:]), np.mean(f2x[-10000:])
        s_f1x, s_f2x = np.std(f1x[-10000:]), np.std(f2x[-10000:])
        pmf_val = -1 * (m_f1x + m_f2x)
        error_val = np.sqrt(s_f1x**2 + s_f2x**2)

        # Conversion to kcal
        pmf_val = convert_energy("eV", "kcal", pmf_val)
        error_val = convert_energy("eV", "kcal", error_val)

        pmf.append(pmf_val)
        errors.append(error_val)

    return pmf, errors


def process(system, name=None, idx=None):
    systems_path = CONFIG_PATH / "systems"
    simulation = systems_path / (system + ".yaml")
    all_sys = list(systems_path.iterdir())
    if simulation in all_sys:
        yaml_file = simulation
        with open(yaml_file, "r") as y:
            config = yaml.safe_load(y)
            sim = config["sim"]
            thermo_format = config["thermo"]
            if name is None:
                job_name = system
            else:
                job_name = name
            if sim in ["equil", "freeze", "npt", "onlyz", "uncoupled", "zfirst"]:
                data_list = read_one_dir(
                    job_name,
                    thermo_format,
                    indices=[0, 1],
                    unit_conversion=None,
                )
                step, results = data_list
                return results, step
            elif sim == "ti":
                results, errors = read_ti(config["job_name"], thermo_format, nframes=23)
                return results, errors
            elif sim == "smd":
                data_list = read_one_dir(
                    config["job_name"], thermo_format, indices=[2, 3]
                )
                step, results = data_list
                # errors = []
                return results, step
            else:
                raise ValueError(f"Unsupported simulation type: {sim}")


if __name__ == "__main__":
    task = "221_lpg"
    process(task)
    results, steps = process(task)
    # print(len(results))
    mean = np.mean(results[-10000:])
    error = np.std(results[-10000:])
    print("Mean", mean)
    print("Error", error)
    annot = f"Mean:{mean:.2f}\nError:{error:.2f}"
    # _, _, t_c, _ = tb_P_inv(steps, results, range(10, 100, 10), 5, 9)
    # print(t_c)

    fig, ax = plt.subplots(1, 1)
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("System Energy (kcal/mol)")
    fig.suptitle(f"Moving average for {task}")
    window = 200
    results = np.convolve(results[8:], np.ones(window), "valid") / window
    steps = [t / 1000000 for t in steps[-len(results) :]]
    ax.plot(steps, results)
    ax.text(
        0.7,
        0.5,
        annot,
        color="black",
        transform=ax.transAxes,
        verticalalignment="top",
        fontsize=12,
        bbox=dict(facecolor="wheat", edgecolor="black", boxstyle="round,pad=1"),
    )
    fig.savefig(str(BASE_PATH / "plots" / f"{task}_energies.png"))
