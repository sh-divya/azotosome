from supplements.time_energy_gen import read_log_or_out
from supplements.bin_tc_error import tb_P_inv
from squid.utils.units import convert_energy
from squid import files, geometry

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
    "c_mempe": "c_mempe",
}

def_idx = {
    "equil": [0, 1],
    "npt": [0, 1],
    "onlyz": [0, 1],
    "uncoupled": [0, 1],
    "zfirst": [0, 1],
    "eq_pe": [0, 1],
    "freeze": [0, 1],
    "nvt": [0, 1],
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
        step = []
        yvals = []
        for row in read_log_or_out(fobj, thermo_format):
            step.append(float(row[0]))
            if len(column_indices) == 2:
                yvals.append(float(row[1]))
            else:
                yvals.append([float(row[i]) for i in column_indices[1:]])
    #     step_and_vals = [
    #         [float(row[i]) for i in column_indices]
    #         for row in read_log_or_out(fobj, thermo_format)
    #     ]
    # step, yvals = zip(*step_and_vals)

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


def process(system, task="log", name=None, idx=None):
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
            if sim in [
                "equil",
                "freeze",
                "npt",
                "onlyz",
                "uncoupled",
                "zfirst",
                "eq_pe",
                "nvt",
            ]:
                if task == "log":
                    if not idx:
                        idx = def_idx[sim]
                    data_list = read_one_dir(
                        job_name,
                        thermo_format,
                        indices=idx,
                        unit_conversion=None,
                    )
                    cols = [thermo_format[c] for c in idx]
                    (
                        step,
                        results,
                    ) = data_list
                    return results, step, cols
                elif task == "geom":
                    ensemble_angles(
                        system,
                        list(range(450, 550)),
                        config["population"]["acrylonitrile"],
                        config["box_dims"],
                    )
            elif sim == "ti":
                results, errors = read_ti(config["job_name"], thermo_format, nframes=23)
                return results, errors, None
            elif sim == "smd":
                data_list = read_one_dir(
                    config["job_name"], thermo_format, indices=[2, 3]
                )
                step, results = data_list
                return results, step, None
            else:
                raise ValueError(f"Unsupported simulation type: {sim}")


def sim_summary(xvals, yvals):
    mean = np.mean(yvals[-10000:])
    error = np.std(yvals[-10000:])
    # _, _, t_c, _ = tb_P_inv(steps, results, range(10, 100, 10), 5, 9)
    # print(t_c)

    return mean, error


def plot_results(xvals, yvals, name, cols, units=None, save=True, fig_ax=None):
    window = 200
    steps = np.convolve(xvals[8:], np.ones(window), "valid") / (window * 1000000)
    x = thermo_dix[cols[0]]
    try:
        yvals = list(zip(*yvals))
        if fig_ax is None:
            fig, ax = plt.subplots(len(yvals), 1, figsize=[15, 12])
        else:
            fig, ax = fig_ax
        fig.supxlabel(x)
        for i, y in enumerate(yvals):
            res = np.convolve(y[8:], np.ones(window), "valid") / window
            ax[i].plot(steps[-len(res) :], res)
            ax[i].set_ylabel(thermo_dix[cols[1:][i]])
    except TypeError:
        yvals = yvals[8:]
        if fig_ax is None:
            fig, ax = plt.subplots(1, 1)
        else:
            fig, ax = fig_ax
        res = np.convolve(yvals, np.ones(window), "valid") / window
        ax.plot(steps[-len(res) :], res)
        fig.supxlabel(x)
        ax.set_ylabel(thermo_dix[cols[1]])

    if save:
        fig.suptitle(f"Moving average for {name}")
        fig.tight_layout()
        fig.savefig(str(BASE_PATH / "plots" / f"{name}_energies.png"))
    else:
        return fig, ax


def get_mvee_orientation(frame):
    acryls = [frame[i : i + 7] for i in range(0, len(frame), 7)]
    centers = []
    radii = []
    for mol in acryls:
        A, c = geometry.mvee(mol)
        _, Q, _ = np.linalg.svd(A)
        r = []
        for i in range(3):
            r.append(1 / (Q[i] ** 0.5))
        radii.append(r)
        centers.append(c)

    return centers, radii


def get_vec_orientation(frame, ang_idx=[2, 3]):
    acryls = [frame[i : i + 7] for i in range(0, len(frame), 7)]
    centers = []
    radii = []
    for mol in acryls:
        c = geometry.get_center_of_mass(mol)
        # print(c)
        vec1 = mol[ang_idx[0]].flatten()
        vec2 = mol[ang_idx[1]].flatten()
        r = vec2 - vec1
        # print(r)

        radii.append(r)
        centers.append(c)
    print(centers)
    print(radii)
    return centers, radii


def ensemble_angles(name, frame_idx, mols, box_dims, save=True):
    xyz_file = STRUCT_DIR / f"{name}.xyz"
    all_frames = files.read_xyz(str(xyz_file))
    frames = [all_frames[idx] for idx in frame_idx]
    box_size = [high - low for high, low in zip(*box_dims)]
    centers = []
    radii = []
    for frame in frames:
        # centers, radii = get_mvee_orientation(frame[: 7 * mols])
        c, r = get_vec_orientation(frame[: 7 * mols])
        centers.append(c)
        radii.append(r)
    centers = np.mean(centers, axis=0)
    radii = np.mean(radii, axis=0)
    fig, ax = plot_sys(centers, radii, box_size)

    if save:
        fig.savefig(BASE_PATH / "plots" / f"{name}_angles.png")
    else:
        return fig, ax



def plot_sys(points, vectors, size):
    print(points)
    print(vectors)
    points = [p[:2] for p in points]
    x, y = zip(*points)
    # for i, p in enumerate(points):
    #     print(p, x[i], y[i])
    vectors = [v[:2] for v in vectors]
    u, v = zip(*vectors)
    
    fig, ax = plt.subplots(1, 1)
    ax.quiver(list(x), list(y), list(u), list(v))
    return fig, ax


def compare_runs(xvals, yvals, sims, cols):
    labels = ["Staggered", "Flat"]
    fig_ax = plt.subplots(1, 1)
    for i, (x, y) in enumerate(zip(xvals, yvals)):
        plot_results(x, y, sims[i], cols, save=False, fig_ax=fig_ax)
    fig, ax = fig_ax
    for line, label in zip(ax.lines, labels):
        line.set_label(label)
    ax.legend()
    fig.savefig(str(BASE_PATH / "plots" / "temp1.png"))


if __name__ == "__main__":
    task1 = "221_js2"
    # process(task, "geom")
    # results, steps, cols = process(task, idx=[0, 1, 2, 3, 4])
    # results, steps, cols = process(task, idx=[0, 1, 2, 3])
    # results, steps, cols = process(task, idx=[0, 1, 4, 5])
    results1, steps1, cols1 = process(task1, idx=[0, 1])
    task2 = "221_js2_flat"
    results2, steps2, cols2 = process(task2, idx=[0, 1])
    # annot = f"Mean:{mean:.2f}\nError:{error:.2f}"
    compare_runs([steps1, steps2], [results1, results2], [task1, task2], cols1)
    # plot_results(steps, results, task, cols)
