from supplements.time_energy_gen import read_log_or_out
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
NEB_DIR = BASE_PATH / "neb"


def read_one_dir(
    job_name, thermo_format, indices=[0, 1], unit_conversion=None, plot_log=True
):
    job_dir = osp.join(LAMMPS_DIR, job_name)
    if plot_log:
        all_out_fptr = [f for f in os.listdir(job_dir) if f.endswith(".log")]
    else:
        all_out_fptr = [f for f in os.listdir(job_dir) if re.match(r".+\.o\d+", f)]

    data_list = []

    for f, fptr in enumerate(all_out_fptr):
        with open(osp.join(job_dir, fptr)) as fobj:
            column_indices = indices[:]
            step_and_vals = [
                [float(row[i]) for i in column_indices]
                for row in read_log_or_out(fobj, thermo_format)
            ]
        step, yvals = zip(*step_and_vals)

        if unit_conversion is not None and len(indices) == 2:
            yvals = [unit_conversion(yval) for yval in yvals]

        data_list = yvals[:]
    return data_list


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


def read_single_neb():
    pass


def neb_sweep(config, filter=True):
    base_cols = [
        "Step",
        "MaxReplicaForce",
        "MaxAtomForce",
        "GradV0",
        "GradV1",
        "GradVc",
        "EBF",
        "EBR",
        "RDT",
    ]

    sweep = yaml.safe_load(open(str(CONFIG_PATH / "neb" / "sweep.yaml")))
    frames = sweep["nframes"]

    job_base = config["job_name"] + "_sweep"
    mech = config["job_name"]

    plt.figure()
    plt.title(mech)
    ecols = [f"E{i}" for i in range(frames)]
    rcols = [f"R{i}" for i in range(frames)]
    cols = ecols + rcols
    master_df = pd.DataFrame(columns=cols)
    for i, d in enumerate(os.listdir(str(NEB_DIR))):
        sim_out = []
        if d.startswith(job_base):
            job_dir = NEB_DIR / d
            log_file = job_dir / "log.lammps"
            for line in open(str(log_file)).readlines()[3:]:
                row = line.split()
                if row[0] == "Climbing" or row[0] == "Step":
                    continue
                data = [float(r) for r in row]
                sim_out.append(data)
            df = pd.DataFrame(np.array(sim_out))

            extra_cols = [
                "RD" + str(f // 2 + 1) if f % 2 == 0 else "PE" + str(f // 2 + 1)
                for f in range(2 * frames)
            ]

            all_cols = base_cols + extra_cols

            df.columns = all_cols

            coords = [extra_cols[i] for i in range(len(extra_cols)) if i % 2 == 0]
            energies = [extra_cols[i] for i in range(len(extra_cols)) if i % 2 != 0]
            rel_energies = (df[energies].iloc[-1] - df["PE1"].iloc[-1]).to_numpy(
                dtype=float
            )
            if filter:
                suff = "filter"
                if np.amin(rel_energies) < 0.0:
                    continue
                elif rel_energies[-1] > rel_energies[0] + 0.05:
                    continue
                elif np.amax(rel_energies) >= 2.0:
                    continue
            else:
                suff = "all"
            rc = df[coords].iloc[-1].to_numpy(dtype=float)
            master_df.loc[i] = list(rel_energies) + list(rc)
            plt.plot(rc, rel_energies)
    plt.ylabel("Energy eV")
    plt.xlabel("Reaction Coordinate")
    plt.savefig(str(BASE_PATH / "plots" / f"{config['job_name']}_{suff}.png"))

    y = master_df[ecols].mean()
    x = master_df[rcols].mean()
    errs = master_df[ecols].std()
    print(y.tolist())
    print(x.tolist())
    print(errs.to_list())

    plt.figure()
    plt.title(mech)
    plt.ylabel("Energy eV")
    plt.xlabel("Reaction Coordinate")
    plt.ylim([0.0, 1.5])
    plt.scatter(x, y)
    plt.errorbar(x, y, yerr=errs, capsize=3.5)
    plt.savefig(str(BASE_PATH / "plots" / f"{config['job_name']}_{suff}_avg.png"))

    return master_df.mean(), master_df.std()


def process(task, sim="md"):
    systems_path = CONFIG_PATH / "systems" / sim
    simulation = systems_path / (task + ".yaml")
    all_sys = list(systems_path.iterdir())
    if simulation in all_sys:
        yaml_file = simulation
        with open(yaml_file, "r") as y:
            config = yaml.safe_load(y)
            sim = config["sim"]
            thermo_format = config["output"]

            if sim == "md" or sim == "single":
                data_list = read_one_dir(
                    config["job_name"],
                    thermo_format,
                    indices=[0, 1],
                    unit_conversion=lambda x: convert_energy("eV", "kcal", x),
                )
                results = data_list
                errors = []
            elif sim == "ti":
                results, errors = read_ti(config["job_name"], thermo_format, nframes=23)
            elif sim == "smd":
                data_list = read_one_dir(
                    config["job_name"], thermo_format, indices=[2, 3]
                )
                results = data_list
                errors = []
            elif sim == "neb":
                neb_sweep(config, False)
                return
            else:
                raise ValueError(f"Unsupported simulation type: {sim}")

        return results, errors


if __name__ == "__main__":
    # task = "md_smaller_interface"
    # results, errors = process(task)
    # print("Mean", np.mean(results[-1000:]))
    # print("Error", np.std(results[-1000:]))
    task = "ti_100_ex"
    process(task, "neb")
