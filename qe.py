import hydra
import os
import time
from omegaconf import DictConfig, OmegaConf
from squid import files, structures
from squid.jobs import slurm
from pathlib import Path
from pymatgen.core.lattice import Lattice
from squid.utils.units import convert_dist
from squid.geometry import get_center_of_mass
import numpy as np
import shutil

BASE_PATH = Path(__file__).parent.resolve()
xyz_path = BASE_PATH / "structs"
data_path = BASE_PATH / "data"
param_path = BASE_PATH / "params"
config_path = BASE_PATH / "config"


def qe_cell_dims(lengths, sg, angles):
    lengths = [a - b for a, b in zip(*lengths)]
    if sg == 0:
        lat = Lattice.from_parameters(
            a=lengths[0], b=lengths[1], c=lengths[2], alpha=angles[0],
            beta=angles[1], gamma=angles[2]
        )
        celldm_str = []
        for line in lat.matrix.tolist():
            line_str = [i if i > 1e-8 else 0.0 for i in line]
            line_str = list(map(lambda x: f"{x:<12f}", line_str))
            celldm_str.append(" ".join(line_str))
        return "\n".join(celldm_str)
    else:
        celldm = [None for _ in range(6)]
        lengths = [convert_dist("Angstrom", "Bohr", l) for l in lengths]
        angles = [a * np.pi / 180 for a in angles]
        celldm[0] = lengths[0]
        celldm[1] = lengths[1] / lengths[0]
        celldm[2] = lengths[2] / lengths[0]
        if sg == 8:
            pass
        elif sg == 14:
            for i in range(3):
                celldm[i + 3] = np.cos(angles[i])
        else:
            raise KeyError("This Quantum Epsresso ibrav val is not yet implemented")
        celldm_str = []
        for i, dm in enumerate(celldm):
            if dm != None:
                celldm_str.append(f"  celldm({str(i + 1)}) = {dm:.2f}")
        return "\n" + "\n".join(celldm_str)



def dump_atoms_data(system):
    atom_dix = {"H": 1, "C": 2, "N": 3}
    atoms = files.read_xyz(str(xyz_path / system))
    com = -1 * get_center_of_mass(atoms)
    atoms_data = []
    for a in atoms:
        a.translate(com)
        label = atom_dix[a.element]
        atoms_data.append(f"{label}\t{a.x: 12.10f}\t{a.y: 12.10f}\t{a.z: 12.10f}")
    pos_str = "\n".join(atoms_data)

    return pos_str, len(atoms)


def write_input_script(inp_base, job_name, sys, restart):
    dims = [str((high - low)) for high, low in zip(*sys.box_dims)]
    dims = " ".join(dims)
    inp_str = inp_base.replace("$DIRNAME$", job_name)
    inp_str = inp_str.replace("$QEDIR$", str(BASE_PATH / "qe"))
    inp_str = inp_str.replace("$PRMDIR$", str(param_path))
    inp_str = inp_str.replace("$IBRAV$", str(sys["cell_shape"]["sg"]))
    if not restart:
        inp_str = inp_str.replace("$NAME$", job_name)
        inp_str = inp_str.replace("$START$", "from_scratch")
    else:
        jobpath = BASE_PATH / "qe" / job_name
        outfiles = [x for x in jobpath.iterdir() if str(x).split(".")[-1] == "out"]
        # inp_str = inp_str.replace("$NAME$", f"{job_name}{len(outfiles)}")
        inp_str = inp_str.replace("$NAME$", f"{job_name}")
        inp_str = inp_str.replace("$START$", "restart")
    pos, nat = dump_atoms_data(sys["xyz"])
    celldim = qe_cell_dims(
        sys["box_dims"], sys["cell_shape"]["sg"], sys["cell_shape"]["angles"]
    )
    inp_str = inp_str.replace("$POSITIONS$", pos)
    inp_str = inp_str.replace("$NUM$", str(nat))

    if sys["cell_shape"]["sg"] != 0:
        inp_str = inp_str.replace("$DIMENSIONS$", celldim)
    else:
        inp_str = inp_str.replace("$DIMENSIONS$", "")
        inp_str += "\nCELL_PARAMETERS angstrom\n"
        inp_str += celldim
    return inp_str


def create_job_dir(job_name, inp_str, clean=True):
    inp_dir = BASE_PATH / "qe" / job_name
    try:
        inp_dir.mkdir(parents=True, exist_ok=clean)
    except FileExistsError:
        for f in inp_dir.iterdir():
            if f.is_dir():
                shutil.rmtree(f)
            else:
                f.unlink()

    inp_file = inp_dir / f"{job_name}.in"
    with open(inp_file, "w") as f:
        f.write(inp_str)
    return inp_dir, inp_file


def write_slurm_inp(job_name, inp_path, cores, k):
    out_path = inp_path.parent
    outfiles = [x for x in out_path.iterdir() if str(x).split(".")[-1] == "out"]
    # out_path = str(out_path / f"{job_name}{len(outfiles)}.out")
    out_path = str(out_path / f"{job_name}.out")

    slurm_inp = f"""
date
ml qe/7.2
mpirun -np {str(cores)} pw.x -npool {str(k)} -in {inp_path} > {out_path}
date
"""
    return slurm_inp


@hydra.main(version_base=None, config_path=str(config_path), config_name="qe")
@hydra.main(version_base=None, config_path=str(config_path), config_name="qe")
@hydra.main(version_base=None, config_path=str(config_path), config_name="qe")
def main(cfg: DictConfig) -> None:
    inp, sys, job_name = cfg.inputs, cfg.systems, cfg.job_name
    inp_str = write_input_script(inp.job_str, job_name, sys, restart=cfg.restart)
    inp_str = inp_str.replace("$KPOINTS$", " ".join(cfg.kpoints.split("-")))
    hms = [int(i) for i in cfg.walltime.split(":")]
    conversion = [3600, 60, 1]
    time = sum([i * j for i, j in zip(hms, conversion)])
    inp_str = inp_str.replace("$TIME$", str(max(10000, time - 10000)))
    inp_dir, inp_fptr = create_job_dir(job_name, inp_str, cfg.restart)
    slurm_inp = write_slurm_inp(job_name, inp_fptr, cfg.cores, cfg.kpoints.split("-")[0])

    print(f"Submitting job..... {job_name}")
    print("-----------------")
    # print(inp_str)
    # raise Exception
    os.chdir(inp_dir)
    slurm.submit_job(
        job_name,
        slurm_inp,
        queue="parallel",
        allocation="pclancy3",
        ntasks=cfg.cores,
        walltime=cfg.walltime,
    )
    os.chdir(BASE_PATH)


if __name__ == "__main__":
    main()
