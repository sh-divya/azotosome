import hydra
import os
import time
from omegaconf import DictConfig, OmegaConf
from squid import files, structures
from squid.jobs import slurm
from pathlib import Path

BASE_PATH = Path(__file__).parent.resolve()
xyz_path = BASE_PATH / "structs"
data_path = BASE_PATH / "data"
param_path = BASE_PATH / "params"
config_path = BASE_PATH / "config"


def write_input_script(inp_base, job_name, sys, restart):
    dims = [str((high - low)) for high, low in zip(*sys.box_dims)]
    dims = " ".join(dims)
    inp_str = inp_base.replace("$JOB_NAME$", job_name)
    inp_str = inp_str.replace("$PRM_PATH$", str(param_path))
    if restart:
        inp_str = inp_str.replace("$XYZ_FILE$", str(xyz_path / sys.xyz))
        inp_str = inp_str.replace("$DIMENSIONS$", dims)
    else:
        inp_str = inp_str.replace("$RESTART$", sys.xyz)
    return inp_str


def create_job_dir(job_name, inp_str, clean=True):
    inp_dir = BASE_PATH / "cp2k" / job_name
    try:
        inp_dir.mkdir(parents=True, exist_ok=not clean)
    except FileExistsError:
        for f in inp_dir.iterdir():
            f.unlink()

    inp_file = inp_dir / f"{job_name}.in"
    with open(inp_file, "w") as f:
        f.write(inp_str)
    return inp_dir, inp_file


def write_slurm_inp(job_name, inp_path, cores):
    out_path = str(inp_path).split(".")[0] + ".out"
    slurm_inp = f"""
date
ml GCC/12.3.0 openmpi/4.1.4 cp2k
mpirun -np {str(cores)} cp2k.popt -in {inp_path} > {out_path}
date
"""
    return slurm_inp


@hydra.main(version_base=None, config_path=str(config_path), config_name="aimd")
def main(cfg: DictConfig) -> None:
    inp, sys, job_name = cfg.inputs, cfg.systems, cfg.job_name
    inp_str = write_input_script(inp.job_str, job_name, sys, restart=cfg.restart)
    inp_dir, inp_fptr = create_job_dir(job_name, inp_str, cfg.restart)
    slurm_inp = write_slurm_inp(job_name, inp_fptr, cfg.cores)

    print(f"Submitting job..... {job_name}")
    print("-----------------")
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
