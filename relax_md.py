from copy import deepcopy
from squid.structures import *
from squid.geometry import get_center_of_mass, packmol
from squid.forcefields.parameters import Parameters
from squid import files
from squid import lammps
from pathlib import Path
import click
import yaml

BASE_PATH = Path(__file__).parent
struct_path = BASE_PATH / "structs"
param_path = BASE_PATH / "params"
config_path = BASE_PATH / "config"
pack_path = BASE_PATH / "sys_packmol"
EPS = 1


def monolayer(num, coords, box, param):
    base_acryl = files.read_cml(str(struct_path / "acrylonitrile.cml"))[0]
    with open(str(param_path / f"acrylonitrile_{param}.labels")) as fobj:
        labels = map(int, fobj.read().split("\n"))
        for atom, label in zip(base_acryl.atoms, labels):
            atom.label = str(label)
    for _ in range(num):
        mol = deepcopy(base_acryl)
        box.add(mol)
    for a, atom in enumerate(box.atoms):
        atom.set_position(coords[a].flatten())
    return box


def solvent(num, coords, box, param):
    base_meth = files.read_cml(str(struct_path / "methane.cml"))[0]
    with open(str(param_path / f"methane_{param}.labels")) as fobj:
        labels = map(int, fobj.read().split("\n"))
        for atom, label in zip(base_meth.atoms, labels):
            atom.label = str(label)

    mols = [deepcopy(base_meth) for _ in range(num)]
    if len(coords) > 0:
        atoms = []
        for m in mols:
            box.add(m)
            for a in m.atoms:
                atoms.append(a)
        for a, c in zip(mols, coords):
            atom.set_position(c.flatten())
    else:
        max_x = box.xhi
        min_x = box.xlo

        max_y = box.yhi
        min_y = box.ylo

        # max_z = max([a.z for a in box.atoms]) - EPS
        # min_z = min([a.z for a in box.atoms]) + EPS
        max_z = 4.0
        min_z = -4.0
        constraint = f"outside box {min_x:.1f} {min_y:.1f} {min_z:.1f} {max_x:.1f} {max_y:.1f} {max_z:.1f}"
        packmol(
            box,
            [base_meth],
            number=num,
            additional=constraint,
            extra_block_at_end="nloop0 10000",
            persist=True,
            tolerance=2.0,
        )

    return box


def set_coords(mol_list, coords, name, param, box_dims, center=False):
    num1, num2 = mol_list
    mol_lens = [7, 4]
    atom_pops = [l * c for l, c in zip(mol_lens, [num1, num2])]
    xlen, ylen, zlen = [a - b for a, b in zip(*box_dims)]
    box = System(name, box_size=[xlen, ylen, zlen], periodic=True)
    acryl = monolayer(num1, coords[: atom_pops[0]], box, param)
    if center:
        cog = -1 * get_center_of_mass(coords[: atom_pops[0]])
        for a, atom in enumerate(acryl.atoms):
            atom.translate(cog)

    if coords[atom_pops[0]:]:
        meth = solvent(num2, coords[atom_pops[0] :], acryl, param)
        if center:
            for atom in meth.atoms[coords[atom_pops[0]] :]:
                atom.translate(cog)
        else:
            pass
    else:
        meth = solvent(num2, [], acryl, param)

    return meth


def write_input(input_script, cfg, sys):
    mol_lens = [7, 4]
    num1 = cfg["population"]["acrylonitrile"]
    num2 = cfg["population"]["methane"]
    atom_pops = [l * c for l, c in zip(mol_lens, [num1, num2])]

    input_str = yaml.safe_load(
        open(str(config_path / "inputs" / f"{input_script}.yaml"))
    )["job_str"]

    sys.set_types(opls_file=str(param_path / cfg["params"]))
        
    og_name = sys.name
    data_path = BASE_PATH / "data"
    sys.name = str(data_path / og_name)
    lammps.write_lammps_data(sys, pair_coeffs_included=False)

    xyz_name = str(struct_path / og_name)
    input_str = input_str.replace("$DATA$", sys.name + ".data")
    input_str = input_str.replace("$NAME$", xyz_name)
    input_str = input_str.replace("$PAIR_COEFFS$", sys.dump_pair_coeffs())
    input_str = input_str.replace("$ELEMENTS$", " ".join(sys.get_elements()))
    input_str = input_str.replace("$THERMO$", " ".join(cfg["thermo"]))
    input_str = input_str.replace("$VEL$", cfg["velocity"])
    input_str = input_str.replace("$ATOMS$", str(atom_pops[0]))
    sys.name = og_name

    return input_str


@click.command()
@click.option("--system", default="221_lpg")
@click.option("--name", default=None)
@click.option("--walltime", default="2:0:0")
@click.option("--cores", default=6)
@click.option("--debug", is_flag=True)
def create_membrane(system, name, walltime, cores, debug):
    systems_path = config_path / "systems"
    all_sys = [str(file.name) for file in systems_path.iterdir()]
    yaml_file = system + ".yaml"
    if yaml_file in all_sys:
        with open(str(systems_path / yaml_file), "r") as yobj:
            config = yaml.safe_load(yobj)
            coords = files.read_xyz(str(struct_path / config["xyz"]))
            counts = config["population"]
            num1, num2 = counts["acrylonitrile"], counts["methane"]
            if name == None:
                config["job_name"] = yaml_file.split(".")[0]
            else:
                config["job_name"] = name
            sys = set_coords(
                [num1, num2],
                coords,
                config["job_name"],
                config["params"].split(".")[0],
                box_dims=config["box_dims"],
                center=config["center"],
            )
            inp_str = write_input(config["sim"], config, sys)
            if debug:
                print(inp_str)
            else:
                lammps.job(
                    sys.name,
                    inp_str,
                    nprocs=cores,
                    walltime=walltime,
                    queue="defq",
                    pair_coeffs_in_data_file=False,
                    prebash="source /data/apps/go.sh\nml cl-lammps\n"
                )


if __name__ == "__main__":
    create_membrane()
