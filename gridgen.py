from pathlib import Path
from copy import deepcopy
import yaml

import numpy as np
from squid import files

from sys_prep import box_info, make_supercell, place_acryls, conv_density
from relax_md import base_create_membrane

BASE_PATH = Path(__file__).parent
struct_path = BASE_PATH / "structs"
systems_path = BASE_PATH / "config/systems"

min_deltax = 1.5
max_deltax = 4.0
xvals = np.linspace(min_deltax, max_deltax + 0.5, 6, endpoint=False)

min_deltay = 2.5
max_deltay = 4.0
yvals = np.linspace(min_deltay, max_deltay + 0.5, 4, endpoint=False)

prefix = "stag_pack"

methane = files.read_xyz(str(struct_path / "tianle_meth.xyz"))
base_yaml = yaml.safe_load(open(systems_path / "221_js2.yaml"))

# for x in xvals:
#     for y in yvals:
for x in [4.0]:
    for y in [3.5]:
        xstr = str(int(x * 10))
        ystr = str(int(y * 10))
        name = f"{prefix}x{xstr}y{ystr}"
        unit = [place_acryls(4, distance=[x, y, -2.058], write=False)]
        solvent = [deepcopy(methane), deepcopy(methane)]
        # unit.extend(solvent)
        cell, name = make_supercell(
            name,
            unit,
            [2, 2, 1],
            [x, y, 0],
            offset=[[0, 0, 0]],  # , [3, 0, 7], [3, 0, -5]],
            write=True,
        )

        a, b, c = box_info(cell)
        meth_num = int(a * b * 15 * conv_density(459))
        sys_yaml = deepcopy(base_yaml)
        sys_yaml["xyz"] = name
        box_dims = [[a / 2, b / 2, 10.5], [-a / 2, -b / 2, -10.5]]
        sys_yaml["box_dims"] = box_dims
        sys_yaml["population"]["methane"] = meth_num
        # # print(sys_yaml)
        fptr = name.split(".")[0]
        fobj = open(str(systems_path / f"{fptr}.yaml"), "w")
        yaml.dump(sys_yaml, fobj)
        base_create_membrane(
            fptr, fptr + "_nvt", walltime="4:0:0", cores=8, debug=False
        )
        # raise Exception
