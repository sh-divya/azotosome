from squid import files, structures, geometry
from pathlib import Path
from copy import deepcopy
import numpy as np

BASE_PATH = Path(__file__).parent
struct_path = BASE_PATH / "structs"
param_path = BASE_PATH / "params"
config_path = BASE_PATH / "config"


def box_info(xyz_file):
    if isinstance(xyz_file, str):
        system = files.read_xyz(str(struct_path / xyz_file))
    else:
        system = deepcopy(xyz_file)

    max_x = max([a.x for a in system])
    min_x = min([a.x for a in system])

    max_y = max([a.y for a in system])
    min_y = min([a.y for a in system])

    max_z = max([a.z for a in system])
    min_z = min([a.z for a in system])

    # print("Lower bounds for XYZ")
    # print(min_x, min_y, min_z)

    # print("Upper bounds for XYZ")
    # print(max_x, max_y, max_z)

    # print("Box side lengths")
    a = round(max_x - min_x, 2)
    b = round(max_y - min_y, 2)
    c = round(max_z - min_z, 2)
    # print(f"a:{a}, b:{b}, c:{c}")

    return float(a), float(b), float(c)


def make_supercell(prefix, xyz, repeat, dist, offset, write=True):
    if not isinstance(xyz, list):
        xyz = [xyz]
        offset = [offset]

    supercell = []
    # rot_mat = geometry.rotation_matrix([0, 1, 0], 90)
    for f, fptr in enumerate(xyz):
        if isinstance(fptr, str):
            og = files.read_xyz(str(struct_path / fptr))
        else:
            og = deepcopy(fptr)
        # og = geometry.rotate_atoms(og, rot_mat)
        displace = offset[f]
        for atoms in og:
            atoms.translate(displace)
        xreps, yreps, zreps = repeat

        for x in range(xreps):
            for y in range(yreps):
                for z in range(zreps):
                    for atom in og:
                        temp = deepcopy(atom)
                        move_vect = [(2 * m + 1) * r for m, r in zip([x, y, z], dist)]
                        temp.translate(move_vect)
                        supercell += [temp]

    repeat = list(map(str, repeat))
    name = str(f"{prefix}_{''.join(repeat)}.xyz")

    if write:
        files.write_xyz(supercell, str(struct_path / name))
    return supercell, name


def place_acryls(
    num, stagger=False, distance=2.0, nx=None, prefix="custom_azotosome", write=True
):
    acryl = files.read_cml(str(struct_path / "scaled_acryl.cml"))[0]
    cog = geometry.get_center_of_mass(acryl.atoms)
    if not nx:
        length = int(num**0.5)
    else:
        length = nx
    monolayer = []
    row = 0
    rot_mat = geometry.rotation_matrix([0, 0, 1], 180)
    acryl.rotate(rot_mat)
    # rot_mat = geometry.rotation_matrix([0, 1, 0], 90)
    # acryl.rotate(rot_mat)
    if isinstance(distance, float):
        distance = [distance, distance, 0]

    for n in range(num):
        z = False
        new_mol = deepcopy(acryl)
        if n % length == 0 and n > 0:
            row += 1
        if row % 2:
            rot_mat = geometry.rotation_matrix([1, 0, 0], 180)
            new_mol.rotate(rot_mat)
            z = not z
        if n % 2:
            rot_mat = geometry.rotation_matrix([1, 0, 0], 180)
            new_mol.rotate(rot_mat)
            z = not z
        if z:
            z = distance[2]
        col = n % length
        new_mol.translate([col * distance[0], row * distance[1], z])
        monolayer.append(new_mol)

    frame = [atom for mol in monolayer for atom in mol.atoms]
    if write:
        files.write_xyz(frame, str(struct_path / f"{prefix}_{num}.xyz"))

    return frame


def conv_density(density=0.526, mol_mass=16):
    Na = 6.022e23
    kg_to_molnum = Na * (1000 / mol_mass)
    m3_to_A3 = 1e30
    return density * kg_to_molnum / m3_to_A3


def spacing(xyz_file):
    system = files.read_xyz(str(struct_path / xyz_file))
    xyz = []
    dist = []
    for a, atom in enumerate(system[::7]):
        single = [system[a + i] for i in range(7)]
        com = geometry.get_center_of_mass(single)
        xyz.append(com)

    x, y, z = [], [], []
    for pos1 in xyz:
        for pos2 in xyz:
            if (pos1 != pos2).all():
                tmp = pos1 - pos2
                x.append(np.abs(tmp[0]))
                y.append(np.abs(tmp[1]))
                z.append(np.abs(tmp[2]))
                tmp = np.linalg.norm(tmp)
                dist.append(tmp)

    print("Min x", min(x))
    print("Min y", min(y))
    print("Min z", min(z))


if __name__ == "__main__":
    # print(conv_density(459))
    a, b, c = box_info("tianle_azo_meth_421.xyz")
    # a, b, c = box_info("tianle_azo_221.xyz")
    # a, b, c = box_info("rahm_pn21a_small_222.xyz")
    # num_methane = int(a * b * 15 * conv_density(459))
    print(a, b, c)
    # print(num_methane)
    prefix = "tianle_azo_meth"
    num = "221"
    make_supercell(
        prefix,
        ["tianle_azo.xyz", "tianle_meth.xyz", "tianle_meth.xyz"],
        [2, 4, 1],
        [3.5, 3.5, 0],
        offset=[[0, 0, 0], [3, 0, 7], [3, 0, -5]],
    )

    prefix = "tianle_azo_vac"
    num = "221"
    make_supercell(
        prefix,
        ["tianle_azo.xyz"],
        [2, 2, 1],
        [3.5, 3.5, 0],
        offset=[[0, 0, 0]],
    )

    # prefix = "stag_azo"
    # num = "221"
    # make_supercell(
    #     prefix,
    #     ["stag_4.xyz"],
    #     [2, 2, 1],
    #     [4.0, 3.5, 0],
    #     offset=[[0, 0, 0]],
    # )

    # prefix = "rahm_pn21a_small"
    # make_supercell(
    #     prefix,
    #     ["pn21a_unit_boxrel.xyz"],
    #     [2, 2, 2],
    #     [3.64, 2.7, 4.22],
    #     offset=[[0, 0, 0]],
    # )

    # prefix = "stag"
    # num = 4
    # place_acryls(num, distance=[4.0, 3.5, -2.058], prefix=prefix)
    # place_acryls(num, distance=[3.5, 3.52, -2.058], prefix="tianle")
