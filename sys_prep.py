from squid import files, structures, geometry
from pathlib import Path
from copy import deepcopy

BASE_PATH = Path(__file__).parent
struct_path = BASE_PATH / "structs"
param_path = BASE_PATH / "params"
config_path = BASE_PATH / "config"


def box_info(xyz_file):
    system = files.read_xyz(str(struct_path / xyz_file))

    max_x = max([a.x for a in system])
    min_x = min([a.x for a in system])

    max_y = max([a.y for a in system])
    min_y = min([a.y for a in system])

    max_z = max([a.z for a in system])
    min_z = min([a.z for a in system])

    print("Lower bounds for XYZ")
    print(min_x, min_y, min_z)

    print("Upper bounds for XYZ")
    print(max_x, max_y, max_z)

    print("Box side lengths")
    print(max_x - min_x, max_y - min_y, max_z - min_z)


def make_supercell(prefix, xyz, repeat, dist, offset):
    if not isinstance(xyz, list):
        xyz = [xyz]
        offset = [offset]
    supercell = []
    rot_mat = geometry.rotation_matrix([0, 1, 0], 90)
    for f, fptr in enumerate(xyz):
        og = files.read_xyz(str(struct_path / fptr))
        og = geometry.rotate_atoms(og, rot_mat)
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

    print(len(supercell))
    repeat = list(map(str, repeat))
    name = str(f"{prefix}_{''.join(repeat)}.xyz")
    return name, files.write_xyz(supercell, str(struct_path / name))


def place_acryls(num, stagger=False, distance=2.0, nx=None, prefix="custom_azotosome"):
    acryl = files.read_cml(str(struct_path / "scaled_acryl.cml"))[0]
    cog = geometry.get_center_of_mass(acryl.atoms)
    if not nx:
        length = int(num**0.5)
    else:
        length = nx
    monolayer = []
    row = 0
    # rot_mat = geometry.rotation_matrix([1, 0, 0], 90)
    # acryl.rotate(rot_mat)
    # rot_mat = geometry.rotation_matrix([0, 0, 1], 90)
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
    files.write_xyz(frame, str(struct_path / f"{prefix}_{num}.xyz"))


def conv_density(density=0.526, mol_mass=16):
    Na = 6.022e23
    kg_to_molnum = Na * (1000 / mol_mass)
    m3_to_A3 = 1e30
    return density * kg_to_molnum / m3_to_A3


if __name__ == "__main__":
    # print(conv_density(438.9))
    # box_info("rahm_unsolv_rehydrogen_reorg.xyz")
    prefix = "azo_meth"
    num = "441"
    make_supercell(
        prefix,
        ["tianle_unit.xyz", "tianle_meth.xyz", "tianle_meth.xyz"],
        [4, 4, 1],
        [3.5, 3.5, 0],
        offset=[[0, 0, 0], [0, 0, 7], [0, 0, -7]],
    )
    # box_info("rahm_unsolv_rehydrogen_reorg_441.xyz")
    # num = 16
    # prefix = "tianle"
    # place_acryls(num, distance=[3.5, 3, -2.058], prefix="tianle")
    box_info(f"{prefix}_{num}.xyz")
