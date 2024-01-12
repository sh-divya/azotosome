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

def make_supercell(xyz, repeat, dist):
    og = files.read_xyz(str(struct_path / xyz))
    supercell = []
    
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
    # raise Exception
    repeat= list(map(str, repeat))
    name = str(f"{xyz.split('.')[0]}_{''.join(repeat)}.xyz")
    return name, files.write_xyz(supercell, str(struct_path / name))

def place_acryls(num, stagger=False, distance=2.0):
    acryl = files.read_cml(str(struct_path / "acrylonitrile.cml"))[0]
    cog = geometry.get_center_of_mass(acryl.atoms)
    length = int(num ** 0.5)
    monolayer = []
    row = 0
    rot_mat = geometry.rotation_matrix([0, 1, 0], 90)
    rot_atoms = geometry.align_centroid(acryl.atoms)

    for atom1, atom2 in zip(acryl.atoms, rot_atoms):
        atom1.set_position(atom2.flatten())
    acryl.rotate(rot_mat)
    for n in range(num):
        new_mol = deepcopy(acryl)
        if n % 2 or row % 2:
            rot_mat = geometry.rotation_matrix([1, 0, 0], 180)
            new_mol.rotate(rot_mat)
        if n == length:
            new_mol.translate(0.0, row * distance, 0.0)
            row += 1
        new_mol.translate([n * distance, 0.0, 0.0])
        monolayer.append(new_mol)
        
    frame = [atom for mol in monolayer for atom in mol.atoms]
    files.write_xyz(frame, str(struct_path / f"custom_azotosome_{num}.xyz"))

if __name__ == "__main__":
    # box_info("rahm_unsolv_rehydrogen_reorg.xyz")
    # make_supercell("rahm_unsolv_rehydrogen_reorg.xyz", [2, 2, 1], [3.2, 4.0, 3.5])
    place_acryls(16)