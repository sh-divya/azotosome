from squid import files
from pathlib import Path

BASE_PATH = Path(__file__).parent
struct_path = BASE_PATH / "structs"

def unwrap(mol, xyz, cell_dims=None, name=None):
    mol = files.read_cml(str(struct_path / mol))[0]
    sys = files.read_xyz(str(struct_path / xyz))
    mol_len = len(mol.atoms)
    for s, sub in enumerate(sys[::mol_len]):
        for bond in mol.bonds:
            a1, a2 = bond.atoms
            b1 = s * mol_len + a1.index - 1
            b2 = s * mol_len + a2.index - 1
            d1 = sys[b1].flatten()
            d2 = sys[b2].flatten()
            delta = (d1 - d2).tolist()
            for i, d in enumerate(delta):
                if abs(d) > cell_dims[i] / 2:
                    if d > 0:
                        d2[i] = d2[i] + cell_dims[i]
                    if d < 0:
                        d2[i] = d2[i] - cell_dims[i]
                    sys[b2].set_position(d2)

    if name:
        files.write_xyz(sys, str(struct_path / name))
    return sys

if __name__ == "__main__":
    unwrap("scaled_acryl.cml", "pn21a_rehydrogen.xyz", [7.06881, 6.37355, 6.90727], name="pn21a_unwrap.xyz")