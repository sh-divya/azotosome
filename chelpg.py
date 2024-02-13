from squid import files, orca
from pathlib import Path

BASE_PATH = Path(__file__).parent
struct_path = BASE_PATH / "structs"

acryl = files.read_cml(str(struct_path / "scaled_acryl.cml"))[0]
orca.job(
    "acrylonitrile",
    "!HF 6-31G* OPT CHELPG",
    atoms=acryl.atoms,
    queue="shared",
    nprocs=8,
    walltime="24:0:0",
    prebash="module load orca/4.1.2",
)
print("Done")
