Python pipeline to submit MD simulations using LAMMPS for predefined systems of azotosome membranes. These are acrylonitrile prosoed to possibly support the evolution of life on Titan, the moon of Saturn. Titan has seas of methane at a temperature of 94K and 1.4 atom pressure, which is just above the freezing point of acrylonitrile, which has also been found on the surface of Titan.

## Requirements

The scripts in this repo were run using `python==3.9.15` You will need the following python libraries and software

- Python Libraries
    - clancylab-squid
    - numpy
    - click
    - omegaconf
    - hydra
    - pymatgen
    - https://github.com/sh-divya/supplements
    - pandas
- Software
    - CP2K (2023.1)
    - LAMMPS (23Jun22)
    - PACKMOL
    - Quantum Espresso (7.2)

## How to use

High level functioning of scripts:

- `sys_prep.py`
- `relax_md.py`
- `cp2k.py`
- `qe.py`
- `post.py`

### Directory Structure

You might find the following directories when using this repository

- `structs`: incldues all the XYZ or CML or LAMMPS dump files that can be used as input or that are output. Currently `Quantum Espresso` and `CP2K` output coordinates are not output to this directory
- `config`: YAML files that include parameters required to submit LAMMPS / QE / CP2K simulations; and files that have the contents of the input files for the respective simulations
- `cp2k`, `qe`, `lammps`, `orca`: submitted simulations are organized under these directories, with each simulation compiled under a single directory with the SLURM job name
- `params`: Pseudo-potential files for CP2K / QE and OPLS parameter files for LAMMPS
- `outputs`: a by-product of using hydra as part of the submission pipeline
- `sys-packmol`: the `PACKMOL` sytem generation takes place within this directory, with its input, outputs and log files

### Working Azotosome Monolayer

**Example for 221 supercell**

1. `sys_prep.py`
    - `make_supercell($file-name-here$, ["tianle_azo.xyz", "tianle_meth.xyz", "tianle_meth.xyz"], [2, 2, 1], [3.5, 3.5, 0], offset=[[0, 0, 0], [3, 0, 7], [3, 0, -5]])`
    - Get box size using `box_info` and prepare yaml files under `config/systems/$system_name$.yaml`
2. `python relax_md.py --system=$system_name$.yaml --name=$job-name-here$`

### Working Pn21a crystal

**Example for 221 supercell**

1. `sys_prep.py`
    - `make_supercell($file-name-here$, ["pn21a_unit_boxrel.xyz"], [2, 2, 2], [3.64, 2.7, 4.22], offset=[[0, 0, 0]])`
2. Repeat as in monolayer case

## Relevant Work