# 0 - Introduction 

When coarsening biomolecules, the identification of the optimal number of sites such that the loss of information starting from an all-atom conformation of biomoleculs is not trivial. So far, a lot of coarse-grained and multi-resolution model have been developed. One of the most promising one is the CANVAS model. CANVAS is the acronym of Coarse-grained Anisotropic Network model for VAriable resolution Simulation.

The CANVAS strategy leverages the blurred and approximate nature of coarse-grained models to identify effective sites based on a user-provided input, and determines the interactions among them based on the molecule’s structure and all-atom force field, making it unnecessary to run reference simulations. This strategy makes the parametrisation of the model practically instantaneous, and allows the modulation of the system’s resolution in a quasi-continuous manner across the structure, from all-atom to (very) coarse-grained. Most notably, the interaction between regions of the system at different resolution (including the solvent) is accounted for and straightforward to set up, allowing the seamless implementation in standard MD software packages (e.g. GROMACS or LAMMPS).

In CANVAS model three levels of resolution are employed: `all-atom` where all the atoms of the system are token in account; `medium-grained` where all the backbone atoms are retained and treated as CG beads; and finally `coarse-grained` where only the C-alpha atoms are kept and modelled as CG beads. 

However, in this approach, it is required the knowledge of system in order to understand which part of it requires a fully atomistic description, namely we must know which part of the system plays a crucial role, in which part of the system the chemical details have a major effect. Asnwering this question could not be easy. 

Recently, a new approach has been developed, named Resolution and Relevance, with the purpose of identifying the level of resolution that optimally balances simplicity and informativeness. The resolution-relevance framework, or critical variable selection, is a recently developed method for identifying important variables without any prior knowledge of, or assumption on, their nature. The idea at the heart of the approach is that the information on the generative model that underlies the elements of an empirical sample is contained in the distribution of their frequencies, that is to say, in the number of times different outcomes occur in the data set. 

Based on this approach, our purpose is to identify the optimal number of sites when using a coarse-grained and, in particular, a multiple resolution description of proteins by using CANVAS model. This project requires a combination of three methods/softwares: Relevance and Resolution, Mapping Entropy, and CANVAS:
* The former software is entirely written in Python and it has the scope of identifying the optimal number of sites when coarse graining biomolecules.
* The Mapping Entropy tool, on the other hand, is written in C. In input requires the optimal number of sites, and returns the selection of the sites the minimizes the loss of information when reducing the number of degrees of a system.
* The CANVAS model whose code is written in Python and available on github repository, allows to model a biomolecule with three levels of resolutions as above explained. In input requires the result come out from ME tool.

In this section the tool for identifying the optimal number of sites is reported. The final output will be, indeed, the Number of Optimal sites. 
Afterwards this number will be the output of the ME tool that gives actually in output the selection of atoms. Finally, another code will be required whose scope will be to find the CANVAS selection sites closest to the output of ME. 

In this way the process of coarsening protein with a CANVAS model will be completely automatized. 

<br/>

# 1 - Usage 

The typical usage of the program consists in a call to `remove_H_atoms.py`, `ResRel-MPI-py` and `Hs-Hk-plot.py` in succession by using Python3: 
* **`remove_H_atoms.py`**: has the preliminary scope of removing all the hydrogen atoms (H, H1, H2, HW, etc...) from both the reference file (usually 'gro' or 'pdf', 'psf', etc...) and the trajectory one ('xtc', 'trr', 'dcd', 'gro', 'lammpstrj', etc...). The reason lies in the fact that, in the calculation of the RMDS map and in particular in the calculation of Resolution and Relevance after keeping a group of atoms, we do not want to consider the hydrogens, as they are not heavy-atoms. If you already have reference and trajectory without hydrogen, this code can be ignored. 

* **`ResRel-MPY.py`**: this is the core program beacuse has the scope of calculating the Relevance-Resolution (changing the number of sites and for different mappings) points after calculating the RSD map among each frame and the other ones (an allignment between a couple of frames is required every time). Basically, it writes a file splitted in 3 rows:
    * 1st row: values of Resolution (Hs)
    * 2nd row: values of Relevance (Hk)
    * 3rd row: number of retained sites for that specific Hs and Hk 
  
* **`Hs-Hk-plot`**: this code has a dual purpose:
    * Drawing a saving different plot regard Resolution and Relevance, slope, and histogram of frequencies 
    * Calculating the optimal number of sites of a biomolcule starting from an atomistic trajectory, such that the loss of information after decimating atoms is minimized.

Details will be provided in Sec. XXX.

Before running the python scripts, read carefully the next section that provides a detailed explaination of each task and argument. Moreover, take care to not moving them outside the main folder (`ResRel-identification-Optimal-N-Sites/`) otherwise a fatal error occurs and it is printed on screen.

<br/>

# 2 - remove_H_atoms.py

This script requires two mandatory files: the coordinate/topology file of all-atom structure of the biomolecule (`gro`, `pdb`, `xyz`, `psf`, ...) and the trajectory file in any format (`lammpstrj`, `dcd`, `trr`, `xtc`, ...). No optional arguments are available. 

In order to launch the **remove_H_atoms.py** scripts, the command-line is the following:

```sh
python3 remove_H_atoms.py -r <Coordinate FILE> -t <Trajectory FILE> 

   or:

python3 remove_H_atoms.py --ref <Coordinate FILE> --traj <Trajectory FILE>
```

The output of the program are the coordinate file (_`Reference_noH.gro`_) and the trajectory (_`Trajectory_noH.xtc`_) after removing all the hydrogen atoms. For further information, please type on terminal `python3 remove_H_atoms.py` or `python3 remove_H_atoms.py -h`. 

Before running the python scripts, read carefully the next section that provides a detailed explaination of each argument.

## 2.1 - Arguments of "remove_H_atoms.py"

As shown in **Sec. 2** the coordinate/topology file of all-atom structure of the biomolecule (`gro`, `pdb`, `xyz`, `psf`, ...), and the trajectory file in any format (`lammpstrj`, `dcd`, `trr`, `xtc`, ...) are always, mandatory. Moreover, no optional arguments are available. A short explaination of the above mentioned files is the following:

* **`Coordinate FILE`**: Mandatory File of atom Coordinates (xyz, gro, pdb, psf, ...) 

* **`Trajectory FILE`**: Mandatory File containing the trajectory of the biomolecule (trr, dcd, lammpstrj, gro, ...)

Examples are reported in **Sec. XX**

<br/>

# 3 - ResRel-MPI.py 
