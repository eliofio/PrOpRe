# 1 - Introduction 

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

# 2- Requirements

* **`Python3`**: it is an interpreted, object-oriented, high-level programming language with dynamic semantics. 
  The installation guide is provided [Here](https://docs.python-guide.org/starting/installation/). 
  If you are working on _Linux_ or _MacOs_ system, Python3 should be already installed. 
  On the other hand, if you are using Windows operating system, it is not certain for its presence.
  Please, be care of working with Python3 (3.7 or 3.9 is the best choice) as the code could return an error if using Python2.
  
* **`Python3 Libraries`**: The most of libraries used in these codes are installed by default by python after first installation. However, four of them usually require a subsequent installation: 

   * [**`MDAnalysis`**](https://www.mdanalysis.org): It is an open source Python library that helps you to quickly write your own analysis algorithm 
                       for studying trajectories produced by the most popular simulation packages. 
         
   * [**`NumPy`**](https://numpy.org): It  stands for _Numerical Python_ and it is a Python library used for working with arrays. 
                  It also has functions for working in domain of linear algebra, fourier transform, and matrices. 
                  NumPy was created in 2005 by Travis Oliphant. It is an open source project and you can use it freely.   
              
   * [**`Matplotlib`**](https://matplotlib.org): It is a low level graph plotting library in python that serves as a visualization utility created by John D. Hunter. 
                       It is open source and we can use it freely. Moreover, Matplotlib is mostly written in python, a few segments are written in C,
                       Objective-C and Javascript for Platform compatibility.
                   
   * [**`SciPy`**](https://scipy.org): It is a free and open-source Python library used for scientific computing and technical computing. 
                  It was created by Travis Oliphant. SciPy contains modules for optimization, linear algebra, integration, interpolation, 
                  special functions, FFT, signal and image processing, ODE solvers and other tasks common in science and engineering.


              
                   
      To install the lastest stable releases with conda do:
      
      ```bash 
      conda config --add channels conda-forge
   
      conda install mdanalysis
      conda install numpy
      conda install matplotlib
      conda install scipy
      ```
   
      On the other hand, to install the latest stable release with pip or pip3 (which should be available in all Python installations) do:

      ```bash
      pip3 install --upgrade MDAnalysis
      pip3 install numpy
      pip3 install matplotlib
      pip3 install scipy
      ```
   
<br/>

# 3 - Usage 

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

# 4 - remove_H_atoms.py

This script requires two mandatory files: the coordinate/topology file of all-atom structure of the biomolecule (`gro`, `pdb`, `xyz`, `psf`, ...) and the trajectory file in any format (`lammpstrj`, `dcd`, `trr`, `xtc`, ...). No optional arguments are available. 

In order to launch the **remove_H_atoms.py** scripts, the command-line is the following:

```sh
python3 remove_H_atoms.py -r <Coordinate FILE> -t <Trajectory FILE> 

   or:

python3 remove_H_atoms.py --ref <Coordinate FILE> --traj <Trajectory FILE>
```

The output of the program are the coordinate file (_`Reference_noH.gro`_) and the trajectory (_`Trajectory_noH.xtc`_) after removing all the hydrogen atoms. A short explaination of arguments is provided by launching the command `python3 remove_H_atoms.py -h` or `python3 remove_H_atoms.py --help`. Alternatively, for printing a short usage message, please type: `python3 remove_H_atoms.py` or `python3 remove_H_atoms.py -u`

Before running the python scripts, read carefully the next section that provides a detailed explaination of each argument.

## 4.1 - Arguments of "remove_H_atoms.py"

As shown in **Sec. 4** the coordinate/topology file of all-atom structure of the biomolecule (`gro`, `pdb`, `xyz`, `psf`, ...), and the trajectory file in any format (`lammpstrj`, `dcd`, `trr`, `xtc`, ...) are always, mandatory. Moreover, no optional arguments are available. A short explaination of the above mentioned files is the following:

* **`Coordinate FILE`**: Mandatory File of atom Coordinates (xyz, gro, pdb, psf, ...) 

* **`Trajectory FILE`**: Mandatory File containing the trajectory of the biomolecule (trr, dcd, lammpstrj, gro, ...)

Examples are reported in **Sec. XX**

<br/>

# 5 - ResRel-MPI.py 

This script requires two mandatory files: the coordinate/topology file of all-atom structure of the biomolecule without hydrogen atoms (`gro`, `pdb`, `xyz`, `psf`, ...) and the trajectory file in any format (`lammpstrj`, `dcd`, `trr`, `xtc`, ...). On the other hand, three arguments are optional: 

* _`Nmappings`_: Number of random mappings at fixed number of sites retained; 
* _`Nframes`_: Number of frames read in our trajectory; 
* _`Nstep`_: Number that describes the step when the number of retained sites is changed;  

In order to launch the **ResRel-MPI.py** scripts, the command-line is the following:

```sh
python3 ResRel-MPI.py -r <Reference_noH.gro> -t <Trajectory_noH.xtc> [-m <NMappings>] [-f <Nframes>] [-s <Nsteps>] 

   or:

python3 remove_H_atoms.py --ref <Reference_noH.gro> --traj <Trajectory_noH.xtc> [--mapp NMappings>] [--frames <Nframes>] [--step <Nsteps>]
```
> **NOTE: Please, take in account that "Reference_noH.gro" and "Trajectory_noH.xtc" - i.e. the reference coordinate file and the trajectory one without hydrogen atoms, respectively - are the output files obtained after launching remove_H_atoms.py. This code does not return an error, if using files with hydrogens. However, in order that the calculation of Resolution and Relevance points for each mapping works good, it is necessary to remove H atoms beacuse the latter are not heavy atoms, and thus they move and rotate too much. Please, take care of it.**

