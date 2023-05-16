# Contents 

[**1 - Introduction**](#1---introduction) <br /><br />
[**2 - Requirements**](#2---requirements) <br /><br />
[**3 - Usage**](#3---usage)  <br /><br />
[**4 - remove_H_atoms.py**](#4---remove_h_atomspy) <br />
&emsp; [4.1 - Scope](#41---scope) 
&emsp;&emsp;&emsp; [4.2 - Requirements](#42---requirements) 
&emsp;&emsp;&emsp; [4.3 - Usage](#43---usage) 
&emsp;&emsp;&emsp; [4.4 - Arguments](#44---arguments)
&emsp;&emsp;&emsp; [4.5 - Output](#45---output)  <br /><br />
[**5 - ResRel-MPI.py**](#5---resrel-mpipy) <br />
&emsp; [5.1 - Scope](#51---scope)
&emsp;&emsp;&emsp; [5.2 - Requirements](#52---requirements)
&emsp;&emsp;&emsp; [5.3 - Usage](#53---usage)
&emsp;&emsp;&emsp; [5.4 - Arguments](#54---arguments)
&emsp;&emsp;&emsp; [5.5 - Output](#55---output)  <br /><br />
[**6 - Hs-Hk-plot.py**](#6---hs-hk-plotpy)  <br />
&emsp; [6.1 - Scope](#61---scope)
&emsp;&emsp;&emsp; [6.2 - Tasks](#62---tasks) 
&emsp;&emsp;&emsp; [**6.3 - "Density" Task**](#63---density-task)
&emsp;&emsp;&ensp; [**6.4 - "Bin" Task**](#64---bin-task)
&emsp;&emsp;&emsp; [6.5 - Output](#65---output)<br />
&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&nbsp;&thinsp; [6.3.1 - Requirements](#631---requirements)
&emsp;&emsp;&nbsp;&nbsp;&thinsp;&thinsp; [6.4.1 - Requirements](#641---requirements)  <br />
&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&nbsp;&thinsp; [6.3.2 - Usage](#632---usage) 
&emsp;&emsp;&emsp;&emsp;&emsp;&emsp; [6.4.2 - Usage](#642---usage) <br />
&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&nbsp;&thinsp; [6.3.3 - Arguments](#633---arguments) 
&emsp;&emsp;&emsp;&emsp; [6.4.3 - Arguments](#643---arguments) <br />



# 1 - Introduction
When coarsening biomolecules, the identification of the optimal number of sites to minimize information loss from an all-atom conformation is a challenging task. Several coarse-grained and multi-resolution models have been developed to tackle this issue, and one promising model is CANVAS (Coarse-grained Anisotropic Network model for VAriable resolution Simulation).

The CANVAS strategy leverages the blurred and approximate nature of coarse-grained models to identify effective sites based on a user-provided input, and determines the interactions among them based on the molecule’s structure and all-atom force field, making it unnecessary to run reference simulations. This strategy makes the parametrisation of the model practically instantaneous, and allows the modulation of the system’s resolution in a quasi-continuous manner across the structure, from all-atom to (very) coarse-grained. Most notably, the interaction between regions of the system at different resolution (including the solvent) is accounted for and straightforward to set up, allowing the seamless implementation in standard MD software packages (e.g. GROMACS or LAMMPS).

In CANVAS model three levels of resolution are employed: `all-atom` where all the atoms of the system are token in account; `medium-grained` where the backbone atoms are retained and treated as CG beads; and finally `coarse-grained` where only the $C_\alpha$ atoms are kept and modelled as CG beads. 

However, this approach requires prior knowledge of the system's chemistry and biology to determine which parts necessitate a fully atomistic description, namely in which part of the system the chemical details have a significant impact. Answering this question can be challenging.

Recently, a new method called **Resolution and Relevance** has been developed to identify the optimal resolution level that balances simplicity and informativeness. This framework, also known as critical variable selection, allows for the identification of important variables without prior knowledge or assumptions about their nature. The core idea behind this approach is that the generative model underlying empirical samples can be inferred from the distribution of their frequencies, i.e., the number of times different outcomes occur in the dataset.

Building upon the aforementioned approach, our goal is to identify the optimal number of sites for multi-resolution protein descriptions using the CANVAS model. This project involves the combination of three methods/software: Relevance and Resolution, Mapping Entropy, and CANVAS: 

* The **`Relevance and Resolution software`**, written in Python, aims to determine the optimal number of sites for biomolecule coarse-graining.

* The **`Mapping Entropy tool`**, written in C, takes the optimal number of sites as input and returns the site selection that minimizes information loss during the reduction of degrees of freedom in a system.

* The **`CANVAS model`**, available on a GitHub repository and implemented in Python, allows for the modeling of biomolecules at three levels of resolution as described earlier. It requires the output from the Mapping Entropy tool as input.


In this section, we present the tool for identifying the optimal number of sites. Subsequently, this number will serve as input for the Mapping Entropy tool, which will provide the atom selection. Finally, an additional code will be necessary to find the CANVAS selection sites that are closest to the output of the Mapping Entropy tool. This automated process will facilitate the coarsening of proteins using the CANVAS model.

<br/>

# 2 - Requirements

* **`Python3`**: it is a powerful interpreted, object-oriented, and high-level programming language known for its dynamic semantics. It is highly recommended to use Python 3.7 or 3.9 as they are the most suitable versions. If you're working on a _Linux_ or _macOS_ system, Python 3 should already be installed. However, if you're using Windows, the presence of Python 3 is not guaranteed. To install Python 3, you can follow the installation guide provided [here](https://docs.python-guide.org/starting/installation/). Please ensure that you are working with Python 3 (preferably 3.7 or 3.9) as executing the code with Python 2 may result in errors or unexpected behavior.
  
* **`Python3 libraries`**: Python 3 comes with a wide range of built-in libraries that are installed by default. However, there are certain libraries that may need to be installed separately. Here are four libraries used in this code that typically require subsequent installation:

   * [**`MDAnalysis`**](https://www.mdanalysis.org): It is an open source Python library that helps to quickly write your own analysis algorithm 
                                                     for studying trajectories produced by the most popular simulation packages. 
         
   * [**`NumPy`**](https://numpy.org): It  stands for _Numerical Python_ and it a fundamental library for numerical computing in Python. 
                                       It provides support for large, multi-dimensional arrays and matrices, along with a collection of mathematical functions 
                                       to operate on these arrays efficiently. It also has functions for working in domain of linear algebra and fourier transform. 
                                       NumPy was created in 2005 by Travis Oliphant. It is an open source project and you can use it freely.            
              
   * [**`Matplotlib`**](https://matplotlib.org): It is a low level graph plotting library in python that serves as a visualization utility created by John D. Hunter. 
                                                 It is open source and we can use it freely. Moreover, Matplotlib is mostly written in python, 
                                                 a few segments are written in C, Objective-C and Javascript for Platform compatibility.
                                  
   * [**`SciPy`**](https://scipy.org): It is a free and open-source Python library used for scientific computing and technical computing. 
                                       It was created by Travis Oliphant. SciPy contains modules for optimization, linear algebra, integration, 
                                       interpolation, special functions, FFT, signal and image processing, ODE solvers and other tasks common 
                                       in science and engineering.

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
* **`remove_H_atoms.py`**: has the preliminary scope of removing all the hydrogen atoms  from both the reference file  and the trajectory one. The reason lies in the fact that, in the calculation of the RMDS map and in particular in the calculation of Resolution and Relevance after keeping a group of atoms, it is preferable to not consider the hydrogens, as they are not heavy-atoms. This code can be ignored if having  reference and trajectory without hydrogen atoms. Details can be found in **Sec. 4**. 

* **`ResRel-MPY.py`**: This is the core program beacuse has the scope of calculating the Relevance-Resolution (changing the number of sites and for different mappings) points after calculating the RSD map among each frame and the other ones (an allignment between a couple of frames is required every time). Basically, it writes a file splitted in 3 rows:
    * 1st row: values of Resolution (Hs); 
    * 2nd row: values of Relevance (Hk); 
    * 3rd row: number of retained sites for that specific Hs and Hk. 
Details can be found in **Sec. 5**. 
  
* **`Hs-Hk-plot`**: this code has a dual purpose:
    * Drawing a saving different plot regard Resolution and Relevance, slope, and histogram of frequencies 
    * Calculating the optimal number of sites of a biomolcule starting from an atomistic trajectory, such that the loss of information after decimating atoms is minimized.
Details will be provided in **Sec. 6**.

Before running the python scripts, read carefully the next section that provides a detailed explaination of each task and argument. Moreover, take care to not moving them outside the main folder (`ResRel-identification-Optimal-N-Sites/`) otherwise a fatal error occurs and it is printed on screen.

<br/>

# 4 - remove_H_atoms.py

## 4.1 - Scope 
This script has the preliminary, but important, scope of removing all the hydrogen atoms (H, H1, H2, HW, etc...) from both the reference file (usually 'gro' or 'pdf', 'psf', etc...) and the trajectory one ('xtc', 'trr', 'dcd', 'gro', 'lammpstrj', etc...). Hydrogen atoms in a protein move and rotate more compared to heavy atoms. This is because hydrogen atoms have a much smaller mass than heavy atoms like carbon, nitrogen, and oxygen. This phenomenon is known as the reduced mass effect and is an important factor in protein dynamics. The mobility of hydrogen atoms can influence the three-dimensional structure and stability of the protein. For this reason, in the calculation of the RSD map, on of the key points for computing Resolution and Relevance after keeping a group of atoms, it is recommended to avoid hydrogens. If you already have reference and trajectory without hydrogen, this code can be ignored.

## 4.2 - Requirements
This script requires two mandatory files: the coordinate/topology file of all-atom structure of the biomolecule (`gro`, `pdb`, `xyz`, `psf`, ...) and the trajectory file in any format (`lammpstrj`, `dcd`, `trr`, `xtc`, ...). No optional arguments are available. 


## 4.3 - Usage
In order to launch the **remove_H_atoms.py** scripts, the command-line is the following:

```sh
python3 remove_H_atoms.py -r <Coordinate FILE> -t <Trajectory FILE> 

   or:

python3 remove_H_atoms.py --ref <Coordinate FILE> --traj <Trajectory FILE>
```

A short explaination of arguments is provided by launching the command `python3 remove_H_atoms.py -h` or `python3 remove_H_atoms.py --help`. Alternatively, for printing a short usage message, please type: `python3 remove_H_atoms.py` or `python3 remove_H_atoms.py -u`

Before running the python scripts, read carefully the next section that provides a detailed explaination of each argument.


## 4.4 - Arguments
As shown in **Sec. 4** the coordinate/topology file of all-atom structure of the biomolecule (`gro`, `pdb`, `xyz`, `psf`, ...), and the trajectory file in any format (`lammpstrj`, `dcd`, `trr`, `xtc`, ...) are always, mandatory. Moreover, no optional arguments are available. A short explaination of the above mentioned files is the following:

* **`Coordinate FILE`**: Mandatory File of atom Coordinates (xyz, gro, pdb, psf, ...) 

* **`Trajectory FILE`**: Mandatory File containing the trajectory of the biomolecule (trr, dcd, lammpstrj, gro, ...)

Examples are reported in **Sec. XX**


## 4.5 - Output
The output of the program are the coordinate file (_`Reference_noH.gro`_) and the trajectory (_`Trajectory_noH.xtc`_) after removing all the hydrogen atoms. 


<br/>

# 5 - ResRel-MPI.py 

## 5.1 - Scope 
This is the core program because has the scope of calculating the Relevance-Resolution (changing the number of retained sites and for different mappings) points. Several steps are needed: 
1. Calculation of the all-atom _RSD map_ among each frame and the other ones (an allignment between a couple of frames is required every time); 
2. Creating the dendogram relative to the all-atom RSD map of trajectory using the average linkage UPGMA algorithm; 
3. Cutting the dendogram in order to identify the value of cutoff such that all atomistic conformation can be distinguished; 
4. Starting with a number of retained sites equals to _Natoms-1_, a random mapping is proposed. Then, based on the cutoff chosen before, it is possible to know the number of clusters at such cutoff in this configuation and, finally, Hs-Hk point is computed; 
5. The procedure in (4) is iterated for a number _M_ of mappings (by default _M = 50_); 
6. Reducing progressively the number _N_ of retained sites, (4) and (5) steps are done, until no atoms are retained.
7. At the end, a complete curve Hs-Hk points is drawn, ready to be analyzed, with the main purpose of calculating the optimal number of sites (next code).

For the sake of clarity, the 4-5-6 step described above are depicted in **Fig. 1** in terms of flux diagram, with the purpose of showing that two nested for-loop are required for the calcution of all Hs-Hk points: the external loop is done on the number of retained sites $N_s$, whereas the more internal one is on the _M_ random mappings at fixed number of retained sites.

<div align="center">

<img src="4-5-6.jpg" alt="Scheme" width="550">
</div>
<div align = "center">
<b>Fig.1</b> - <i> Schematic representation of Relevance and Resolution points calculation. Two nested for loop are required: the external one is done on the number of retained sites (N<sub>s</sub>), the more internal one is on the "M" random mappings at fixed number of retained sites.</i>
</div>


## 5.2 - Requirements

This script requires two mandatory files: the coordinate/topology file of all-atom structure of the biomolecule without hydrogen atoms (`gro`, `pdb`, `xyz`, `psf`, ...) and the trajectory file in any format (`lammpstrj`, `dcd`, `trr`, `xtc`, ...). On the other hand, four arguments are optional: 

* _`Nmappings`_: Number of random mappings at fixed number of sites retained; 
* _`Nframes`_: Number of frames read in our trajectory; 
* _`Nstep`_: Number that describes the step when the number of retained sites is changed;  
* _`ncpu`_: Number of cpus that will be used for parallizing the calculation of RSD map for each mapping. 

## 5.3 - Usage 

In order to launch the **ResRel-MPI.py** scripts, the command-line is the following:

```sh
python3 ResRel-MPI.py -r <Reference_noH.gro> -t <Trajectory_noH.xtc> [-m <NMappings>] [-f <Nframes>] [-s <Nsteps>] 

   or:

python3 remove_H_atoms.py --ref <Reference_noH.gro> --traj <Trajectory_noH.xtc> [--mapp NMappings>] [--frames <Nframes>] [--step <Nsteps>]
```

> **NOTE: Please, take in account that "Reference_noH.gro" and "Trajectory_noH.xtc" - i.e. the reference coordinate file and the trajectory one without hydrogen atoms, respectively - are the output files obtained after launching remove_H_atoms.py. This code does not return an error, if using files with hydrogens. However, in order that the calculation of Resolution and Relevance points for each mapping works good, it is necessary to remove H atoms beacuse the latter are not heavy atoms, and thus they move and rotate too much. Please, take care of it.**

A short explaination of arguments is provided by launching the command `python3 ResRel-MPI.py -h` or `python3 ResRel-MPI.py --help`. Alternatively, for printing a short usage message, please type: `python3 ResRel-MPI.py` or `python3 ResRel-MPI.py -u`.

Before running the python scripts, read carefully the next section that provides a detailed explaination of each argument.


## 5.4 - Arguments

As shown in **Sec. 5** the coordinate file and the trajectory of all-atom structure of the biomolecule without hydrogens (_Reference_noH.gro_ and _Trajectory_noH.xtc_, respectively) are, always, mandatory. On the other hand, the number of mappings at fixed number of sites (Nmappings), the number of frames read from trajectory (Nframes), and the step corresponing at the variation of the number of sites (Nstep) are optional arguments. A short explaination of the above mentioned files is the following:

* **`Coordinate FILE noH`**: Mandatory File (`-r/--ref`) of atom Coordinates **without** hydrogen atoms (xyz, gro, pdb, psf, ...). If using "remove_H_atoms,py" code, then the file is by default called _Reference_noH.gro_; 

* **`Trajectory FILE noH`**: Mandatory File (`-t/--traj`) containing the trajectory of the biomolecule **without** hydrogens (trr, dcd, lammpstrj, gro, ...). If using "remove_H_atoms,py" code, then the file is by default called _Trajectory_noH.gro_; 

* **`NMappings`**: Optional argument (`-m/--mapp`) that specifies the number of random mappings _M_ at fixed number of sites retained. Any integer number higher than 0 is accepted. The default value is _M = 50_. For a fixed number of sites the code choses randomly _M_ combinations with respect the total number of atoms. Changing such value, more or less combinations are chosen; 

* **`Nframes`**: Optional arguments (`-f/--frames`) that specifies the number of frames _F_ read in our trajectory. In order to guarantee this precise number of frames and the spanning of entire trajectory, an initial number of frames will be discarded and the trajectory is read every _nSteps_. The default value is _F = 1000_. Any integer number (less than the original number of frames) is accepted. If the string "_all_" is set, then every frame of trajectory is read (`-f all`). A trajectory, in general, could contain more than 1000 frames. However, since the calculation of the RSD map goes as Nframes squared ($N^2$), increasing the number of frames would require much more time. Higher the number of cores containing your cluster in a single node, more frames can be token in account. Be careful choosing such value or selecting all frames. 

* **`Nstep`**: Optional argument (`-s/--step`) that corresponds at a integer number which describes the variation of the number of retained sites in the calculation of 'Nmappings' Resolution and Relevance points. It means that, starting from the total number of atoms minus 1 (Natoms-1) 'Nmappings' Resolution and Relevance points are computed for each  mapping. Then, according with 'Nstep' value, the number of retained sites decreases of 'Nstep' (Natoms-1-Nstep) and other 'Nmappings' Resolution and Relevance points are computed for each mapping, and so on... until 3 retained sites are reached. The user have two ways of defining this argument:
    * As percentage of the total number of atoms: the number of retained sites in the for-loop employed in the calculation of Relevance and Resolution will decrease of such percentage. The default value is 0.5%. Therefore, it means that if the total number of atoms (Natoms) is 10'000, then the 0.5% of such number is 50. In each for-loop the number of retained sites starts from 10'000, and decreases by 50, until 3 is reached (Natoms, Natoms-50, Natoms-100, ..., 3). According with percentage definition, Nstep is an integer and float numbers between 0 and 100 followed by the '%' symbol i.e. <INT/FLOAT>% (no spaces between)
In case this percentage returns a step = 0, an error is printed on screen. Default value = 0.5%. 
    * It is also possible to write directly the step without calculating it as percentage of the total number of atoms. In such case, the step must be an INTEGER number between 1 and Natoms-1, otherwise an error is returned. For example if Natoms=10000 and Nstep=100, then starting from Natoms, in each for-loop the number of retained sites wll decrease by 100, until 3 is reached (Natoms, Natoms-100, Natoms-200, ..., 3); 
    
* **`NumberCpu`**: Optional argument (`-n/--ncpu`) that corresponds at the number of cpus that you would like to employ for parallizing the calculation of RSD map for each mapping. By default, the code will be parallelized by employing the maximum number of cores available in a single node in your laptop/cluster. However, such number can be easily changed using `-n/--ncpu <nCPU>`, and the code will be parallelized by employing nCPU cores. As just mentioned, on the other hand, in case `-n/--npu <nCPU>` is **not** set, the code will be parallelized by employing the maximum number of allowed cores in your laptop/cluster. 



## 5.5 - Output 

The output of this code is twofold: 

* **`trace_${ProteinName}.txt`**: it simply gives a trace of how many Resolution & Relevance points has already been calculated and the time required for calculating a single point at fixed number of retained sites (lower the retained sites, lower the time for calculating RSD-map and consequently for computing Hs-Hk point). Moreover, also an estimation of the remainder total time is done; 
    
* **`Hs-Hk-Nsites-${ProteinName}.txt`**: it is a more important file, as it contains all the value of resolution (Hs), relevance (Hk), and the number of sites retained associated to Hs and Hk. This file is organized in 3 rows. It will be useful for drawing plots an computing the optimal number of sites, that is the purpose of the third script, namely "Hs-Hk-plots.py". 
that Number that describes the variation of the number of retained sites
A short explaination of arguments is provided by launching the command `python3 ResRel-MPI.py -h` or `python3 ResRel-MPI.py --help`. Alternatively, for printing a short usage message, please type: `python3 ResRel-MPI.py` or `python3 ResRel-MPI.py -u`

# 6 - Hs-Hk-plot.py 

## 6.1 - Scope 
This code has dual purpose: 
1. It creates and saves four plots:
    * Resolution & Relevance plot: same colors are indicative of Hs-Hk points come out from same number of retained sites and different mappings.
        Moreover a zoom of the region of interest (where the slope is close to -1) is also shown on the same plot;      
    *  Zoom of Relevance & Resolution curve (Hk-Hs) in the windows where the slope is -1; 
    *  Slope against an increasing index (1 to N points); 
    *  Histogram of Frequencies, that is the number of sites with more occourrences.

2. Calculating the optimal number of sites of a biomolcule starting from an atomistic trajectory, such that the loss of information after decimating atoms is minimized; the result will be reported in "Opt-number-of-sites.txt". 

## 6.2 - Tasks 
The Relevance and Resolution plot is made up _N_ points (by default _N about 10000_ points). Thus, in order to find the optimal number of sites is necessary to "simplify" the ResRel curve. A TALE SCOPO average values of Resolution ($\bar H_s$) and Relevance ($\bar H_k$) are computed. This calculation can be done in two different ways:  

* **`density`**: The x-axes is divided in _X_ intervals such that each of one contains the same number of points _D_ (by default _D = 100_). Thus, the interval lenght is not fixed. What is fixed is the number (density) of Hs-Hk points in each interval. Then, in each interval the average values for Hs and Hk are computed (Hs_avg and Hk_avg). If using this option (recommended) the calculation of the average values for Resolution (Hs) and Relevance (Hk) is based on the same density of points. Indeed, every N points (the default value is 100) the average calculation for Hs and Hk is performed. In this way, the computation of average values is more fair since the lenght of interval is not fixed, but it is chosen in order that the number of points (Hs-Hk) is always the same. The default value for N is 100 points; however, such value con be changed using the flag [-d]. 
                                        
   <div align="center">
   <img src="density.jpg" alt="Scheme" width="800">
   </div>
   <div align = "center">
   <b>Fig.2</b> - <i> Pictorial representation.</i>
   </div>
                                        
                  

* **`bin`**: If using this option the x-axes is divided in _W_ windows (intervals) having the same lenght (by defaults _W = 50_). Thus, fixing the interval, the x-axis (Resolution) that goes from 0 to 1, by definition of Resolution, is divided in _W_ windows and the **bin** is thus defined as _1/W_. Then, in each interval the average values of Hs (Hs_avg) are the central values of bin lenght while Hk_avg is computed. Be careful using this option instead of "density" one, because its choice could be not ideal: indeed, the density of points along the curve is different: there will be windows very dense of points, and other ones with few points. In this way,  the computation of average values could be not fair and not precise because of different dennsity of points. Use this option with caution. The default value for _W_ is 50 windows; however, such value can be changed using the flag [-w]. 

   <div align="center">
   <img src="bin.jpg" alt="Scheme" width="800">
   </div>
   <div align = "center">
   <b>Fig.3</b> - <i> Pictorial representation.</i>
   </div>

According with one of the two options, follow Sec. 6.3 if the choice is "density", the Sec. 6.4 in case of "bin" choice. 

## 6.3 - "Density" Task 

### 6.3.1 - Requirements 
_`density`_ task requires one mandatory file, that is **`Hs-Hk-Nsites-${ProteinName}.txt`**, the core file that contains all the value of resolution (Hs), relevance (Hk), and the number of sites retained associated to Hs and Hk. On the other hand, two arguments are optional: 

* _`DensityPoints`_: integer number that specify the fixed number of points in each (variable lenght) interval.
* _`SlopeRange`_: it specifies how to find the best interval such that new average curve made up of Hs_avg and Hk_avg is close to -1.

All the details of the arguments just described can be found in **Sec.XXX**.

### 6.3.2 - Usage 
In order to launch the **density** task the command-line is the following:

```sh
python3 Hs-Hk-plot.py density -f <Hs-Hk-Nsites-${ProteinName}.txt> [-d <density>] [-s <range>] 

   or:
   
python3 Hs-Hk-plot.py density --file <Hs-Hk-Nsites-${ProteinName}.txt> [--DensityPoints <density>] [--SlopeRange <range>] 
```
> **NOTE: Please, take in account that "Hs-Hk-Nsites-${ProteinName}.txt" is the output of "ResRel-MPI.py" described in details in Sec.5. Such file contains all the value of resolution (Hs), relevance (Hk), and the number of sites retained associated to Hs and Hk**

For further information, please type on terminal `python3 Hs-Hk-plot.py density`


### 6.3.3 - Arguments 
As explained in Sec. 6.3.1, "density" task requires one mandatory file, and two optional files: 

* _`Hs-Hk-Nsites-${ProteinName}.txt`_: Mandatory file that contains the Resolution (Hs), Relevance (Hk) and the number of sites (N) in the 1st, 2nd and 3rd row, respectively. This is the ouptut of 'ResRel-MPI.py' program. It is organized in 3 rows. Each one contains the value of Hs, Hk, and Nsites separated by one space.
  ```
    ----------------------------------------------------
    | Hs-1       Hs-2       Hs-3       .....  Hs-N     |
    | Hk-1       Hk-2       Hk-3       .....  Hk-N     |
    | Nsites-1   Nsites-2   Nsites-3   .....  Nsites-N |
    ----------------------------------------------------
  ```
  
* _`DensityPoints`_: Optional argument (`-d/--DensityPoints`) that specifies the number of points _D_ that fall in each interval (of variable lenght) as shown in Fig.2. By default, such value is _D = 100_. In interval of same density, the calulation of Hs_avg and Hk_avg is computed. Using this flag the dafult value can be easily changed. 

* _`SlopeRange`_: Optional argument (`-s/--SlopeRange`) that specifies how to find the best interval of this average new curve for computing the optimal number of sites. From literature we know that the slope $\mu = -1$ is associated to the point of optimal tradeoff between parsimony of the representation (low resolution) and its informativeness (high relevance). Therefore, specifically, after computing the average values of Hs and Hk, whatever the option chosen, two ways of defining this argument are possible:
    * It is possible to define a range close to -1 in terms of percentage (default: from -1.10 to -0.90 that corresponds at 10% and taking the rightmost value in such range (i.e. with higher Resolution). Higher range values could be not good, since we want a slope very close to -1. Default value = 10%; 
    * Another possibility is to take the closest value of slope to -1. In this case 'closest' string has to be used [-s/--SlopeRange closest]

## 6.4 - "bin" Task 

### 6.4.1 - Requirements 
_`bin`_ task requires one mandatory file, that is **`Hs-Hk-Nsites-${ProteinName}.txt`**, the core file that contains all the value of resolution (Hs), relevance (Hk), and the number of sites retained associated to Hs and Hk. On the other hand, two arguments are optional: 

* _`NumberWindows`_: integer number that specifies the fixed lenght of interval in which the x-axis (Resolution) is divided, whereas the density of points in each interval is variable. 
* _`SlopeRange`_: it specifies how to find the best interval such that new average curve made up of Hs_avg and Hk_avg is close to -1.

All the details of the arguments just described can be found in **Sec.6.4.3**.

### 6.4.2 - Usage 
In order to launch the **bin** task the command-line is the following:

```sh
python3 Hs-Hk-plot.py bin -f <Hs-Hk-Nsites-${ProteinName}.txt> [-w <nWindows>] [-s <range>] 

   or:
   
python3 Hs-Hk-plot.py bin --file <Hs-Hk-Nsites-${ProteinName}.txt> [--NumeberWindows <nWindows>] [--SlopeRange <range>] 
```
> **NOTE: Please, take in account that "Hs-Hk-Nsites-${ProteinName}.txt" is the output of "ResRel-MPI.py" described in details in Sec.5. Such file contains all the value of resolution (Hs), relevance (Hk), and the number of sites retained associated to Hs and Hk**

For further information, please type on terminal `python3 Hs-Hk-plot.py bin`.

### 6.4.3 - Arguments 
As explained in Sec. 6.4.1, "density" task requires one mandatory file, and two optional files: 

* _`Hs-Hk-Nsites-${ProteinName}.txt`_: Mandatory file that contains the Resolution (Hs), Relevance (Hk) and the number of sites (N) in the 1st, 2nd and 3rd row, respectively. This is the ouptut of 'ResRel-MPI.py' program. It is organized in 3 rows. Each one contains the value of Hs, Hk, and Nsites separated by one space.
  ```
    ----------------------------------------------------
    | Hs-1       Hs-2       Hs-3       .....  Hs-N     |
    | Hk-1       Hk-2       Hk-3       .....  Hk-N     |
    | Nsites-1   Nsites-2   Nsites-3   .....  Nsites-N |
    ----------------------------------------------------
  ```
  
* _`NumberWindows`_: Optional argument (`-w/--NumberWindows`) that specifies the number of intervals into which the x-axes (Resolution) is divided. Thus, the lenght of the interval is fixed, whereas the number of points that fall inside each interval is variable, as shown schematically in Fig.3. By default, _W = 50_. Knowing the value of _W_ it is possible to define the bin as _1/W_ (at numeratore "1" is the total lenght of x-axis: indeed, the resolution goes between 0 and 1, and thus the lenght is actually 1). In interval of same bin lenght, the calulation of Hk_avg is computed, whereas by construction, the value of Hs_avg corresponds at the center value of binlenght. Using this flag (-w) the default value can be easily changed. 

* _`SlopeRange`_: Optional argument (`-s/--SlopeRange`) that specifies how to find the best interval of this average new curve for computing the optimal number of sites. From literature we know that the slope $\mu = -1$ is associated to the point of optimal tradeoff between parsimony of the representation (low resolution) and its informativeness (high relevance). Therefore, specifically, after computing the average values of Hs and Hk, whatever the option chosen, two ways of defining this argument are possible:
    * It is possible to define a range close to -1 in terms of percentage (default: from -1.10 to -0.90 that corresponds at 10% and taking the rightmost value in such range (i.e. with higher Resolution). Higher range values could be not good, since we want a slope very close to -1. Default value = 10%; 
    * Another possibility is to take the closest value of slope to -1. In this case 'closest' string has to be used [-s/--SlopeRange closest]


## 6.5 - Output
This code returns 4 plots in PDF format, and a TXT file: 
* _`Reso.pdf`_: Resolution & Relevance curve: same colors are indicative of Hs-Hk points come out from same number of retained sites and different mappings.
        Moreover a zoom of the region of interest (where the slope $\mu$ is close to -1) is also shown on the same plot;      
* _`Zoom-Reso.pdf`_: Zoom of Relevance & Resolution curve (Hk-Hs) in the windows where the slope is -1; 
* _`slope.pdf`_: Slope against an increasing index (_1_ to _N_ points); 
* _`histo_Nsites.pdf`_ Histogram of Frequencies, that is the number of sites with more occourrences.

* _`Opt-number-of-sites.txt`_: file showing a summary of the arguments employed and, more important, the optimal number of sites of a biomolcule starting from an atomistic trajectory, such that the loss of information after decimating atoms is minimized.



