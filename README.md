# Contents 

**[1 - Introduction](#1---introduction)** <br /><br />
**[2 - Requirements](#2---requirements)** <br /><br />
**[3 - Usage](#3---usage)**  <br /><br />
**[4 - remove_H_atoms.py](#4---remove_h_atomspy)** <br />
&emsp; [4.1 - Scope](#41---scope) 
&emsp;&emsp;&emsp; [4.2 - Requirements](#42---requirements) 
&emsp;&emsp;&emsp; [4.3 - Usage](#43---usage) 
&emsp;&emsp;&emsp; [4.4 - Arguments](#44---arguments)
&emsp;&emsp;&emsp; [4.5 - Output](#45---output)  <br /><br />
**[5 - ResRel-MPI.py](#5---resrel-mpipy)** <br />
&emsp; [5.1 - Scope](#51---scope)
&emsp;&emsp;&emsp; [5.2 - Requirements](#52---requirements)
&emsp;&emsp;&emsp; [5.3 - Usage](#53---usage)
&emsp;&emsp;&emsp; [5.4 - Arguments](#54---arguments)
&emsp;&emsp;&emsp; [5.5 - Output](#55---output)  <br /><br />
**[6 - Hs-Hk-plot.py](#6---hs-hk-plotpy)**  <br />
&emsp; [6.1 - Scope](#61---scope)
&emsp;&emsp;&emsp; [6.2 - Tasks](#62---tasks) 
&emsp;&emsp;&emsp; **[6.3 - "Density" Task](#63---density-task)**
&emsp;&emsp;&ensp; **[6.4 - "Bin" Task](#64---bin-task)**
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

The typical usage of the program consists in a call to `remove_H_atoms.py`, `ResRel-MPI.py` and `Hs-Hk-plot.py` in succession by using Python3: 

* **`remove_H_atoms.py`**: It has the preliminary purpose of removing all hydrogen atoms from both the reference file and the trajectory file. The reason for this step is that, in the calculation of the Root Square Deviation (RSD) map, which is a key component for computing the Resolution and Relevance, it is preferable to exclude hydrogen atoms as they are not heavy atoms. It is important to note that if the reference and trajectory files already exclude hydrogen atoms, you can ignore this code. Additional details can be found in **[Section 4](#4---remove_h_atomspy)**.

* **`ResRel-MPI.py`**: This is the core program beacuse has the scope of calculating the Relevance and Resolution points (changing the number of sites and exploring different mappings) by analyzing the RSD map among each frame and the other ones. The program generates an output file with three rows of data:
    * 1<sup>st</sup> row: values of Resolution ($H_s$); 
    * 2<sup>nd</sup> row: values of Relevance ($H_k$); 
    * 3<sup>rd</sup> row: number of retained sites for that specific $H_s$ and $H_k$ values. 
For a more comprehensive understanding of the program's functionality and implementation details, please refer to **[Section 5](#5---resrel-mpipy)** in this documentation.
 
* **`Hs-Hk-plot`**: this code serves a dual purpose:
    * Drawing a saving various plots related to Resolution and Relevance, slope, and histogram of frequencies. 
    * Calculating the optimal number of sites for a biomolecule based on an atomistic trajectory, with the aim of minimizing the loss of information when atoms are decimated.
Further details regarding the functionality and usage of this code will be provided in **[Section 6](#6---hs-hk-plotpy)**.

Before running the Python scripts, it is important to read the next section carefully, as it provides a detailed explanation of each task and argument. Additionally, it is crucial to avoid moving the scripts outside the main folder (`ResRel-identification-Optimal-N-Sites/`) otherwise a fatal error occurs which will be displayed on the screen.

<br/>

# 4 - remove_H_atoms.py

## 4.1 - Scope 
This script serves a preliminary yet significant purpose of removing all hydrogen atoms (`H`, `H1`, `H2`, `HW`,...) from both the reference file (typically `gro`, `pdb`, `xyz`, `psf`, etc.) and the trajectory file (`xtc`, `trr`, `dcd`, `gro`, `lammpstrj`,...). Hydrogen atoms in proteins exhibit greater movement and rotation compared to heavy atoms such as carbon, nitrogen, and oxygen. This behavior is attributed to the significantly smaller mass of hydrogen atoms. Known as the reduced mass effect, this phenomenon plays a crucial role in protein dynamics. The mobility of hydrogen atoms can influence the three-dimensional structure and stability of the protein. Therefore, in the computation of the RSD map and the subsequent determination of Resolution and Relevance, it is advisable to exclude hydrogen atoms. If your reference and trajectory files already lack hydrogen atoms, you can ignore this code.

## 4.2 - Requirements
This script requires two mandatory files: the coordinate/topology file epresenting the all-atom structure of the biomolecule (`gro`, `pdb`, `xyz`, `psf`, ...) and the trajectory file in any format (`lammpstrj`, `dcd`, `trr`, `xtc`, ...). No optional arguments are available. 


## 4.3 - Usage
In order to launch the **remove_H_atoms.py** scripts, the command-line is the following:

```sh
python3 remove_H_atoms.py -r <Coordinate FILE> -t <Trajectory FILE> 

   or:

python3 remove_H_atoms.py --ref <Coordinate FILE> --traj <Trajectory FILE>
```

To obtain a brief explanation of the arguments, you can execute the command `python3 remove_H_atoms.py -h` or `python3 remove_H_atoms.py --help`. Additionally, if you wish to print a concise usage message, you can use either `python3 remove_H_atoms.py` or `python3 remove_H_atoms.py -u`.

However, it is strongly recommended to read the following section attentively before running the Python scripts as it provides a comprehensive explanation of each argument.


## 4.4 - Arguments
As described in **Section 4**, both the coordinate/topology file of the all-atom structure of the biomolecule and the trajectory file in any format  are always mandatory. There are no optional arguments available. Here is a brief explanation of the aforementioned files:


* **`Coordinate FILE`**: This file is mandatory and it that contains the atom coordinates. It can be in formats such as _xyz_, _gro_, _pdb_, _psf_, etc. It provides information on the positions of the atoms and their connectivity.

* **`Trajectory FILE`**: This file is mandatory and it contains the trajectory information of the biomolecule. It can be in formats such as _trr_, _dcd_, _lammpstrj_, _gro_, etc.  It provides information on how the biomolecule moves over time.


## 4.5 - Output
The program produces two output files:

1. _`Reference_noH.gro`_: This file contains the coordinate data of the biomolecule's all-atom structure after eliminating all the hydrogen atoms. It represents the updated coordinates of the remaining heavy atoms.
2. _`Trajectory_noH.xtc`_: This file comprises the trajectory data of the biomolecule after the removal of hydrogen atoms. It captures the time-dependent movement of the remaining heavy atoms over the course of the simulation.





<br/>

# 5 - ResRel-MPI.py 

## 5.1 - Scope 
This program serves as the core component for calculating Resolution-Relevance (H<sub>s</sub>-H<sub>k</sub>) points, which involves varying the number of retained sites and different mappings. The calculation process involves several steps:

1. Calculation of the all-atom _RSD map_ between each frame and the other frames, requiring alignment between each frame pair.
2. Construction of a dendrogram based on the all-atom RSD map of the trajectory using the average linkage UPGMA algorithm.
3. Cutting the dendrogram to determine the cutoff value that allows for the distinction of all atomistic conformations.
4. Starting with a number of retained sites equal to $N_{\text{atoms}}-1$, a random mapping is proposed. Based on the pre-determined cutoff, the number of clusters at that cutoff is enstablished for this configuration, and the H<sub>s</sub>-H<sub>k</sub> point is computed.
5. Steps 4 is repeated for a specified number _M_ of mappings (by default, _M = 50_).
6. Gradually reducing the number _N_ of retained sites, steps 4 and 5 are iterated until no atoms are retained.
7. Finally, a complete curve of H<sub>s</sub>-H<sub>k</sub> points is drawn, ready for analysis, with the primary goal of calculating the optimal number of sites (look **[Section 6](#6---hs-hk-plotpy)**).

For better comprehension, **Figure 1** illustrates the flux diagram representing steps 4-5-6, demonstrating that two nested for-loops are required to calculate all the Hs-Hk points: the outer loop iterates over the number of retained sites ($N_s$), while the inner loop performs $M$ random mappings at a fixed number of retained sites.

<div align="center">

<img src="4-5-6.jpg" alt="Scheme" width="550">
</div>
<div align = "center">
<b>Figure 1</b> - <i> Schematic representation of the calculation process for Relevance and Resolution points. This process involves two nested loops: the outer loop iterates over the number of retained sites (N<sub>s</sub>), while the inner loop iterates over 'M' random mappings at a fixed number of retained sites.</i>
</div>


## 5.2 - Requirements
To run this script, two mandatory files are required: the coordinate/topology file of the biomolecule without hydrogen atoms `gro`, `pdb`, `xyz`, `psf`, ...) and the trajectory file in any format (`lammpstrj`, `dcd`, `trr`, `xtc`, ...). Additionally, four optional arguments can be specified:

* _`Nmappings`_: number of random mappings generated at each fixed number of retained sites.
* _`Nframes`_: number of frames to be read in the trajectory file.
* _`Nstep`_: step that describes the decrement in the number of sites to be retained, starting from $N_{\text{atoms}}-1$, during the calculation. 
* _`ncpu`_: number of CPUs to be used for parallelizing the calculation of the RSD map for each mapping.


## 5.3 - Usage 

To run the **ResRel-MPI.py** script, the command-line is the following:

```sh
python3 ResRel-MPI.py -r <Reference_noH.gro> -t <Trajectory_noH.xtc> [-m <NMappings>] [-f <Nframes>] [-s <Nsteps>] 

   or:

python3 remove_H_atoms.py --ref <Reference_noH.gro> --traj <Trajectory_noH.xtc> [--mapp NMappings>] [--frames <Nframes>] [--step <Nsteps>]
```
> **NOTE: Please note that the "Reference_noH.gro" and "Trajectory_noH.xtc" files mentioned here refer to the output files obtained after running the `remove_H_atoms.py` script. It is crucial to remove hydrogen atoms from the files for accurate calculation of Resolution and Relevance points. Although this code does not throw an error if hydrogen atoms are present, their excessive movement and rotation can adversely affect the calculation. Please ensure that you have removed hydrogen atoms from the files before proceeding with the calculation.

To obtain a brief explanation of the arguments, you can execute the command `python3 ResRel-MPI.py -h` or `python3 ResRel-MPI.py --help`. Additionally, if you wish to print a concise usage message, you can use either `python3 ResRel-MPI.py` or `python3 ResRel-MPI.py -u`.

However, it is strongly recommended to read the following section attentively before running the Python scripts as it provides a comprehensive explanation of each argument.



## 5.4 - Arguments
In **Section 5**, it is emphasized that the coordinate file (_Reference_noH.gro_) and the trajectory (_Trajectory_noH.xtc_) of the all-atom structure of the biomolecule without hydrogen atoms are always required inputs. These files provide the necessary information for the calculation.  On the other hand, the number of mappings at fixed number of sites (_Nmappings_), the number of frames to be read from trajectory (_Nframes_), and the step that describes the decrement in the number of sites to be retained (_Nstep_) are optional arguments. The following is the summary of the different files and input parameters required by the code:

* **`Coordinate FILE noH`**: This is a mandatory file (`-r/--ref`) containing the atom coordinates of the biomolecule _without_ hydrogen atoms (in formats such as xyz, gro, pdb, psf, etc.). If the "remove_H_atoms.py" script is used, the default name for this file is _Reference_noH.gro_.

* **`Trajectory FILE noH`**: This is another mandatory file (`-t/--traj`) containing the trajectory of the biomolecule _without_ hydrogen atoms (in formats such as trr, dcd, lammpstrj, gro, etc.). If the "remove_H_atoms.py" script is used, the default name for this file is _Trajectory_noH.gro_.

* **`NMappings`**: This is an optional argument (`-m/--mapp`) that specifies the number of random mappings $M$, namely the number of combinations that will be chosen randomly, at a fixed number of retained sites. Each mapping represents a unique combination of atoms. By default, the value of _Nmappings_ $M$ is set to 50, meaning that 50 random mappings will be chosen. However, you have the flexibility to adjust this value according to your specific needs. Increasing the value of Nmappings will result in a greater number of random mappings being generated, while decreasing it will yield fewer mappings. Choosing a higher value for _Nmappings_ can provide a more comprehensive exploration of different atom combinations, but it will also increase the computational time required for the calculation. Conversely, selecting a lower value will reduce the computational burden but may result in a less exhaustive sampling of mappings. By adjusting the Nmappings argument, you can strike a balance between computational efficiency and the level of exploration of different mappings that suits your specific requirements. The default value of Nmappings, which is set to 50, ($M = 50$) serves as a good compromise between the two factors. You can adjust the Nmappings value according to your specific needs and the available computational resources.


* **`Nframes`**: This is an optional parameter (`-f/--frames`) that allows you to specify the number of frames $F$ to be read from the trajectory. The program ensures that this exact number of frames is included in the analysis, spanning the entire trajectory. To achieve this, an initial number of frames will be discarded, and the trajectory will be read at regular intervals. The default value for _Nframes_ is set to 1000, meaning that 1000 frames will be considered for the analysis. However, you have the flexibility to adjust this value according to your specific requirements. Any integer number less than the original number of frames in the trajectory is accepted. If you set the _Nframes_ argument to the string "_all_", indicated by `-f all`, the program will read every frame available in the trajectory. While this option allows for a comprehensive analysis of the entire trajectory, it's important to note that the calculation of the RSD (Root Squared Deviation) map involves a computational complexity proportional to the square of the number of frames (Nframes<sup>2</sup>). Consequently, increasing the number of frames will significantly increase the computation time. The ability to read more frames depends on the number of cores available in a single node of your computing environment. If you have a higher number of cores, you can process more frames efficiently. However, it's crucial to exercise caution when choosing the value of _Nframes_ or selecting the "all" option, as the computational resources required can grow substantially. Consider your specific analysis needs, the computational resources at your disposal, and the desired trade-off between computational time and analysis comprehensiveness when selecting the appropriate value for _Nframes_.



* **`Nstep`**: This is an optional argument (`-s/--step`) that describes the decrement in the number of sites to be retained, starting from $N_{\text{atoms}}-1$ during the calculation of Resolution and Relevance points for each mapping. In other words, this parameter allows for controlling the granularity of the reduction in the number of retained sites during the calculation, considering the specific number of mappings conducted. Here are the two ways you can define this argument:

    * Percentage of the total number of atoms: By specifying Nstep as a percentage, the number of retained sites in the for-loop used for calculating Relevance and Resolution will decrease by that percentage. The default value is 0.5%. For example, if the total number of atoms (Natoms) is 10,000, then 0.5% of this number is 50. In each iteration of the for-loop, the number of retained sites will start from 9999 ($N_{\text{atoms}}-1$) and decrease by 50 until reaching 3 ($N_{\text{atoms}}-1$, $N_{\text{atoms}}-51$, $N_{\text{atoms}}-101$, ..., 3). When using the percentage definition, _Nstep_ should be specified as an integer or float between 0 and 100, followed by the '%' symbol, without any spaces (e.g. _0.5%_ and **not** _0.5 %_). If the provided percentage results in _Nstep = 0_, an error will be displayed. 
    * Directly specifying the step: Alternatively, you can directly define the step without calculating it as a percentage of the total number of atoms. In this case, _Nstep_ should be an integer between $1$ and $N_{\text{atoms}}-1$. If a value outside this range is provided, an error will be raised. For example, if _Natoms_ is 10000 and _Nstep_ is set to $100$, the number of retained sites will decrease by $100$ in each iteration of the for-loop until reaching 3 ($N_{\text{atoms}}-1$, $N_{\text{atoms}}-101$, $N_{\text{atoms}}-201$, ..., 3).

  These options allow you to customize the decrement of the number of retained sites during the calculation process, providing flexibility in exploring different resolutions and optimizing computational efficiency.
  
* **`NumberCpu`**: This is an optional parameter (-n/--ncpu) that determines the number of CPUs used for parallelizing the calculation of the RSD (Root Square Deviation) map for each mapping. By default, the code will automatically utilize the maximum number of available cores in a single node of your laptop or cluster for parallelization. This means that if the `-n/--ncpu` option is not set the code will distribute the computational workload across all the available cores for efficient processing. However, if you want to manually specify the number of CPUs to be used, you can provide the `-n/--ncpu` option followed by the desired number of CPUs, for example, `-n 8` to use 8 CPUs. In this case, the code will parallelize the calculation by employing the specified number of CPUs. The purpose of parallelization is to accelerate the computation process by dividing the workload among multiple processors. By utilizing multiple CPUs, you can potentially reduce the overall processing time for calculating the RSD map. Note: The actual number of CPUs available for parallelization may depend on the hardware specifications of your system or the constraints set by your cluster environment.



## 5.5 - Output 

The output of the code includes two files:

* **`trace_${ProteinName}.txt`**: This file serves as a log or progress report during the execution of the code. It provides updates on the number of Resolution and Relevance points that have been calculated so far. This information helps monitor the progress of the calculation, especially if it takes a long time to complete. Additionally, it includes the time required to calculate a single point at each fixed number of retained sites. This timing information can be useful for performance analysis and optimization. The file may also estimate the remaining total time based on the current progress, giving you an idea of how much time is left for the calculation to complete.


* **`Hs-Hk-Nsites-${ProteinName}.txt`**: This file contains the actual results of the Resolution and Relevance calculations for different numbers of retained sites. It provides three rows of data. The first row corresponds to the values of Resolution ($H_s$); the second row contains the values of Relevance ($H_k$), whereas the third one specifies the corresponding number of retained sites for each Resolution and Relevance point. This information allows you to analyze the relationship between the number of retained sites and the quality of the calculated Resolution and Relevance values. By examining these values, you can identify the optimal number of retained sites that strikes a balance between capturing structural information and essential dynamics. (look **Section 6**). 

A short explaination of arguments is provided by launching the command `python3 ResRel-MPI.py -h` or `python3 ResRel-MPI.py --help`. Alternatively, for printing a short usage message, please type: `python3 ResRel-MPI.py` or `python3 ResRel-MPI.py -u`


# 6 - Hs-Hk-plot.py 

## 6.1 - Scope 
This code serves a dual purpose:

1. **Plotting data**. The code generates and saves four plots:

  * Resolution & Relevance Plot: This plot displays the Resolution (Hs) and Relevance (Hk) points, with the same color representing points obtained from the same number of retained sites but different mappings. Additionally, a zoomed-in region of interest is shown where the slope is close to -1, providing a detailed view of that area.
  * Zoomed Relevance & Resolution Curve (Hk-Hs): This plot focuses specifically on the region where the slope is -1, providing a closer look at the relationship between Relevance and Resolution in the region of our interest. 
  * Slope vs. Index: This plot shows the slope values plotted against an increasing index ranging from 1 to the total number of points. This visualization helps identify any patterns or trends in the slope values.
  * Histogram of Frequencies: This plot displays the frequency distribution of the number of sites with different occurrences. It provides insights into the distribution of retained sites and their frequencies.
  
2. **Optimal Number of Sites Calculation**. The code also calculates the optimal number of sites for a biomolecule based on an Resolution and relevance plot analysis. The goal is to minimize the loss of information after reducing the number of atoms. The result of this calculation is reported in the file "Opt-number-of-sites.txt". This information is valuable for determining the appropriate number of retained sites that balances the preservation of essential structural information with the reduction in computational complexity.


## 6.2 - Tasks 
The Relevance and Resolution plot is made up _N_ points (by default _N about 10000_ points). Thus, in order to find the optimal number of sites is necessary to "simplify" the ResRel curve. A TALE SCOPO average values of Resolution ($\bar H_s$) and Relevance ($\bar H_k$) are computed. This calculation can be done in two different ways:  

* **`density`**: The x-axes (Resolution Hs) is divided in _X_ intervals such that each of one contains the same number of points _D_ (by default _D = 100_). Thus, the interval lenght is not fixed. What is fixed is the number (density) of Hs-Hk points in each interval. Then, in each interval the average values for Hs and Hk are computed (Hs_avg and Hk_avg). If using this option (recommended) the calculation of the average values for Resolution (Hs) and Relevance (Hk) is based on the same density of points. Indeed, every N points (the default value is 100) the average calculation for Hs and Hk is performed. In this way, the computation of average values is more fair since the lenght of interval is not fixed, but it is chosen in order that the number of points (Hs-Hk) is always the same. The default value for N is 100 points; however, such value con be changed using the flag [-d]. 
                                        
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



