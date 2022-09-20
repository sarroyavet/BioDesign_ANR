# Bio--inspired shape--design algorithm for two--dimensional contact profiles under frictionless contact conditions.
This algorithm generates contact profiles with uniform contact pressure
subjected to maximum von Mises shear restrictions.
It was developed using the finite element package 
[CODE_ASTER 14.4](https://code-aster.org/) and the meshing algorithm
[GMSH 4.9.3](https://gmsh.info/doc/texinfo/gmsh.html).

## Requisites
    Code_Aster 14.4 (and all its prerequisites)
    Gmsh 4.9.3
    Python 3.8.10 and libraries os, sys, json, getopt, time, gmsh, math, numpy

## Usage
For the execution of the algorithm use in the terminal:

    python3 RunAnalysis.py -i parameter_file.json

In parameter_file.json, all the parameters necessary are provided.
An example of it can be found in [PARAMS.json](Codes/PARAMS.json).
An overview of the codes is provided in
[Scripts.pdf](Documentation/StableVersion/Scripts.pdf)

### Warning
In line 12 of [BioInspiredGrowth.dumm](Codes/BioInspiredGrowth.dumm),
the path of the folder [PythonUtilities](Codes/PythonUtilities) in
your computer must be provided.

## Author
David Hernández--Aristizábal
[ORCiD](https://orcid.org/0000-0003-4041-3785)
