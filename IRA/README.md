# Description

This is the version 1.5.0 of the shape matching algorithm 
Iterative Rotations and Assignments (IRA), described in the publication
[[1]](#1). It is also the main subject of the dissertation [[2]](#2), where a workflow
inserting IRA into an off-lattice kMC algorithm is developed.
This directory contains the IRA and CShDA algorithms implemented in Fortran.
Further descriptions of these algorithms are given in [[1]](#1) and [[2]](#2). 

# Iterative Rotations and Assignments (IRA)
The Iterative Rotations and Assignments (IRA) algorithm is a shape matching
algorithm for matching generic atomic structures, including structures with
different number of atoms.

# Constrained Shortest Distance Assignment (CShDA)
The Constrained Shortest Distance Assignment (CShDA) algorithm is a Linear
Assignment Problem (LAP) solver, which is used by the IRA algorithm. 

# Files
Routine files:

 - cshda.f90 contains the CShDA routine for pbc and non-pbc cases.

 - ira_routines.f90 contains all the routines specific to IRA.
 
 - set_candidate.f90 contains routines for selecting candidate central atoms and their vectors.

 - read_typ.f90 contains a function to read atomic types from xyz format.

 - sorting_module.f90 contains the mergesort algorithm.

Example program files:

 - f90_program.f90 is an example of a f90 program that calls main IRA routine, followed by application of SVD.
 
 - python_program.py is an example of a python3 program that calls the routine which is a wrapper for the main IRA call and SVD in one call.
 
Example structures:

 - example_inputs/ folder contains some example input files that can be used for the f90 and python3 programs.
 
# Compile and run
TO COMPILE: (you need lapack library, see the Makefile)

   `make all`


TO RUN:

    `f90_program.x   <   example_inputs/f90_input1.xyz`
or:

    `python3   python_program.py`


## References
<a id="1">[1]</a> 
Gunde M., Salles N., Hemeryck A., Martin Samos L.
*IRA: A shape matching approach for recognition and comparison of generic atomic patterns*,
Journal of Chemical Information and Modeling (2021), DOI: [https://doi.org/10.1021/acs.jcim.1c00567](https://doi.org/10.1021/acs.jcim.1c00567), HAL: [hal-03406717](https://hal.laas.fr/hal-03406717), arXiv: [2111.00939](https://export.arxiv.org/abs/2111.00939)

<a id="2">[2]</a>
Gunde M.: *Development of IRA: a shape matching algorithm, its implementation
and utility in a general off-lattice kMC kernel*, PhD dissertation,
November 2021.
[PDF link](http://thesesups.ups-tlse.fr/5109/1/2021TOU30132.pdf) 
