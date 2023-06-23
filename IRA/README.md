# Description

This is the version 1.5.0 of the shape matching algorithm 
Iterative Rotations and Assignments (IRA), described in the publication
[[1]](#1). It is also the main subject of the dissertation [[2]](#2), where a workflow
inserting IRA into an off-lattice kMC algorithm is developed.
This directory contains the IRA and CShDA algorithms implemented in Fortran.
Further descriptions of these algorithms are given in [[1]](#1) and [[2]](#2). 

The default compilation also produces a shared library `shlib_ira.so`, which contains wrapper routines to the main IRA and CShDA routines, such that they are bound to the C-namespace via `bind(C)` and `use iso_c_binding`.

# Iterative Rotations and Assignments (IRA)

The Iterative Rotations and Assignments (IRA) algorithm is a shape matching
algorithm for matching generic atomic structures, including structures with
different number of atoms.

# Constrained Shortest Distance Assignment (CShDA)

The Constrained Shortest Distance Assignment (CShDA) algorithm is a Linear
Assignment Problem (LAP) solver, which is used by the IRA algorithm. 

## Files

Routine files:

 - cshda.f90          : the CShDA routine for pbc and non-pbc cases;
 - ira_routines.f90   : all the routines specific to IRA;
 - set_candidate.f90  : routines for selecting candidate central atoms and their vectors;
 - read_typ.f90       : a function to read atomic types from xyz format;
 - sorting_module.f90 : the mergesort algorithm;
 - library_ira.f90    : C-bound wrappers to the main routines of IRA, compiles into `shlib_ira.so`;
 - ira_mod.py         : python module which calls the interface routines from `shlib_ira.so`, via the `ctypes` module.

Example program files:

 - f90_program.f90 is an example of a f90 program that calls main IRA routine, followed by application of SVD.
 - python_program.py demonstrates some simple use of the `ira_mod` python module, including matching structures with equal and nonequal numbers of atoms.

Example structures:

 - example_inputs/ folder contains some example input files, used for the f90, C, and python3 example programs.

## Compile and run

TO COMPILE: (you need lapack library, see the Makefile)

    make all

To run fortran example program:

    ./f90_program.x   <   example_inputs/input1.xyz
    ./f90_program.x   <   example_inputs/input2.xyz

to run C example program:

    ./c_program.x   example_inputs/input1.xyz
    ./c_program.x   example_inputs/input2.xyz

to run python3 example program:

    python3  python_program.py


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
