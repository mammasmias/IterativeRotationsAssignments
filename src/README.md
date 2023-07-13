# Description

<!-- The shape matching algorithm Iterative Rotations and Assignments (IRA) is described in the publication -->
<!-- [[1]](#1). It is also the main subject of the dissertation [[2]](#2), where a workflow -->
<!-- inserting IRA into an off-lattice kMC algorithm is developed. -->
This directory contains the source code to IRA, CShDA, and SOFI algorithms, implemented in Fortran.
Further descriptions of these algorithms are given in [[1]](#1), [[2]](#2), and [[3]](#3).

The default compilation also produces a shared library `libira.so`, which contains wrapper routines to the main IRA, CShDA, and SOFI routines, such that they are bound to the C-namespace via `bind(C)` and `use iso_c_binding`.

# IRA: Iterative Rotations and Assignments

The Iterative Rotations and Assignments (IRA) algorithm is a shape matching
algorithm for matching generic atomic structures, including structures with
different number of atoms.

# CShDA: Constrained Shortest Distance Assignment

The Constrained Shortest Distance Assignment (CShDA) algorithm is a Linear
Assignment Problem (LAP) solver, which is used by the IRA algorithm.

# SOFI: Symmetry Operations FInder

The Symmetry Operations FInder (SOFI) algorithm is a modification of IRA, that
is designed to find the symmetry operations of a given structure, by solving the
degenerate shape-matching problem. It also uses the CShDA algorithm.

## Files

Source files:

 - cshda.f90          : the CShDA routine for pbc and non-pbc cases;
 - ira_routines.f90   : all the routines specific to IRA;
 - set_candidate.f90  : routines for selecting candidate central atoms and their vectors;
 - read_typ.f90       : a function to read atomic types from xyz format;
 - sorting_module.f90 : the mergesort algorithm;
 - sofi_tools.f90     : misc routines and definitions for SOFI;
 - sofi_routines      : all routines specific to SOFI;
 - library_ira.f90    : C-bound wrappers to the main routines of IRA, compiles into `libira.so`;
 - library_sofi.f90   : C-bound wrappers to the main routines of SOFI, compiles into `libira.so`

## Compile

TO COMPILE: (you need lapack library, see the Makefile)

    make all


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

<a id="3">[3]</a>
Gunde M., Salles N., Grisanti L., Hemeryck A., Martin Samos L.
*SOFI: Finding point group symmetries in atomic clusters as finding the set of degenerate solutions in a shape-matching problem*,
Journal of Chemical Information and Modeling (submitted)
