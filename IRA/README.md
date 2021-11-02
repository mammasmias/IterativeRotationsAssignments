# Description
This directory contains the IRA and CShDA algorithms implemented in Fortran.
Further description of these algorithms is given in [[1]](#1). 

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

 - read_typ.f90 contains a function to read atomic types from xyz format.

 - sorting_module.f90 contains the mergesort algorithm.

Example program files:

 - ira_eq.f90 contains the main program for IRA shape matching when two
   structures contain the same number of atoms.

 - ira_noneq.f90 contains the main program for IRA shape matching when two
   structures contain a different number of atoms.

 - ira_general.f90 contains the main program for IRA which decides which of the
   above cases to do, based on number of atoms.

# Compile and run
TO COMPILE: (you need lapack library)

   `sh compile.sh`


TO RUN:

   `./ira_eq.x      < input_file_eq.xyz      > output_file`

   `./ira_noneq.x   < input_file_noneq.xyz   > output_file`

   `./ira_general.x < input_file_general.xyz > output_file`

The input file should contain the two structures to be matched, both in the xyz
format. The input_file_eq.xyz should contain structures with equal number of
atoms. The input_file_noneq.xyz should contain structures with different number
of atoms, such that the first structure contains less atoms that the second
structure. The input_file_general.xyz con be either of the above two cases.


## References
<a id="1">[1]</a> 
Gunde M., Salles N., Hemeryck A., Martin Samos L.
*IRA: A shape matching approach for recognition and comparison of generic atomic patterns*,
Journal of Chemical Information and Modeling (2021), DOI: [https://doi.org/10.1021/acs.jcim.1c00567](https://doi.org/10.1021/acs.jcim.1c00567), HAL: [hal-03406717](https://hal.laas.fr/hal-03406717), arXiv: [2111.00939](https://export.arxiv.org/abs/2111.00939)

