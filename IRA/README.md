# Description
This directory contains the IRA and CShDA algorithms implemented in Fortran.

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
