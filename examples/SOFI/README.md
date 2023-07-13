# Description

This directory contains three programs which utilise the SOFI algorithm:

 `f90_program.f90` is an example program in Fortran;
 `c_program.c` is an example program in C;
 `python_program.py` is an example program in python3.



# Compile and run

The f90 and C programs need to be compiled, type (see `Makefile`):

    make all

The f90 program is run as follows:

    ./f90_program.x < example_inputs/"filename"

The C program is run as:

    ./c_program.x < example_inputs/"filename"

The python program needs to know the location of the `ira_mod` module, to do this type (you might also put this into your `.bashrc`):

    export PYTHONPATH=$PYTHONPATH:/path/to/IterativeRotationsAssignments/interface

with the correct path to the `/interface` directory where `ira_mod.py` is located.
Then run the python program as:

    python3 python_program.py example_inputs/"filename"

# Example inputs

The `example_inputs` directory contains some structures written in .xyz format. The f90 and C programs expect
a structure with atomic types given as integers, while the python program can take also strings.

The structure `S6_D3d.xyz` has some inexact symmetries, if you set the `sym_thr=0.02` you should obtain the S6 point group, and if you set `sym_thr=0.03` it should be D3d.
