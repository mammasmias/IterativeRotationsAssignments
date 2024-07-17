# Description

This directory contains three example programs which utilise the SOFI algorithm:

 `f90_program.f90` is an example program in Fortran;
 `c_program.c` is an example program in C;
 `python_program.py` is an example program in python3.

> NOTE: these are only simple examples intended to show how SOFI could be called from either of the three languages, i.e. they show how to manage the memory, how to compile the caller programs, and show some simple calls to the library. In order to change the input/output format, the programs should be modified as needed.

# Compile and run

The f90 and C programs need to be compiled, type (see `Makefile`):

    make all

The f90 program is run as follows:

    ./f90_program.x < example_inputs/"filename"

The threshold `sym_thr` is hard-coded in `f90_program.f90`, thus the program should be recompiled upon changing the value.

The C program is run as:

    ./c_program.x < example_inputs/"filename"

In the `c_program.c`, the threshold `sym_thr` is likewise hard-coded.

The python program needs to know the location of the `ira_mod` module, to do this type (you might also put this into your `.bashrc`):

    export PYTHONPATH=$PYTHONPATH:/path/to/IterativeRotationsAssignments/interface

with the correct path to the `/interface` directory where `ira_mod.py` is located.
Then run the python program as:

    python3 python_program.py example_inputs/"filename"

# Example inputs

The `example_inputs` directory contains some structures written in .xyz format. The C program expects
a structure with atomic types given as integers, while the f90 and python programs can take also strings.
Note that not all example inputs are compatible to the example C program.

The structure `S6_D3d.xyz` has some inexact symmetries, if you set the `sym_thr=0.02` you should obtain the S6 point group, and if you set `sym_thr=0.03` it should be D3d. Likewise for some other structures, changing `sym_thr` can change the outcome.
