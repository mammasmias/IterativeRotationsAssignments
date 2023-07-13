# Description

This directory contains three programs which utilise the IRA/CShDA algorithms:

 `f90_program.f90` is an example program in Fortran;
 `c_program.c` is an example program in C;
 `python_program.py` is an example program in python3.



# Compile and run

The f90 and C programs need to be compiled, type (see `Makefile`):

    make all

The f90 program is run as follows:

    ./f90_program.x < example_inputs/input1.xyz
    ./f90_program.x < example_inputs/input2.xyz
    ./f90_program.x < example_inputs/input3.xyz

The C program is run as:

    ./c_program.x  example_inputs/input1.xyz
    ./c_program.x  example_inputs/input2.xyz
    ./c_program.x  example_inputs/input3.xyz

The python program needs to know the location of the `ira_mod` module, to do this type (you might also put this into your `.bashrc`):

    export PYTHONPATH=$PYTHONPATH:/path/to/IterativeRotationsAssignments/interface

with the correct path to the `/interface` directory where `ira_mod.py` is located.
Then run the python program as:

    python3 python_program.py
