gfortran -fdefault-real-8 -c read_typ.f90
gfortran -fdefault-real-8 -c randomize.f90
gfortran -fdefault-real-8 -o randomize.x read_typ.o randomize.o
