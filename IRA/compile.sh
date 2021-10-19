gfortran -fdefault-real-8 -c sorting_module.f90
gfortran -fdefault-real-8 -c cshda.f90
gfortran -fdefault-real-8 -c ira_routines.f90
gfortran -fdefault-real-8 -c read_typ.f90
gfortran -fdefault-real-8 -o ira_eq.x sorting_module.o cshda.o ira_routines.o read_typ.o ira_eq.f90 -llapack
gfortran -fdefault-real-8 -o ira_noneq.x sorting_module.o cshda.o ira_routines.o read_typ.o ira_noneq.f90 -llapack
gfortran -fdefault-real-8 -o ira_general.x sorting_module.o cshda.o ira_routines.o read_typ.o ira_general.f90 -llapack
ar -rcv lib_ira.a sorting_module.o cshda.o ira_routines.o
