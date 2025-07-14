# set module directory to build/Obj
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/Obj)

set(MY_MODS   ira_precision.f90
              timer.f90
              sofi_tools.f90
              err_module.f90
              version.f90
              sorting_module.f90
              set_candidate.f90
              cshda.f90 
              ira_routines.f90 
              sofi_routines.f90
	      library_ira.f90
	      library_sofi.f90
              CACHE INTERNAL "")

# specify shared library libira.so from MY_MODS
add_library(ira SHARED "${MY_MODS}")

# link lapack
target_link_libraries(ira ${LAPACK_LIBRARIES})
