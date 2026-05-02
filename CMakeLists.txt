cmake_minimum_required(VERSION 3.21)

file(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/version.txt vv )
list(GET vv 0 sem_version)
list(GET vv 1 date_version)
message(STATUS "IRA version: ${sem_version}")

project(IRA
  VERSION "${sem_version}"
	DESCRIPTION
  "Algorithms for Iterative Rotations and Assignments (IRA); Symmetry Operations FInder (SOFI); and Constrained Shortest Distance Assignment (CShDA)."
  LANGUAGES Fortran C)




# find lapack/blas
find_package(LAPACK REQUIRED)
## TODO: if(NOT FOUND LAPACK) compile and link local

# set module directory to root/include
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include")

# sources
set(ira_src
  src/version.f90
  src/ira_precision.f90
  src/timer.f90
  src/sofi_tools.f90
  src/err_module.f90
  src/sorting_module.f90
  src/set_candidate.f90
  src/cshda.f90
  src/ira_routines.f90
  src/sofi_routines.f90
	src/library_ira.f90
	src/library_sofi.f90
  CACHE INTERNAL "")



include(CheckFortranCompilerFlag)

# find some available flags
set(flagstr "")
#
check_fortran_compiler_flag( "-ffree-line-length-0" ff0 )
if(ff0)
  list(APPEND flagstr "-ffree-line-length-0" )
else()
endif()
#
check_fortran_compiler_flag("-march=native" march_native )
if(march_native)
  list(APPEND flagstr "-march=native")
else()
  check_fortran_compiler_flag("-xHost" xhost )
  if(xhost)
    list(APPEND flagstr "-xHost")
  endif()
endif()
#
check_fortran_compiler_flag( "-funroll-loops" unroll_loops )
if(unroll_loops)
  list(APPEND flagstr "-funroll-loops" )
endif()
#
check_fortran_compiler_flag( "-ffast-math" fast_math )
if(fast_math)
  list(APPEND flagstr "-ffast-math" )
endif()
#
check_fortran_compiler_flag( "-Ofast" Ofast )
if(Ofast)
  list(APPEND flagstr "-Ofast" )
endif()
message(STATUS "Using compiler flags: ${flagstr}")


# specify shared library libira
add_library(ira SHARED ${ira_src})
target_compile_options( ira
  PRIVATE
  -DSTRVERSION="${sem_version}" -DDATEVERSION=${date_version}
  $<$<COMPILE_LANGUAGE:Fortran>:${flagstr}>
  )

set_target_properties(ira PROPERTIES
  Fortran_PREPROCESS ON
  Fortran_FORMAT FREE
  POSITION_INDEPENDENT_CODE ON
  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib"
  )

# link lapack
target_link_libraries(ira PUBLIC ${LAPACK_LIBRARIES})



# Install rules for Python packaging (scikit-build-core)
install(TARGETS ira LIBRARY DESTINATION ira_mod/lib)
install(FILES
    interface/__init__.py
    interface/ira_mod.py
    DESTINATION ira_mod
)

