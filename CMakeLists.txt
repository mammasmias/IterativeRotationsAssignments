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

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 9.0)
    message(FATAL_ERROR
      "gfortran >= 9.0 is required, found: ${CMAKE_Fortran_COMPILER_VERSION}
       If executing the command `gfortran --version` returns a different version,
       you can try to specify it manually to cmake:
          rm -rf <builddir>
          cmake -B <builddir> -DCMAKE_Fortran_COMPILER=$(which gfortran)
       Stopping.")
  endif()
endif()

# optionallt build tests
option( BUILD_TESTS "build basic IRA tests" OFF )

# specify static library libira
option( BUILD_SHARED_LIBS "build shared libraries" OFF )

option( INSTALL_PYTHON "install the python-specific files" OFF )
if(INSTALL_PYTHON)
  # shared lib is needed for python interf, managed by pyproject.toml
  set(BUILD_SHARED_LIBS ON)
endif()

# set module directory to root/include
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include")

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


# define ira lib
add_library(ira ${ira_src})
set_target_properties(ira PROPERTIES
  Fortran_PREPROCESS ON
  Fortran_FORMAT FREE
  POSITION_INDEPENDENT_CODE ON
  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib" # for shared lib
  ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib" # for static lib, i.e. 'archive'
  )
# compile flags
target_compile_options( ira PRIVATE
  $<$<COMPILE_LANGUAGE:Fortran>:${flagstr}>
  )
# preprocessor defs
target_compile_definitions( ira PRIVATE
  -DSTRVERSION="${sem_version}" -DDATEVERSION=${date_version}
  )
# propagate include dirs for different callers:
# PUBLIC means the dir is used for compilation, and has to be included when linking
# INTERFACE means the dir is not used for compilation, but has to be included when linking
target_include_directories(ira
  PUBLIC    $<$<COMPILE_LANGUAGE:Fortran>:${CMAKE_Fortran_MODULE_DIRECTORY}>
  INTERFACE $<$<COMPILE_LANGUAGE:C>:${CMAKE_CURRENT_SOURCE_DIR}/interface>
  INTERFACE $<$<COMPILE_LANGUAGE:CXX>:${CMAKE_CURRENT_SOURCE_DIR}/interface>
  )


# find lapack/blas
find_package(LAPACK)

# link lapack
if(${LAPACK_FOUND})
  target_link_libraries(ira PUBLIC ${LAPACK_LIBRARIES})
else()
  # compile and link the local lapack files
  message(STATUS "IRA :: LAPACK library not found. Compiling local version (un-optimised).")
  add_library(lap_local OBJECT ${CMAKE_CURRENT_SOURCE_DIR}/lap_local/lap.f)
  target_compile_options(lap_local PRIVATE -O0)
  set_target_properties(lap_local PROPERTIES Fortran_FORMAT FIXED POSITION_INDEPENDENT_CODE ON)
  target_link_libraries(ira lap_local)
endif()


if(BUILD_TESTS)
  add_subdirectory(tests)
endif()


if(INSTALL_PYTHON)
  # Install rules for Python packaging (scikit-build-core):
  # name of the directory is the name of python package
  install(FILES
    interface/__init__.py
    interface/ira_mod.py
    DESTINATION ira_mod )
  install( FILES $<TARGET_FILE:ira>
    DESTINATION ira_mod/lib
    )
else()
  # Regular install
  include(GNUInstallDirs)
  install(TARGETS ira
    EXPORT IRATargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  install(FILES
    interface/iralib_interf.h
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
endif()

