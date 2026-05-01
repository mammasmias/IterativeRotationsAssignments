
# compiler flags (for this dir and below)
# add_compile_options(-ffree-line-length-none -march=native -funroll-loops -ffast-math -Ofast)

# find lapack/blas
find_package(LAPACK REQUIRED)
## TODO: if(NOT FOUND LAPACK) compile and link local

# set module directory to root/include
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_SOURCE_DIR}/include)

set(MY_MODS
  version.f90
  ira_precision.f90
  timer.f90
  sofi_tools.f90
  err_module.f90
  sorting_module.f90
  set_candidate.f90
  cshda.f90
  ira_routines.f90
  sofi_routines.f90
	library_ira.f90
	library_sofi.f90
  CACHE INTERNAL "")

# specify shared library libira.so from MY_MODS
add_library(ira SHARED ${MY_MODS})
target_compile_options( ira
  PRIVATE
  -DSTRVERSION="${sem_version}" -DDATEVERSION=${date_version}
  $<$<COMPILE_LANGUAGE:Fortran>:-ffree-line-length-0>
  )

set_target_properties(ira PROPERTIES
  Fortran_PREPROCESS ON
  Fortran_FORMAT FREE
  POSITION_INDEPENDENT_CODE ON
  LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/lib"
  NO_SONAME 1
  UNITY_BUILD ON
  )

# link lapack
target_link_libraries(ira PUBLIC ${LAPACK_LIBRARIES})
