
# compiler flags (for this dir and below)
add_compile_options(-ffree-line-length-none -march=native -funroll-loops -ffast-math -Ofast)


# set module directory to root/include
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_SOURCE_DIR}/include)

set(MY_MODS
  ${CMAKE_CURRENT_SOURCE_DIR}/ira_precision.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/timer.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/sofi_tools.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/err_module.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/version.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/sorting_module.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/set_candidate.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/cshda.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/ira_routines.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/sofi_routines.f90
	${CMAKE_CURRENT_SOURCE_DIR}/library_ira.f90
	${CMAKE_CURRENT_SOURCE_DIR}/library_sofi.f90
  CACHE INTERNAL "")

# specify shared library libira.so from MY_MODS
add_library(ira SHARED "${MY_MODS}")
set_target_properties(ira PROPERTIES
  Fortran_PREPROCESS ON
  Fortran_FORMAT FREE
  POSITION_INDEPENDENT_CODE ON
  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/lib"
  NO_SONAME 1
  )

# link lapack
target_link_libraries(ira ${LAPACK_LIBRARIES})
