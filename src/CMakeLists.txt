
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


include(CheckFortranCompilerFlag)

# find some available flags
set(flagstr "")
#
check_fortran_compiler_flag( "-ffree-line-length-0" ff0 )
if(ff0)
  list(APPEND flagstr ";-ffree-line-length-0" )
else()
endif()
#
check_fortran_compiler_flag("-march=native" march_native )
if(march_native)
  list(APPEND flagstr ";-march=native")
else()
  check_fortran_compiler_flag("-xHost" xhost )
  if(xhost)
    list(APPEND flagstr ";-xHost")
  endif()
endif()
#
check_fortran_compiler_flag( "-funroll-loops" unroll_loops )
if(unroll_loops)
  list(APPEND flagstr ";-funroll-loops" )
endif()
#
check_fortran_compiler_flag( "-ffast-math" fast_math )
if(fast_math)
  list(APPEND flagstr ";-ffast-math" )
endif()
#
check_fortran_compiler_flag( "-Ofast" Ofast )
if(Ofast)
  list(APPEND flagstr ";-Ofast" )
endif()
message(STATUS ${flagstr})


# specify shared library libira.so from MY_MODS
add_library(ira SHARED ${MY_MODS})
target_compile_options( ira
  PRIVATE
  -DSTRVERSION="${sem_version}" -DDATEVERSION=${date_version}
  $<$<COMPILE_LANGUAGE:Fortran>:${flagstr}>
  # $<$<COMPILE_LANGUAGE:Fortran>:-march=native>
  # $<$<COMPILE_LANGUAGE:Fortran>:-funroll-loops>
  # $<$<COMPILE_LANGUAGE:Fortran>:-ffast-math>
  # $<$<COMPILE_LANGUAGE:Fortran>:-Ofast>
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


