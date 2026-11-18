cmake_minimum_required(VERSION 3.21)

file(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/version.txt vv )
list(GET vv 0 sem_version)
list(GET vv 1 date_version)
message(STATUS "IRA version: ${sem_version}")

project(ira
  VERSION "${sem_version}"
	DESCRIPTION
  "Algorithms for Iterative Rotations and Assignments (IRA); Symmetry Operations FInder (SOFI); and Constrained Shortest Distance Assignment (CShDA)."
  LANGUAGES Fortran C)

# optimize for native machine
option(IRA_ENABLE_NATIVE_OPTIMIZATION
  "optimize IRA for the build machine (for example, '-march=native')" ON)
option(IRA_ENABLE_FAST_MATH
  "allow non-IEEE transformations (for example, '-Ofast -ffast-math')" ON)

# Prefer FlexiBLAS, then any BLAS/LAPACK implementation, and finally the
# collision-free bundled implementation.  Set this to FLEXIBLAS, LAPACK, or
# BUNDLED to require a particular path.
set(IRA_LINALG_PROVIDER "AUTO" CACHE STRING
  "linear algebra provider: AUTO, FLEXIBLAS, LAPACK, or BUNDLED")
set_property(CACHE IRA_LINALG_PROVIDER PROPERTY STRINGS
  AUTO FLEXIBLAS LAPACK BUNDLED)

# optionally build tests
option( IRA_BUILD_TESTS "build basic IRA tests" OFF )

option( IRA_INSTALL_PYTHON "install the python-specific files" OFF )
if(IRA_INSTALL_PYTHON)
  # shared lib is needed for python interf, managed by pyproject.toml
  set(BUILD_SHARED_LIBS ON)
endif()

# include variables for standard install directory names
include(GNUInstallDirs)

# source files
set(ira_src
  src/m_ira_version.f90
  src/ira_precision.f90
  src/timer.f90
  src/m_ira_tools.f90
  src/m_sofi_tools.f90
  src/m_ira_error.f90
  src/sorting_module.f90
  src/m_set_candidate.f90
  src/m_cshda.f90
  src/m_ira_routines.f90
  src/m_sofi_routines.f90
  src/ira_mod.f90
  src/library_ira.f90
  src/library_sofi.f90
  src/ira_external.f90
  src/sofi_external.f90
  )


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
check_fortran_compiler_flag( "-funroll-loops" unroll_loops )
if(unroll_loops)
  list(APPEND flagstr "-funroll-loops" )
endif()
#
if(IRA_ENABLE_FAST_MATH)
  check_fortran_compiler_flag( "-Ofast" Ofast )
  if(Ofast)
    list(APPEND flagstr "-Ofast" )
  else()
    check_fortran_compiler_flag( "-O3" O3 )
    if(O3)
      list(APPEND flagstr "-O3" )
    endif()
  endif()
  check_fortran_compiler_flag( "-ffast-math" fast_math )
  if(fast_math)
    list(APPEND flagstr "-ffast-math" )
  endif()
else()
  check_fortran_compiler_flag( "-O3" O3 )
  if(O3)
    list(APPEND flagstr "-O3" )
  endif()
endif()
#
if( IRA_ENABLE_NATIVE_OPTIMIZATION )
  check_fortran_compiler_flag("-march=native" march_native )
  if(march_native)
    list(APPEND flagstr "-march=native")
  else()
    check_fortran_compiler_flag("-xHost" xhost )
    if(xhost)
      list(APPEND flagstr "-xHost")
    endif()
  endif()
endif()
message(STATUS "IRA :: Using compiler flags: ${flagstr}")

# define ira lib
set(ira_includedir "${CMAKE_CURRENT_BINARY_DIR}/include")
add_library(${PROJECT_NAME} ${ira_src})
add_library(ira::ira ALIAS ira)  # alias target name
set_target_properties( ${PROJECT_NAME}
  PROPERTIES
  Fortran_PREPROCESS ON
  Fortran_FORMAT FREE
  POSITION_INDEPENDENT_CODE ON
  Fortran_MODULE_DIRECTORY "${ira_includedir}"
  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib" # for shared lib
  ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib" # for static lib, i.e. 'archive'
  OUTPUT_NAME "${PROJECT_NAME}"
  VERSION "${PROJECT_VERSION}"
  SOVERSION "${PROJECT_VERSION_MAJOR}"
  )
# compile flags
target_compile_options( ${PROJECT_NAME}
  PRIVATE
  ${flagstr}
  )
# preprocessor defs
target_compile_definitions( ${PROJECT_NAME}
  PRIVATE
  STRVERSION="${sem_version}" DATEVERSION=${date_version}
  )
# propagate include dirs for different callers:
# PUBLIC means the dir is used for compilation, and has to be included when linking
# INTERFACE means the dir is not used for compilation, but has to be included when linking
target_include_directories( ${PROJECT_NAME}
  PUBLIC
  $<BUILD_INTERFACE:${ira_includedir}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/interface>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}>
  )
# string for targets
set(alltargets "")
list(APPEND alltargets "${PROJECT_NAME}" )


string(TOUPPER "${IRA_LINALG_PROVIDER}" _ira_linalg_provider)
set(_ira_linalg_choices AUTO FLEXIBLAS LAPACK BUNDLED)
if(NOT _ira_linalg_provider IN_LIST _ira_linalg_choices)
  message(FATAL_ERROR
    "IRA_LINALG_PROVIDER must be AUTO, FLEXIBLAS, LAPACK, or BUNDLED")
endif()

set(_ira_external_linalg FALSE)
set(IRA_PC_PRIVATE_LIBS "")
set(IRA_CONFIG_LINALG_PROVIDER "BUNDLED")
if(_ira_linalg_provider STREQUAL "AUTO"
    OR _ira_linalg_provider STREQUAL "FLEXIBLAS")
  if(DEFINED BLA_VENDOR)
    set(_ira_saved_bla_vendor "${BLA_VENDOR}")
    set(_ira_had_bla_vendor TRUE)
  else()
    set(_ira_had_bla_vendor FALSE)
  endif()
  set(BLA_VENDOR FlexiBLAS)
  find_package(BLAS QUIET)
  find_package(LAPACK QUIET)
  if(BLAS_FOUND AND LAPACK_FOUND)
    set(_ira_external_linalg TRUE)
    set(_ira_linalg_name "FlexiBLAS")
    set(IRA_PC_PRIVATE_LIBS "-lflexiblas")
    set(IRA_CONFIG_LINALG_PROVIDER "FLEXIBLAS")
  elseif(_ira_linalg_provider STREQUAL "FLEXIBLAS")
    message(FATAL_ERROR
      "FlexiBLAS was requested but its BLAS and LAPACK interfaces were not found")
  endif()
  if(_ira_had_bla_vendor)
    set(BLA_VENDOR "${_ira_saved_bla_vendor}")
  else()
    unset(BLA_VENDOR)
  endif()
endif()

if(NOT _ira_external_linalg
    AND (_ira_linalg_provider STREQUAL "AUTO"
         OR _ira_linalg_provider STREQUAL "LAPACK"))
  unset(BLAS_FOUND)
  unset(BLAS_LIBRARIES)
  unset(LAPACK_FOUND)
  unset(LAPACK_LIBRARIES)
  find_package(BLAS QUIET)
  find_package(LAPACK QUIET)
  if(BLAS_FOUND AND LAPACK_FOUND)
    set(_ira_external_linalg TRUE)
    set(_ira_linalg_name "BLAS/LAPACK")
    set(IRA_PC_PRIVATE_LIBS "-llapack -lblas")
    set(IRA_CONFIG_LINALG_PROVIDER "LAPACK")
  elseif(_ira_linalg_provider STREQUAL "LAPACK")
    message(FATAL_ERROR
      "BLAS/LAPACK was requested but both interfaces were not found")
  endif()
endif()

if(_ira_external_linalg)
  message(STATUS
    "IRA :: Using ${_ira_linalg_name}: ${LAPACK_LIBRARIES};${BLAS_LIBRARIES}")
  target_compile_definitions(${PROJECT_NAME} PRIVATE IRA_USE_EXTERNAL_LAPACK)
  target_link_libraries(${PROJECT_NAME} PRIVATE LAPACK::LAPACK BLAS::BLAS)
else()
  # Every bundled symbol is prefixed ira_, so this fallback cannot interpose on
  # another library's BLAS. Keep dlamch/disnan IEEE-safe and position-independent.
  message(STATUS "IRA :: Using bundled prefixed BLAS/LAPACK fallback")
  add_library(lap_local OBJECT "${CMAKE_CURRENT_SOURCE_DIR}/lap_local/lap.f")
  target_compile_options(lap_local PRIVATE -O2)
  set_target_properties(lap_local PROPERTIES
    Fortran_FORMAT FIXED
    POSITION_INDEPENDENT_CODE ON
    )
  target_link_libraries(${PROJECT_NAME} PRIVATE lap_local)
  # Static libira records this object library as a usage requirement.
  list(APPEND alltargets "lap_local")
endif()

get_target_property(IRA_LIB_TYPE ${PROJECT_NAME} TYPE)
if(IRA_LIB_TYPE STREQUAL "SHARED_LIBRARY")
  add_custom_command(
    TARGET ${PROJECT_NAME} POST_BUILD
    COMMAND /bin/sh "${CMAKE_CURRENT_SOURCE_DIR}/tests/link_external_blas.sh"
            --symbols-only $<TARGET_FILE:${PROJECT_NAME}>
    VERBATIM)
endif()


if(IRA_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()




# install procedures
if(IRA_INSTALL_PYTHON)
  # Install rules for Python packaging (scikit-build-core):
  # name of the directory is the name of python package
  install(FILES
    "interface/__init__.py"
    "interface/ira_mod.py"
    DESTINATION "ira_mod" )
  install(TARGETS
    ${alltargets}
    EXPORT
    "${PROJECT_NAME}-targets"
    LIBRARY DESTINATION  "ira_mod/lib"
    INCLUDES DESTINATION "ira_mod/include"
    )

  # install( FILES $<TARGET_FILE:ira>
  #   DESTINATION ira_mod/lib
  # )
else()
  # Regular install

  # which targets to install, add for export
  install(TARGETS
    ${alltargets}
    EXPORT
    "${PROJECT_NAME}-targets"
    LIBRARY DESTINATION  "${CMAKE_INSTALL_LIBDIR}"
    ARCHIVE DESTINATION  "${CMAKE_INSTALL_LIBDIR}"
    INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}"
    )
  # the module directory
  install(DIRECTORY
    "${ira_includedir}/"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}"
    )
  # the C interface header
  install(FILES
    "interface/iralib_interf.h"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}"
    )
  # export the targets
  install(EXPORT
    "${PROJECT_NAME}-targets"
    NAMESPACE "${PROJECT_NAME}::"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
    )
  # write the config file for find_package()
  include(CMakePackageConfigHelpers)
  configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/template.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config.cmake"
    INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
    )
  # write the version file
  write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake"
    VERSION "${PROJECT_VERSION}"
    COMPATIBILITY SameMinorVersion
    )
  # export the config & version files
  install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
    )
  # configure the pkg-config template
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/template.pc"
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.pc"
    @ONLY
    )
  # install ira.pc
  install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.pc"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig"
  )

  # export config file during build, for include()
  export(EXPORT
    "${PROJECT_NAME}-targets"
    FILE "${CMAKE_CURRENT_BINARY_DIR}/cmake/${PROJECT_NAME}/${PROJECT_NAME}-BuildConfig.cmake"
    NAMESPACE "${PROJECT_NAME}::"
    )
endif()

