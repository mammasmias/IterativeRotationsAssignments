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
  "optimize IRA for native CPU (this adds flags: '-march=native -ffast-math' or equivalent)" ON)

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
check_fortran_compiler_flag( "-Ofast" Ofast )
if(Ofast)
  list(APPEND flagstr "-Ofast" )
endif()
#
if( IRA_ENABLE_NATIVE_OPTIMIZATION )
  check_fortran_compiler_flag( "-ffast-math" fast_math )
  if(fast_math)
    list(APPEND flagstr "-ffast-math" )
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


# find lapack/blas
find_package(LAPACK)
# set(LAPACK_FOUND False) ### test lapack not found

# string for pkgconfig link to lapack in template.pc
set(LAPACKLIB_STR "")
# The bundled lap.f is an unoptimized reference BLAS. Exporting its
# dgemm_ makes a later dlopen bind that symbol instead of a tuned BLAS.
# Hiding is the default. -DIRA_HIDE_BUNDLED_BLAS=OFF keeps the old export.
# Callers that resolved BLAS from that export must link a BLAS themselves.
# GNU ld gets src/libira.map. A linker that rejects --version-script does
# not keep the export: lap.f is archived and linked with --exclude-libs.
# gfortran -fvisibility=hidden does not mark those procedures hidden, and
# objcopy --localize-symbol leaves them in .dynsym (binutils 2.42).
# If both linker flags are rejected, a shared build that would compile
# the bundled BLAS stops at configure.
option(IRA_HIDE_BUNDLED_BLAS
  "Hide reference BLAS symbols in the shared library. Callers that used the bundled export must link a BLAS themselves."
  ON)
set(IRA_VERSION_SCRIPT_APPLIED FALSE)
set(IRA_HIDE_WITH_EXCLUDE_LIBS FALSE)
set(IRA_BLAS_HIDDEN FALSE)
if(IRA_HIDE_BUNDLED_BLAS)
  get_target_property(IRA_LIB_TYPE ${PROJECT_NAME} TYPE)
  if(IRA_LIB_TYPE STREQUAL "SHARED_LIBRARY")
    include(CheckLinkerFlag)
    set(IRA_VERSION_SCRIPT "${CMAKE_CURRENT_SOURCE_DIR}/src/libira.map")
    check_linker_flag(Fortran
      "-Wl,--version-script=${IRA_VERSION_SCRIPT}"
      IRA_LINKER_SUPPORTS_VERSION_SCRIPT)
    if(IRA_LINKER_SUPPORTS_VERSION_SCRIPT)
      target_link_options(${PROJECT_NAME} PRIVATE
        "-Wl,--version-script=${IRA_VERSION_SCRIPT}")
      set(IRA_VERSION_SCRIPT_APPLIED TRUE)
      set(IRA_BLAS_HIDDEN TRUE)
    else()
      check_linker_flag(Fortran
        "-Wl,--exclude-libs,liblap_local.a"
        IRA_LINKER_SUPPORTS_EXCLUDE_LIBS)
      if(IRA_LINKER_SUPPORTS_EXCLUDE_LIBS)
        set(IRA_HIDE_WITH_EXCLUDE_LIBS TRUE)
        set(IRA_BLAS_HIDDEN TRUE)
      endif()
    endif()
  endif()
endif()

# link lapack
if(${LAPACK_FOUND})
  target_link_libraries( ${PROJECT_NAME} PUBLIC ${LAPACK_LIBRARIES})
  # the "-l" prefix is set here, so the @LAPACK_STR@ in template.pc is empty when local lapack is used
  set(LAPACKLIB_STR "-l${LAPACK_LIBRARIES}")
else()
  # compile and link the local lapack files
  message(STATUS "IRA :: LAPACK library not found. Compiling local version (un-optimised).")
  if(IRA_HIDE_BUNDLED_BLAS AND IRA_LIB_TYPE STREQUAL "SHARED_LIBRARY" AND NOT IRA_BLAS_HIDDEN)
    message(FATAL_ERROR
      "IRA error: linker rejected --version-script and --exclude-libs. Bundled dgemm_ would stay exported. Set IRA_HIDE_BUNDLED_BLAS=OFF to keep that export.")
  endif()
  if(IRA_HIDE_WITH_EXCLUDE_LIBS)
    message(WARNING
      "IRA warning: linker rejected --version-script. Hiding bundled reference BLAS with --exclude-libs instead of exporting dgemm_.")
    add_library(lap_local STATIC "${CMAKE_CURRENT_SOURCE_DIR}/lap_local/lap.f")
    target_link_options(${PROJECT_NAME} PRIVATE "-Wl,--exclude-libs,liblap_local.a")
  else()
    add_library(lap_local OBJECT "${CMAKE_CURRENT_SOURCE_DIR}/lap_local/lap.f")
  endif()
  target_compile_options(lap_local PRIVATE -O0)
  set_target_properties(lap_local PROPERTIES Fortran_FORMAT FIXED POSITION_INDEPENDENT_CODE ON)
  target_link_libraries( ${PROJECT_NAME} lap_local)
  # add lap_local target to the exported targets
  list(APPEND alltargets "lap_local" )
endif()
if(IRA_BLAS_HIDDEN)
  add_custom_command(
    TARGET ${PROJECT_NAME} POST_BUILD
    COMMAND /bin/sh "${CMAKE_CURRENT_SOURCE_DIR}/tests/link_external_blas.sh"
            $<TARGET_FILE:${PROJECT_NAME}>
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

