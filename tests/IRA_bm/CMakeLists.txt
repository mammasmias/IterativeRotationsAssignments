add_executable( ira_bm src/main.f90 src/routines.f90 )
set_target_properties( ira_bm PROPERTIES
  Fortran_FORMAT FREE
  Fortran_PREPROCESS ON )
# ira_test_dbg provides dbg.mod and dbg.h (defined in tests/CMakeLists.txt)
target_link_libraries( ira_bm PRIVATE ira ira_test_dbg )
target_compile_options( ira_bm
  PRIVATE
  $<$<Fortran_COMPILER_ID:GNU>:-ffree-line-length-0>
)

# data is not in git (Zenodo); point to it explicitly
set( IRA_BM_DATA_DIR "${PROJECT_SOURCE_DIR}/benchmark_test/ira_bm/data"
     CACHE PATH "Root of IRA benchmark datasets (<dataset>/xyz/*.xyz)" )
option( IRA_BM_INCLUDE_LARGE "Also run large datasets (lj_clusters)" OFF )

if( NOT IS_DIRECTORY "${IRA_BM_DATA_DIR}" )
  message( WARNING "IRA_BM_DATA_DIR='${IRA_BM_DATA_DIR}' not found; no ira_bm tests added. Download the structures and specify path manually." )
  return()
endif()

# every <dataset>/xyz directory becomes one test
file( GLOB _bm_xyz_dirs LIST_DIRECTORIES true "${IRA_BM_DATA_DIR}/*/xyz" )

foreach( _xyz_dir IN LISTS _bm_xyz_dirs )
  get_filename_component( _dataset_dir "${_xyz_dir}" DIRECTORY )
  get_filename_component( _dataset "${_dataset_dir}" NAME )

  if( _dataset STREQUAL "lj_clusters" AND NOT IRA_BM_INCLUDE_LARGE )
    continue()
  endif()

  # absolute paths (like `realpath`), re-glob when files are added
  file( GLOB _files CONFIGURE_DEPENDS "${_xyz_dir}/*.xyz" )
  if( NOT _files )
    continue()
  endif()
  list( SORT _files COMPARE NATURAL )   # equivalent of `sort -V`

  list( LENGTH _files _n )
  message( STATUS "ira_bm: ${_dataset} (${_n} files)" )

  # write the file list to the build dir, keeps the test command short
  set( _listfile "${CMAKE_CURRENT_BINARY_DIR}/${_dataset}.list" )
  list( JOIN _files "\n" _content )
  file( WRITE "${_listfile}" "${_content}\n" )

  add_test( NAME ira_bm_${_dataset}
            # COMMAND ira_bm --log-file --log-each --list "${_listfile}" )
            # COMMAND ira_bm --log-file --list "${_listfile}" )
            COMMAND ira_bm --list "${_listfile}" )
  set_tests_properties( ira_bm_${_dataset} PROPERTIES
            LABELS "benchmark"
            TIMEOUT 600 )
endforeach()
