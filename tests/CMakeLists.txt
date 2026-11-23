enable_testing()

# shared debug helpers: dbg module (dbg.mod) + ASSERT macro (dbg.h)
# compiled once, so targets do not race writing the same dbg.mod
add_library( ira_test_dbg STATIC dbg.f90 )
set_target_properties( ira_test_dbg PROPERTIES
  Fortran_PREPROCESS ON
  Fortran_FORMAT FREE
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/dbg_mod" )
target_include_directories( ira_test_dbg PUBLIC
  "${CMAKE_CURRENT_SOURCE_DIR}"            # dbg.h
  "${CMAKE_CURRENT_BINARY_DIR}/dbg_mod" )  # dbg.mod
target_link_libraries( ira_test_dbg PUBLIC ira )

# basic program with f90 caller
add_executable( f90main main.f90 )
set_target_properties( f90main PROPERTIES
  Fortran_PREPROCESS ON
  Fortran_FORMAT FREE )
target_compile_options( f90main PRIVATE
  -ffree-line-length-0)
target_link_libraries( f90main PRIVATE ira ira_test_dbg )
add_test( NAME f90main COMMAND f90main )

# basic test ira match
add_executable( test_ira test_ira.f90 )
set_target_properties( test_ira PROPERTIES
  Fortran_PREPROCESS ON
  Fortran_FORMAT FREE )
target_compile_options( test_ira PRIVATE
  -ffree-line-length-0)
target_link_libraries( test_ira PRIVATE ira ira_test_dbg )
add_test( NAME test_ira COMMAND test_ira )

# basic program with C caller
add_executable( cmain main.c )
target_link_libraries( cmain PRIVATE ira )
add_test( NAME cmain COMMAND cmain )


enable_language(CXX)
add_executable( cppmain main.cpp )
target_link_libraries( cppmain PRIVATE ira )
add_test( NAME cppmain COMMAND cppmain )

# libira defines ira_dgemm_, not dgemm_. This caller brings its own BLAS.
get_target_property(ira_type ira TYPE)
if(ira_type STREQUAL "SHARED_LIBRARY")
  add_test(
    NAME link_external_blas
    COMMAND /bin/sh "${CMAKE_CURRENT_SOURCE_DIR}/link_external_blas.sh" $<TARGET_FILE:ira>
    )
endif()

if(IRA_BUILD_BENCHMARK)
  add_subdirectory( IRA_bm )
endif()
