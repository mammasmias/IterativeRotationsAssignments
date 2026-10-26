enable_testing()


# basic program with f90 caller
add_executable( f90main main.f90 )
set_target_properties( f90main PROPERTIES
  Fortran_PREPROCESS ON
  Fortran_FORMAT FREE )
target_compile_options( f90main PRIVATE
  -ffree-line-length-0)
target_link_libraries( f90main PRIVATE ira )
add_test( NAME f90main COMMAND f90main )

# basic test ira match
add_executable( test_ira test_ira.f90 )
set_target_properties( test_ira PROPERTIES
  Fortran_PREPROCESS ON
  Fortran_FORMAT FREE )
target_compile_options( test_ira PRIVATE
  -ffree-line-length-0)
target_link_libraries( test_ira PRIVATE ira )
add_test( NAME test_ira COMMAND test_ira )

# basic program with C caller
add_executable( cmain main.c )
target_link_libraries( cmain PRIVATE ira )
add_test( NAME cmain COMMAND cmain )


enable_language(CXX)
add_executable( cppmain main.cpp )
target_link_libraries( cppmain PRIVATE ira )
add_test( NAME cppmain COMMAND cppmain )

# The shared library no longer exports the bundled reference BLAS.
# link_external_blas.c is a caller that must name a BLAS itself.
if(IRA_VERSION_SCRIPT_APPLIED)
  add_test(
    NAME link_external_blas
    COMMAND /bin/sh "${CMAKE_CURRENT_SOURCE_DIR}/link_external_blas.sh" $<TARGET_FILE:ira>
    )
endif()
